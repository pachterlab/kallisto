/*
 * Copyright (c) 2020-2025, NVIDIA CORPORATION.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <test_utils.hpp>

#include <cuco/static_map.cuh>

#include <cuda/functional>
#include <cuda/std/tuple>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>

#include <limits>

template <std::size_t ValidSize, typename Ref>
__global__ void shared_memory_test_kernel(Ref* maps,
                                          typename Ref::key_type const* const insterted_keys,
                                          typename Ref::mapped_type const* const inserted_values,
                                          size_t number_of_elements,
                                          bool* const keys_exist,
                                          bool* const keys_and_values_correct)
{
  // Each block processes one map
  const size_t map_id = blockIdx.x;
  const size_t offset = map_id * number_of_elements;

  __shared__ typename Ref::value_type sm_buffer[ValidSize];

  auto g          = cuco::test::cg::this_thread_block();
  auto insert_ref = maps[map_id].make_copy(g, sm_buffer, cuco::thread_scope_block);
  auto find_ref   = insert_ref.rebind_operators(cuco::op::find);

  for (int i = g.thread_rank(); i < number_of_elements; i += g.size()) {
    auto found_pair_it = find_ref.find(insterted_keys[offset + i]);

    if (found_pair_it != find_ref.end()) {
      keys_exist[offset + i] = true;
      if (found_pair_it->first == insterted_keys[offset + i] and
          found_pair_it->second == inserted_values[offset + i]) {
        keys_and_values_correct[offset + i] = true;
      } else {
        keys_and_values_correct[offset + i] = false;
      }
    } else {
      keys_exist[offset + i]              = false;
      keys_and_values_correct[offset + i] = true;
    }
  }
}

TEMPLATE_TEST_CASE_SIG("static_map shared memory tests",
                       "",
                       ((typename Key, typename Value), Key, Value),
                       (int32_t, int32_t),
                       (int32_t, int64_t),
                       (int64_t, int32_t),
                       (int64_t, int64_t))
{
  constexpr std::size_t number_of_maps  = 1000;
  constexpr std::size_t elements_in_map = 500;
  constexpr std::size_t map_capacity    = 2 * elements_in_map;

  using extent_type = cuco::extent<std::size_t, map_capacity>;
  using map_type    = cuco::static_map<Key,
                                       Value,
                                       extent_type,
                                       cuda::thread_scope_device,
                                       cuda::std::equal_to<Key>,
                                       cuco::linear_probing<1, cuco::default_hash_function<Key>>>;

  // one array for all maps, first elements_in_map element belong to map 0, second to map 1 and so
  // on
  thrust::device_vector<Key> d_keys(number_of_maps * elements_in_map);
  thrust::device_vector<Value> d_values(number_of_maps * elements_in_map);

  thrust::sequence(thrust::device, d_keys.begin(), d_keys.end());
  thrust::sequence(thrust::device, d_values.begin(), d_values.end(), 1);

  // using std::unique_ptr because static_map does not have copy/move constructor/assignment
  // operator yet
  std::vector<std::unique_ptr<map_type>> maps;
  for (std::size_t map_id = 0; map_id < number_of_maps; ++map_id) {
    maps.push_back(std::make_unique<map_type>(
      extent_type{}, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{-1}));
  }

  thrust::device_vector<bool> d_keys_exist(number_of_maps * elements_in_map);
  thrust::device_vector<bool> d_keys_and_values_correct(number_of_maps * elements_in_map);

  using ref_type = typename map_type::template ref_type<cuco::op::insert_tag>;

  SECTION("Keys are all found after insertion.")
  {
    auto pairs_begin = thrust::make_transform_iterator(
      thrust::counting_iterator{0},
      cuda::proclaim_return_type<cuco::pair<Key, Value>>(
        [d_keys = d_keys.data(), d_values = d_values.data()] __device__(int idx) {
          return cuco::pair<Key, Value>(d_keys[idx], d_values[idx]);
        }));
    std::vector<ref_type> h_refs;
    for (std::size_t map_id = 0; map_id < number_of_maps; ++map_id) {
      const std::size_t offset = map_id * elements_in_map;

      map_type* map = maps[map_id].get();
      map->insert(pairs_begin + offset, pairs_begin + offset + elements_in_map);
      h_refs.push_back(map->ref(cuco::op::insert));
    }
    thrust::device_vector<ref_type> d_refs(h_refs);

    // maybe_unused to silence false positive "variable set but not used" warning
    [[maybe_unused]] auto constexpr valid_size =
      cuco::make_valid_extent<typename ref_type::probing_scheme_type,
                              typename ref_type::storage_ref_type>(extent_type{});

    shared_memory_test_kernel<valid_size.value(), ref_type>
      <<<number_of_maps, 64>>>(d_refs.data().get(),
                               d_keys.data().get(),
                               d_values.data().get(),
                               elements_in_map,
                               d_keys_exist.data().get(),
                               d_keys_and_values_correct.data().get());

    REQUIRE(d_keys_exist.size() == d_keys_and_values_correct.size());
    auto zip = thrust::make_zip_iterator(
      cuda::std::tuple{d_keys_exist.begin(), d_keys_and_values_correct.begin()});

    REQUIRE(cuco::test::all_of(zip,
                               zip + d_keys_exist.size(),
                               cuda::proclaim_return_type<bool>([] __device__(auto const& z) {
                                 return cuda::std::get<0>(z) and cuda::std::get<1>(z);
                               })));
  }

  SECTION("No key is found before insertion.")
  {
    std::vector<ref_type> h_refs;
    for (std::size_t map_id = 0; map_id < number_of_maps; ++map_id) {
      h_refs.push_back(maps[map_id].get()->ref(cuco::op::insert));
    }
    thrust::device_vector<ref_type> d_refs(h_refs);

    // maybe_unused to silence false positive "variable set but not used" warning
    [[maybe_unused]] auto constexpr valid_size =
      cuco::make_valid_extent<typename ref_type::probing_scheme_type,
                              typename ref_type::storage_ref_type>(extent_type{});

    shared_memory_test_kernel<valid_size.value(), ref_type>
      <<<number_of_maps, 64>>>(d_refs.data().get(),
                               d_keys.data().get(),
                               d_values.data().get(),
                               elements_in_map,
                               d_keys_exist.data().get(),
                               d_keys_and_values_correct.data().get());

    REQUIRE(cuco::test::none_of(d_keys_exist.begin(), d_keys_exist.end(), cuda::std::identity{}));
  }
}

auto constexpr cg_size     = 1;
auto constexpr bucket_size = 1;

template <std::size_t ValidSize>
__global__ void shared_memory_hash_table_kernel(bool* key_found)
{
  using Key       = int32_t;
  using Value     = int32_t;
  using slot_type = cuco::pair<Key, Value>;

  __shared__ slot_type map[ValidSize];

  using extent_type      = cuco::extent<std::size_t, ValidSize>;
  using storage_ref_type = cuco::bucket_storage_ref<slot_type, bucket_size, extent_type>;

  auto raw_ref =
    cuco::static_map_ref{cuco::empty_key<Key>{-1},
                         cuco::empty_value<Value>{-1},
                         cuda::std::equal_to<Key>{},
                         cuco::linear_probing<cg_size, cuco::default_hash_function<Key>>{},
                         cuco::thread_scope_block,
                         storage_ref_type{extent_type{}, map}};

  auto const block = cooperative_groups::this_thread_block();
  raw_ref.initialize(block);

  auto const index = threadIdx.x + blockIdx.x * blockDim.x;
  auto const rank  = block.thread_rank();

  // insert {thread_rank, thread_rank} for each thread in thread-block
  auto insert_ref = raw_ref.rebind_operators(cuco::op::insert);
  insert_ref.insert(slot_type{rank, rank});
  block.sync();

  auto find_ref             = insert_ref.rebind_operators(cuco::op::find);
  auto const retrieved_pair = find_ref.find(rank);
  block.sync();

  if (retrieved_pair != find_ref.end() && retrieved_pair->second == rank) {
    key_found[index] = true;
  }
}

TEST_CASE("static map shared memory slots.", "")
{
  constexpr std::size_t N = 256;
  // maybe_unused to silence false positive "variable set but not used" warning
  [[maybe_unused]] auto constexpr valid_size =
    cuco::make_valid_extent<cg_size, bucket_size>(cuco::extent<std::size_t, N>{});

  thrust::device_vector<bool> key_found(N, false);
  shared_memory_hash_table_kernel<valid_size.value()><<<8, 32>>>(key_found.data().get());
  CUCO_CUDA_TRY(cudaDeviceSynchronize());

  REQUIRE(cuco::test::all_of(key_found.begin(), key_found.end(), cuda::std::identity{}));
}
