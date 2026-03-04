/*
 * Copyright (c) 2023-2025, NVIDIA CORPORATION.
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
#include <cuco/utility/reduction_functors.cuh>

#include <cuda/atomic>
#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/zip_iterator.h>

#include <catch2/catch_template_test_macros.hpp>

#include <cstdint>

using size_type = std::size_t;

template <bool HasInit, typename Map, typename Init>
void test_insert_or_apply(Map& map, size_type num_keys, size_type num_unique_keys, Init init)
{
  REQUIRE((num_keys % num_unique_keys) == 0);

  using Key   = typename Map::key_type;
  using Value = typename Map::mapped_type;

  // Insert pairs
  auto pairs_begin = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>(0),
    cuda::proclaim_return_type<cuco::pair<Key, Value>>([num_unique_keys] __device__(auto i) {
      return cuco::pair<Key, Value>{i % num_unique_keys, 1};
    }));

  auto constexpr plus_op = cuco::reduce::plus{};
  if constexpr (HasInit) {
    map.insert_or_apply(pairs_begin, pairs_begin + num_keys, init, plus_op);
  } else {
    map.insert_or_apply(pairs_begin, pairs_begin + num_keys, plus_op);
  }

  REQUIRE(map.size() == num_unique_keys);

  thrust::device_vector<Key> d_keys(num_unique_keys);
  thrust::device_vector<Value> d_values(num_unique_keys);
  map.retrieve_all(d_keys.begin(), d_values.begin());

  REQUIRE(cuco::test::equal(d_values.begin(),
                            d_values.end(),
                            thrust::make_constant_iterator<Value>(num_keys / num_unique_keys),
                            cuda::std::equal_to<Value>{}));
}

template <bool HasInit, typename Map, typename Init>
void test_insert_or_apply_shmem(Map& map, size_type num_keys, size_type num_unique_keys, Init init)
{
  REQUIRE((num_keys % num_unique_keys) == 0);

  using Key   = typename Map::key_type;
  using Value = typename Map::mapped_type;

  using KeyEqual         = typename Map::key_equal;
  using ProbingScheme    = typename Map::probing_scheme_type;
  using Allocator        = typename Map::allocator_type;
  auto constexpr cg_size = Map::cg_size;

  int32_t constexpr shmem_block_size = 1024;

  using shmem_size_type = int32_t;

  shmem_size_type constexpr cardinality_threshold   = shmem_block_size;
  shmem_size_type constexpr shared_map_num_elements = cardinality_threshold + shmem_block_size;
  float constexpr load_factor                       = 0.7;
  shmem_size_type constexpr shared_map_size =
    static_cast<shmem_size_type>((1.0 / load_factor) * shared_map_num_elements);

  using extent_type     = cuco::extent<shmem_size_type, shared_map_size>;
  using shared_map_type = cuco::static_map<Key,
                                           Value,
                                           extent_type,
                                           cuda::thread_scope_block,
                                           KeyEqual,
                                           ProbingScheme,
                                           Allocator,
                                           cuco::storage<1>>;

  using shared_map_ref_type = typename shared_map_type::template ref_type<>;
  auto constexpr valid_extent =
    cuco::make_valid_extent<typename shared_map_ref_type::probing_scheme_type,
                            typename shared_map_ref_type::storage_ref_type>(extent_type{});

  // Insert pairs
  auto pairs_begin = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>(0),
    cuda::proclaim_return_type<cuco::pair<Key, Value>>([num_unique_keys] __device__(auto i) {
      return cuco::pair<Key, Value>{i % num_unique_keys, 1};
    }));

  auto const shmem_grid_size = cuco::detail::grid_size(num_keys, cg_size, 1, shmem_block_size);

  cuda::stream_ref stream{cudaStream_t{nullptr}};

  // launch the shmem kernel
  cuco::detail::static_map_ns::
    insert_or_apply_shmem<HasInit, cg_size, shmem_block_size, shared_map_ref_type>
    <<<shmem_grid_size, shmem_block_size, 0, stream.get()>>>(pairs_begin,
                                                             num_keys,
                                                             init,
                                                             cuco::reduce::plus{},
                                                             map.ref(cuco::op::insert_or_apply),
                                                             valid_extent);

  REQUIRE(map.size() == num_unique_keys);

  thrust::device_vector<Key> d_keys(num_unique_keys);
  thrust::device_vector<Value> d_values(num_unique_keys);
  map.retrieve_all(d_keys.begin(), d_values.begin());

  REQUIRE(cuco::test::equal(d_values.begin(),
                            d_values.end(),
                            thrust::make_constant_iterator<Value>(num_keys / num_unique_keys),
                            cuda::std::equal_to<Value>{}));
}

/*
TEMPLATE_TEST_CASE_SIG(
  "static_map insert_or_apply tests",
  "",
  ((typename Key, typename Value, cuco::test::probe_sequence Probe, int CGSize),
   Key,
   Value,
   Probe,
   CGSize),
  (int32_t, int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int32_t, int64_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, int64_t, cuco::test::probe_sequence::linear_probing, 2))
{
  constexpr size_type num_keys{10'000};
  constexpr size_type num_unique_keys{100};

  using probe = std::conditional_t<
    Probe == cuco::test::probe_sequence::linear_probing,
    cuco::linear_probing<CGSize, cuco::murmurhash3_32<Key>>,
    cuco::double_hashing<CGSize, cuco::murmurhash3_32<Key>, cuco::murmurhash3_32<Key>>>;

  using map_type = cuco::static_map<Key,
                                    Value,
                                    cuco::extent<size_type>,
                                    cuda::thread_scope_device,
                                    cuda::std::equal_to<Key>,
                                    probe,
                                    cuco::cuda_allocator<cuda::std::byte>,
                                    cuco::storage<2>>;

  SECTION("sentinel equals init; has_init = true")
  {
    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply<true>(map, num_keys, num_unique_keys, static_cast<Value>(0));
  }
  SECTION("sentinel equals init; has_init = false")
  {
    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply<false>(map, num_keys, num_unique_keys, static_cast<Value>(0));
  }
  SECTION("sentinel not equals init; has_init = true")
  {
    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply<true>(map, num_keys, num_unique_keys, static_cast<Value>(-1));
  }
  SECTION("sentinel not equals init; has_init = false")
  {
    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply<false>(map, num_keys, num_unique_keys, static_cast<Value>(-1));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_map insert_or_apply all unique keys tests", "", ((typename Key)), (int32_t), (int64_t))
{
  using Value = Key;

  constexpr size_type num_keys = 100;

  using map_type = cuco::static_map<Key,
                                    Value,
                                    cuco::extent<size_type>,
                                    cuda::thread_scope_device,
                                    cuda::std::equal_to<Key>,
                                    cuco::linear_probing<2, cuco::murmurhash3_32<Key>>,
                                    cuco::cuda_allocator<cuda::std::byte>,
                                    cuco::storage<2>>;

  SECTION("sentinel equals init; has_init = true")
  {
    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply<true>(map, num_keys, num_keys, static_cast<Value>(0));
  }
  SECTION("sentinel equals init; has_init = false")
  {
    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply<false>(map, num_keys, num_keys, static_cast<Value>(0));
  }
  SECTION("sentinel not equals init; has_init = true")
  {
    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply<true>(map, num_keys, num_keys, static_cast<Value>(-1));
  }
  SECTION("sentinel not equals init; has_init = false")
  {
    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply<false>(map, num_keys, num_keys, static_cast<Value>(-1));
  }
}
*/

TEMPLATE_TEST_CASE_SIG(
  "static_map insert_or_apply shared memory", "", ((typename Key)), (int32_t), (int64_t))
{
  using Value = Key;

  using map_type = cuco::static_map<Key,
                                    Value,
                                    cuco::extent<size_type>,
                                    cuda::thread_scope_device,
                                    cuda::std::equal_to<Key>,
                                    cuco::linear_probing<1, cuco::murmurhash3_32<Key>>,
                                    cuco::cuda_allocator<cuda::std::byte>,
                                    cuco::storage<2>>;

  SECTION("duplicate keys")
  {
    constexpr size_type num_keys        = 10'000;
    constexpr size_type num_unique_keys = 100;

    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply_shmem<true>(map, num_keys, num_unique_keys, static_cast<Value>(0));
  }

  SECTION("unique keys")
  {
    constexpr size_type num_keys        = 10'000;
    constexpr size_type num_unique_keys = num_keys;

    auto map = map_type{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{0}};
    test_insert_or_apply_shmem<true>(map, num_keys, num_unique_keys, static_cast<Value>(0));
  }
}
