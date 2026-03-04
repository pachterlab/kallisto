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
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sort.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = int32_t;

template <typename Map>
void test_unique_sequence(Map& map, size_type num_keys)
{
  using Key   = typename Map::key_type;
  using Value = typename Map::mapped_type;

  auto keys_begin  = thrust::counting_iterator<Key>{0};
  auto pairs_begin = thrust::make_transform_iterator(
    thrust::make_counting_iterator<size_type>(0),
    cuda::proclaim_return_type<cuco::pair<Key, Value>>(
      [] __device__(auto i) { return cuco::pair<Key, Value>{i, i}; }));
  thrust::device_vector<bool> d_contained(num_keys);

  auto zip_equal = cuda::proclaim_return_type<bool>(
    [] __device__(auto const& p) { return cuda::std::get<0>(p) == cuda::std::get<1>(p); });
  auto is_even =
    cuda::proclaim_return_type<bool>([] __device__(auto const& i) { return i % 2 == 0; });

  SECTION("Non-inserted keys should not be contained.")
  {
    REQUIRE(map.size() == 0);

    REQUIRE(map.count(keys_begin, keys_begin + num_keys) == 0);

    map.contains(keys_begin, keys_begin + num_keys, d_contained.begin());
    REQUIRE(cuco::test::none_of(d_contained.begin(), d_contained.end(), cuda::std::identity{}));
  }

  SECTION("All conditionally inserted keys should be contained")
  {
    auto const inserted = map.insert_if(
      pairs_begin, pairs_begin + num_keys, thrust::counting_iterator<std::size_t>(0), is_even);
    REQUIRE(inserted == num_keys / 2);
    REQUIRE(map.size() == num_keys / 2);

    map.contains(keys_begin, keys_begin + num_keys, d_contained.begin());
    REQUIRE(cuco::test::equal(
      d_contained.begin(),
      d_contained.end(),
      thrust::counting_iterator<std::size_t>(0),
      cuda::proclaim_return_type<bool>([] __device__(auto const& idx_contained, auto const& idx) {
        return ((idx % 2) == 0) == idx_contained;
      })));
  }

  map.insert(pairs_begin, pairs_begin + num_keys);

  SECTION("All inserted keys should be contained.")
  {
    REQUIRE(map.count(keys_begin, keys_begin + num_keys) == num_keys);

    map.contains(keys_begin, keys_begin + num_keys, d_contained.begin());
    REQUIRE(cuco::test::all_of(d_contained.begin(), d_contained.end(), cuda::std::identity{}));
  }

  SECTION("Conditional contains should return true on even inputs.")
  {
    map.contains_if(keys_begin,
                    keys_begin + num_keys,
                    thrust::counting_iterator<std::size_t>(0),
                    is_even,
                    d_contained.begin());
    auto gold_iter =
      thrust::make_transform_iterator(thrust::counting_iterator<std::size_t>(0), is_even);
    auto zip = thrust::make_zip_iterator(cuda::std::tuple{d_contained.begin(), gold_iter});
    REQUIRE(cuco::test::all_of(zip, zip + num_keys, zip_equal));
  }

  SECTION("All inserted key-values should be properly retrieved")
  {
    thrust::device_vector<Value> d_values(num_keys);

    auto const [keys_end, values_end] = map.retrieve_all(keys_begin, d_values.begin());
    REQUIRE(std::distance(keys_begin, keys_end) == num_keys);
    REQUIRE(std::distance(d_values.begin(), values_end) == num_keys);

    thrust::sort(thrust::device, d_values.begin(), values_end);
    REQUIRE(cuco::test::equal(d_values.begin(),
                              values_end,
                              thrust::make_counting_iterator<Value>(0),
                              cuda::std::equal_to<Value>{}));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_map: contains + retrieve_all tests",
  "",
  ((typename Key, typename Value, cuco::test::probe_sequence Probe, int CGSize),
   Key,
   Value,
   Probe,
   CGSize),
  (int32_t, int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, int64_t, cuco::test::probe_sequence::linear_probing, 2))
{
  constexpr size_type num_keys{400};

  // XXX: testing static extent is intended, DO NOT CHANGE
  using extent_type = cuco::extent<size_type, num_keys>;
  using probe       = std::conditional_t<
          Probe == cuco::test::probe_sequence::linear_probing,
          cuco::linear_probing<CGSize, cuco::murmurhash3_32<Key>>,
          cuco::double_hashing<CGSize, cuco::murmurhash3_32<Key>, cuco::murmurhash3_32<Key>>>;

  auto map = cuco::static_map<Key,
                              Value,
                              extent_type,
                              cuda::thread_scope_device,
                              cuda::std::equal_to<Key>,
                              probe,
                              cuco::cuda_allocator<cuda::std::byte>,
                              cuco::storage<2>>{
    extent_type{}, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{-1}};

  test_unique_sequence(map, num_keys);
}
