/*
 * Copyright (c) 2024, NVIDIA CORPORATION.
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
#include <cuda/std/iterator>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/sort.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = int32_t;

int32_t constexpr SENTINEL = -1;

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

  thrust::device_vector<cuco::pair<Key, Value>> d_results(num_keys);
  auto output_begin = d_results.begin();

  SECTION("Non-inserted keys have empty retrieval output")
  {
    auto const count = map.count(keys_begin, keys_begin + num_keys);
    thrust::device_vector<cuco::pair<Key, Value>> d_results(num_keys);

    REQUIRE(count == 0);

    auto const [_, output_end] = map.retrieve(
      keys_begin, keys_begin + num_keys, thrust::discard_iterator{}, d_results.begin());
    auto const size = std::distance(d_results.begin(), output_end);

    REQUIRE(size == 0);
  }

  map.insert(pairs_begin, pairs_begin + num_keys);

  SECTION("Total count should be equal to the number of inserted pairs.")
  {
    // Count matching keys
    auto const count = map.count(keys_begin, keys_begin + num_keys);

    REQUIRE(count == num_keys);

    auto [_, output_end] =
      map.retrieve(keys_begin, keys_begin + num_keys, thrust::discard_iterator{}, output_begin);
    auto const size = cuda::std::distance(output_begin, output_end);

    REQUIRE(size == num_keys);

    // sort before compare
    thrust::sort(
      thrust::device,
      d_results.begin(),
      d_results.end(),
      [] __device__(const cuco::pair<Key, Value>& lhs, const cuco::pair<Key, Value>& rhs) {
        return lhs.first < rhs.first;
      });

    REQUIRE(
      cuco::test::equal(pairs_begin,
                        pairs_begin + num_keys,
                        output_begin,
                        [] __device__(cuco::pair<Key, Value> lhs, cuco::pair<Key, Value> rhs) {
                          return lhs.first == rhs.first and lhs.second == rhs.second;
                        }));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_map: retrieve tests",
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
  constexpr size_type num_keys{1'000};

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
                              cuco::storage<1>>{
    extent_type{}, cuco::empty_key<Key>{SENTINEL}, cuco::empty_value<Value>{SENTINEL}};

  test_unique_sequence(map, num_keys);
}
