/*
 * Copyright (c) 2025, NVIDIA CORPORATION.
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

#include <cuco/static_multimap.cuh>

#include <cuda/std/functional>
#include <cuda/std/iterator>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/sort.h>

#include <catch2/catch_template_test_macros.hpp>

struct custom_key_eq {
  template <typename T>
  __device__ bool operator()(T const& lhs, T const& rhs) const
  {
    return lhs % 2 == 0 ? lhs == rhs : false;
  }
};

template <typename Map>
void test_retrieve(Map& map, std::size_t num_items)
{
  using Key   = typename Map::key_type;
  using Value = typename Map::mapped_type;

  auto const num_gold = num_items / 2;

  auto const keys_begin = thrust::counting_iterator<Key>{0};
  // multiplicity = 2
  auto const pairs_begin = thrust::make_transform_iterator(
    keys_begin, cuda::proclaim_return_type<cuco::pair<Key, Value>>([] __device__(auto i) {
      return cuco::pair<Key, Value>{i / 2, i / 2};
    }));

  thrust::device_vector<cuco::pair<Key, Value>> d_results(num_gold);
  auto output_begin = d_results.begin();

  map.insert(pairs_begin, pairs_begin + num_items);

  SECTION("Total count should be equal to the number of inserted pairs.")
  {
    // Count matching keys
    auto const num =
      map.count(keys_begin, keys_begin + num_items, custom_key_eq{}, map.hash_function());

    REQUIRE(num == num_gold);

    auto [_, output_end]   = map.retrieve(keys_begin,
                                        keys_begin + num_items,
                                        custom_key_eq{},
                                        map.hash_function(),
                                        thrust::discard_iterator{},
                                        output_begin);
    std::size_t const size = cuda::std::distance(output_begin, output_end);

    REQUIRE(size == num_gold);

    // sort before compare
    thrust::sort(
      thrust::device,
      d_results.begin(),
      d_results.end(),
      [] __device__(const cuco::pair<Key, Value>& lhs, const cuco::pair<Key, Value>& rhs) {
        if (lhs.first != rhs.first) { return lhs.first < rhs.first; }
        return lhs.second < rhs.second;
      });

    auto const gold_begin = thrust::make_transform_iterator(
      keys_begin, cuda::proclaim_return_type<cuco::pair<Key, Value>>([] __device__(auto i) {
        return cuco::pair<Key, Value>{(i / 2) * 2, (i / 2) * 2};
      }));
    REQUIRE(
      cuco::test::equal(gold_begin,
                        gold_begin + num_gold,
                        output_begin,
                        [] __device__(cuco::pair<Key, Value> lhs, cuco::pair<Key, Value> rhs) {
                          return lhs.first == rhs.first and lhs.second == rhs.second;
                        }));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_multimap retrieve tests",
  "",
  ((typename T, cuco::test::probe_sequence Probe, int CGSize), T, Probe, CGSize),
  (int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, cuco::test::probe_sequence::double_hashing, 4),
  (int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, cuco::test::probe_sequence::double_hashing, 4),
  (int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, cuco::test::probe_sequence::linear_probing, 4),
  (int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, cuco::test::probe_sequence::linear_probing, 4))
{
  constexpr std::size_t num_items{1'000};

  using probe = std::conditional_t<
    Probe == cuco::test::probe_sequence::linear_probing,
    cuco::linear_probing<CGSize, cuco::default_hash_function<T>>,
    cuco::double_hashing<CGSize, cuco::default_hash_function<T>, cuco::default_hash_function<T>>>;

  auto map = cuco::experimental::static_multimap<T,
                                                 T,
                                                 cuco::extent<std::size_t>,
                                                 cuda::thread_scope_device,
                                                 cuda::std::equal_to<T>,
                                                 probe,
                                                 cuco::cuda_allocator<cuda::std::byte>,
                                                 cuco::storage<2>>{
    num_items * 2, cuco::empty_key<T>{-1}, cuco::empty_value<T>{-1}};

  test_retrieve(map, num_items);
}
