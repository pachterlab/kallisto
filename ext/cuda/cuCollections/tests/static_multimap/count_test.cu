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

#include <cuco/static_multimap.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/generate.h>
#include <thrust/sequence.h>
#include <thrust/shuffle.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = int32_t;

static size_type constexpr multiplicity = 3;

template <typename Map>
void test_multiplicity_count(Map& map, size_type num_keys)
{
  using Key   = typename Map::key_type;
  using Value = typename Map::mapped_type;

  auto const keys_begin = thrust::counting_iterator<Key>{0};

  SECTION("Count of empty map should be zero.")
  {
    auto const count = map.count(keys_begin, keys_begin + num_keys);
    REQUIRE(count == 0);
  }

  SECTION("Count of n unique keys should be n.")
  {
    auto const pairs_begin = thrust::make_transform_iterator(
      thrust::make_counting_iterator<size_type>(0),
      cuda::proclaim_return_type<cuco::pair<Key, Value>>(
        [] __device__(auto i) { return cuco::pair<Key, Value>{i, i}; }));
    map.insert(pairs_begin, pairs_begin + num_keys);

    auto const count = map.count(keys_begin, keys_begin + num_keys);
    REQUIRE(count == num_keys);
  }

  SECTION("Count of n unique keys should be n x multiplicity.")
  {
    auto const pairs_begin = thrust::make_transform_iterator(
      thrust::make_counting_iterator<size_type>(0),
      cuda::proclaim_return_type<cuco::pair<Key, Value>>(
        [] __device__(auto i) { return cuco::pair<Key, Value>{i / multiplicity, i}; }));
    map.insert(pairs_begin, pairs_begin + num_keys * multiplicity);

    auto const count = map.count(keys_begin, keys_begin + num_keys);
    REQUIRE(count == num_keys * multiplicity);
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_multimap count tests",
  "",
  ((typename T, cuco::test::probe_sequence Probe, int CGSize), T, Probe, CGSize),
  (int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, cuco::test::probe_sequence::linear_probing, 2))
{
  constexpr size_type num_keys{1'000};

  using probe = std::conditional_t<
    Probe == cuco::test::probe_sequence::linear_probing,
    cuco::linear_probing<CGSize, cuco::default_hash_function<T>>,
    cuco::double_hashing<CGSize, cuco::default_hash_function<T>, cuco::default_hash_function<T>>>;

  auto map = cuco::experimental::static_multimap<T,
                                                 T,
                                                 cuco::extent<size_type>,
                                                 cuda::thread_scope_device,
                                                 cuda::std::equal_to<T>,
                                                 probe,
                                                 cuco::cuda_allocator<cuda::std::byte>,
                                                 cuco::storage<2>>{
    num_keys * multiplicity, cuco::empty_key<T>{-1}, cuco::empty_value<T>{-1}};

  test_multiplicity_count(map, num_keys);
}
