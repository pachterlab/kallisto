/*
 * Copyright (c) 2024-2025, NVIDIA CORPORATION.
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

#include <cuco/static_multiset.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = int32_t;

template <typename Set>
void test_unique_sequence(Set& set, size_type num_keys)
{
  using Key = typename Set::key_type;

  auto keys_begin = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>{0},
    cuda::proclaim_return_type<Key>([] __device__(auto i) { return Key{i}; }));

  SECTION("Count of empty set should be zero.")
  {
    auto const count = set.count(keys_begin, keys_begin + num_keys);
    REQUIRE(count == 0);
  }

  set.insert(keys_begin, keys_begin + num_keys);

  SECTION("Count of n unique keys should be n.")
  {
    auto const count = set.count(keys_begin, keys_begin + num_keys);
    REQUIRE(count == num_keys);
  }

  auto constexpr multiplicity = 3;
  auto query_begin            = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>{0},
    cuda::proclaim_return_type<Key>([] __device__(auto i) { return Key{i / multiplicity}; }));

  SECTION("Count of 3n unique keys should be 3n.")
  {
    auto const count = set.count(query_begin, query_begin + num_keys * multiplicity);
    REQUIRE(count == num_keys * multiplicity);
  }
}

template <typename Set>
void test_count_each(Set& set, size_type num_keys)
{
  using Key = typename Set::key_type;

  thrust::device_vector<size_type> d_counts(num_keys);
  auto const counts_begin = d_counts.begin();

  auto keys_begin = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>{0},
    cuda::proclaim_return_type<Key>([] __device__(auto i) { return Key{i}; }));

  set.clear();

  SECTION("Count_each of empty set should be all zeros.")
  {
    set.count_each(
      keys_begin, keys_begin + num_keys, set.key_eq(), set.hash_function(), counts_begin);
    REQUIRE(cuco::test::all_of(
      d_counts.begin(),
      d_counts.end(),
      cuda::proclaim_return_type<bool>([] __device__(size_type count) { return count == 0; })));
  }

  set.insert(keys_begin, keys_begin + num_keys);

  SECTION("Count_each of n unique keys should be all ones.")
  {
    set.count_each(
      keys_begin, keys_begin + num_keys, set.key_eq(), set.hash_function(), counts_begin);
    REQUIRE(cuco::test::all_of(
      d_counts.begin(),
      d_counts.end(),
      cuda::proclaim_return_type<bool>([] __device__(size_type count) { return count == 1; })));
  }

  set.clear();

  auto constexpr multiplicity = 3;
  auto duplicate_keys_begin   = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>{0},
    cuda::proclaim_return_type<Key>([] __device__(auto i) { return Key{i / multiplicity}; }));
  set.insert(duplicate_keys_begin, duplicate_keys_begin + num_keys);

  auto const query_begin = thrust::counting_iterator<size_type>{0};
  auto const query_size  = num_keys / multiplicity;
  SECTION("Count_each with duplicates should return correct counts.")
  {
    set.count_each(
      query_begin, query_begin + query_size, set.key_eq(), set.hash_function(), counts_begin);
    REQUIRE(cuco::test::all_of(d_counts.begin(),
                               d_counts.begin() + query_size,
                               cuda::proclaim_return_type<bool>([] __device__(size_type count) {
                                 return count == multiplicity;
                               })));
  }
}

template <typename Set>
void test_count_each_outer(Set& set, size_type num_keys)
{
  using Key = typename Set::key_type;

  thrust::device_vector<size_type> d_counts(num_keys);
  auto const counts_begin = d_counts.begin();

  auto keys_begin = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>{0},
    cuda::proclaim_return_type<Key>([] __device__(auto i) { return Key{i}; }));

  set.clear();

  SECTION("Count_each_outer of empty set should be all ones.")
  {
    set.count_each_outer(
      keys_begin, keys_begin + num_keys, set.key_eq(), set.hash_function(), counts_begin);
    REQUIRE(cuco::test::all_of(
      d_counts.begin(),
      d_counts.end(),
      cuda::proclaim_return_type<bool>([] __device__(size_type count) { return count == 1; })));
  }

  set.insert(keys_begin, keys_begin + num_keys);

  SECTION("Count_each_outer of n unique keys should be all ones.")
  {
    set.count_each_outer(
      keys_begin, keys_begin + num_keys, set.key_eq(), set.hash_function(), counts_begin);
    REQUIRE(cuco::test::all_of(
      d_counts.begin(),
      d_counts.end(),
      cuda::proclaim_return_type<bool>([] __device__(size_type count) { return count == 1; })));
  }

  set.clear();

  auto constexpr multiplicity = 3;
  auto duplicate_keys_begin   = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>{0},
    cuda::proclaim_return_type<Key>([] __device__(auto i) { return Key{i / multiplicity}; }));
  set.insert(duplicate_keys_begin, duplicate_keys_begin + num_keys);

  auto const query_size  = num_keys / multiplicity;
  auto const query_begin = thrust::counting_iterator<size_type>{0};

  SECTION("Count_each_outer with duplicates should return correct counts.")
  {
    set.count_each_outer(
      query_begin, query_begin + query_size, set.key_eq(), set.hash_function(), counts_begin);
    REQUIRE(cuco::test::all_of(d_counts.begin(),
                               d_counts.begin() + query_size,
                               cuda::proclaim_return_type<bool>([] __device__(size_type count) {
                                 return count == multiplicity;
                               })));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_multiset count tests",
  "",
  ((typename Key, cuco::test::probe_sequence Probe, int CGSize), Key, Probe, CGSize),
  (int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, cuco::test::probe_sequence::linear_probing, 2))
{
  constexpr size_type num_keys{666};

  using probe = std::conditional_t<Probe == cuco::test::probe_sequence::linear_probing,
                                   cuco::linear_probing<CGSize, cuco::default_hash_function<Key>>,
                                   cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>>;

  auto set =
    cuco::static_multiset{num_keys, cuco::empty_key<Key>{-1}, {}, probe{}, {}, cuco::storage<2>{}};

  test_unique_sequence(set, num_keys);
  test_count_each(set, num_keys);
  test_count_each_outer(set, num_keys);
}
