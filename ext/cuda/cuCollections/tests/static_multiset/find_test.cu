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
#include <cuda/std/tuple>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = int32_t;

int32_t constexpr SENTINEL = -1;

template <typename Set>
void test_unique_sequence(Set& set, size_type num_keys)
{
  using Key = typename Set::key_type;

  auto keys_begin = thrust::counting_iterator<Key>{0};
  thrust::device_vector<bool> d_contained(num_keys);

  auto zip_equal = cuda::proclaim_return_type<bool>(
    [] __device__(auto const& p) { return cuda::std::get<0>(p) == cuda::std::get<1>(p); });

  SECTION("Non-inserted keys have no matches")
  {
    thrust::device_vector<Key> d_results(num_keys);

    set.find(keys_begin, keys_begin + num_keys, d_results.begin());
    auto zip = thrust::make_zip_iterator(cuda::std::tuple{
      d_results.begin(), thrust::constant_iterator<Key>{set.empty_key_sentinel()}});

    REQUIRE(cuco::test::all_of(zip, zip + num_keys, zip_equal));
  }

  set.insert(keys_begin, keys_begin + num_keys);

  SECTION("All inserted keys should be correctly recovered during find")
  {
    thrust::device_vector<Key> d_results(num_keys);

    set.find(keys_begin, keys_begin + num_keys, d_results.begin());
    auto zip = thrust::make_zip_iterator(cuda::std::tuple{d_results.begin(), keys_begin});

    REQUIRE(cuco::test::all_of(zip, zip + num_keys, zip_equal));
  }

  SECTION("Conditional find should return valid values on even inputs.")
  {
    auto found_results = thrust::device_vector<Key>(num_keys);
    auto is_even =
      cuda::proclaim_return_type<bool>([] __device__(auto const& i) { return i % 2 == 0; });
    auto gold_fn = cuda::proclaim_return_type<Key>(
      [] __device__(auto const& i) { return i % 2 == 0 ? static_cast<Key>(i) : Key{SENTINEL}; });

    set.find_if(keys_begin,
                keys_begin + num_keys,
                thrust::counting_iterator<std::size_t>{0},
                is_even,
                found_results.begin());

    REQUIRE(cuco::test::equal(
      found_results.begin(),
      found_results.end(),
      thrust::make_transform_iterator(thrust::counting_iterator<Key>{0}, gold_fn),
      cuda::proclaim_return_type<bool>(
        [] __device__(auto const& found, auto const& gold) { return found == gold; })));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_multiset find tests",
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
  constexpr size_type num_keys{1'000};

  using probe = std::conditional_t<Probe == cuco::test::probe_sequence::linear_probing,
                                   cuco::linear_probing<CGSize, cuco::default_hash_function<Key>>,
                                   cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>>;

  auto set = cuco::static_multiset{
    num_keys, cuco::empty_key<Key>{SENTINEL}, {}, probe{}, {}, cuco::storage<2>{}};

  test_unique_sequence(set, num_keys);
}
