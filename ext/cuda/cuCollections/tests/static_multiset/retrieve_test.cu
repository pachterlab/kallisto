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

#include <cuco/static_multiset.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/distance.h>
#include <thrust/functional.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/sort.h>

#include <catch2/catch_template_test_macros.hpp>

#include <limits>

template <class Container>
void test_multiplicity(Container& container, std::size_t num_keys, std::size_t multiplicity)
{
  using key_type                = typename Container::key_type;
  auto const empty_key_sentinel = container.empty_key_sentinel();

  container.clear();

  auto const num_unique_keys = num_keys / multiplicity;
  REQUIRE(num_unique_keys > 0);
  auto const num_actual_keys = num_unique_keys * multiplicity;
  REQUIRE(num_actual_keys <= num_keys);

  thrust::device_vector<key_type> probed_keys(num_actual_keys);
  thrust::device_vector<key_type> matched_keys(num_actual_keys);

  auto const keys_begin = thrust::make_transform_iterator(
    thrust::counting_iterator<key_type>(0),
    cuda::proclaim_return_type<key_type>([multiplicity] __device__(auto const& i) {
      return static_cast<key_type>(i / multiplicity);
    }));

  container.insert(keys_begin, keys_begin + num_actual_keys);
  REQUIRE(container.size() == num_actual_keys);

  SECTION("All inserted keys should be contained.")
  {
    auto const [probed_end, matched_end] = container.retrieve(
      keys_begin, keys_begin + num_actual_keys, probed_keys.begin(), matched_keys.begin());
    thrust::sort(probed_keys.begin(), probed_end);
    thrust::sort(matched_keys.begin(), matched_end);
    REQUIRE(cuco::test::equal(
      probed_keys.begin(), probed_keys.end(), keys_begin, cuda::std::equal_to<key_type>{}));
    REQUIRE(cuco::test::equal(
      matched_keys.begin(), matched_keys.end(), keys_begin, cuda::std::equal_to<key_type>{}));
  }
}

template <class Container>
void test_outer(Container& container, std::size_t num_keys)
{
  using key_type                = typename Container::key_type;
  auto const empty_key_sentinel = container.empty_key_sentinel();

  container.clear();

  auto const keys_begin = thrust::counting_iterator<key_type>{0};
  auto const query_size = num_keys * 2ull;

  thrust::device_vector<key_type> probed_keys(num_keys * 2ull);
  thrust::device_vector<key_type> matched_keys(num_keys * 2ull);

  SECTION("Non-inserted keys should output sentinels.")
  {
    auto const [probed_end, matched_end] = container.retrieve_outer(keys_begin,
                                                                    keys_begin + query_size,
                                                                    container.key_eq(),
                                                                    container.hash_function(),
                                                                    probed_keys.begin(),
                                                                    matched_keys.begin());
    REQUIRE(static_cast<std::size_t>(std::distance(probed_keys.begin(), probed_end)) ==
            num_keys * 2ull);
    REQUIRE(static_cast<std::size_t>(std::distance(matched_keys.begin(), matched_end)) ==
            num_keys * 2ull);
    REQUIRE(cuco::test::all_of(
      matched_keys.begin(),
      matched_keys.end(),
      cuda::proclaim_return_type<bool>([empty_key_sentinel] __device__(auto const& k) {
        return static_cast<bool>(k == static_cast<key_type>(empty_key_sentinel));
      })));
  }

  container.insert(keys_begin, keys_begin + num_keys);

  SECTION("All inserted keys should be contained.")
  {
    auto const [probed_end, matched_end] = container.retrieve_outer(keys_begin,
                                                                    keys_begin + query_size,
                                                                    container.key_eq(),
                                                                    container.hash_function(),
                                                                    probed_keys.begin(),
                                                                    matched_keys.begin());
    thrust::sort_by_key(
      probed_keys.begin(), probed_end, matched_keys.begin(), cuda::std::less<key_type>());

    REQUIRE(cuco::test::equal(
      probed_keys.begin(), probed_keys.end(), keys_begin, cuda::std::equal_to<key_type>{}));
    REQUIRE(cuco::test::equal(matched_keys.begin(),
                              matched_keys.begin() + num_keys,
                              keys_begin,
                              cuda::std::equal_to<key_type>{}));
    REQUIRE(cuco::test::all_of(
      matched_keys.begin() + num_keys,
      matched_keys.end(),
      cuda::proclaim_return_type<bool>([empty_key_sentinel] __device__(auto const& k) {
        return static_cast<bool>(k == static_cast<key_type>(empty_key_sentinel));
      })));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_multiset retrieve tests",
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
  constexpr std::size_t num_keys{400};
  constexpr double desired_load_factor = 0.5;
  constexpr auto empty_key_sentinel    = std::numeric_limits<Key>::max();

  using probe = std::conditional_t<Probe == cuco::test::probe_sequence::linear_probing,
                                   cuco::linear_probing<CGSize, cuco::default_hash_function<Key>>,
                                   cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>>;

  auto set = cuco::static_multiset{
    num_keys, desired_load_factor, cuco::empty_key<Key>{empty_key_sentinel}, {}, probe{}};

  test_multiplicity(set, num_keys, 1);  // unique sequence
  test_multiplicity(set, num_keys, 2);  // each key occurs twice
  test_multiplicity(set, num_keys, 11);
  test_outer(set, num_keys);
}
