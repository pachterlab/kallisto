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

#include <cuco/static_set.cuh>

#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

#include <catch2/catch_template_test_macros.hpp>

static constexpr int key_sentinel = -1;

template <typename Set>
void test_unique_sequence(Set& set, std::size_t num_keys)
{
  using Key = typename Set::key_type;

  thrust::device_vector<Key> keys(num_keys);
  thrust::device_vector<Key> matched_keys(num_keys);

  auto iter = thrust::counting_iterator<Key>{0};

  SECTION("Non-inserted keys should not be contained.")
  {
    REQUIRE(set.size() == 0);

    REQUIRE(set.count(iter, iter + num_keys) == 0);

    auto const [probe_end, matched_end] =
      set.retrieve(iter, iter + num_keys, keys.begin(), matched_keys.begin());
    REQUIRE(std::distance(keys.begin(), probe_end) == 0);
    REQUIRE(std::distance(matched_keys.begin(), matched_end) == 0);
  }

  set.insert(iter, iter + num_keys);

  SECTION("All inserted key/value pairs should be contained.")
  {
    REQUIRE(set.count(iter, iter + num_keys) == num_keys);

    auto const [probe_end, matched_end] =
      set.retrieve(iter, iter + num_keys, keys.begin(), matched_keys.begin());
    thrust::sort(keys.begin(), probe_end);
    thrust::sort(matched_keys.begin(), matched_end);
    REQUIRE(cuco::test::equal(
      keys.begin(), keys.end(), thrust::counting_iterator<Key>(0), cuda::std::equal_to<Key>{}));
    REQUIRE(cuco::test::equal(matched_keys.begin(),
                              matched_keys.end(),
                              thrust::counting_iterator<Key>(0),
                              cuda::std::equal_to<Key>{}));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_set retrieve tests",
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
  constexpr double desired_load_factor = 1.;

  using probe = std::conditional_t<Probe == cuco::test::probe_sequence::linear_probing,
                                   cuco::linear_probing<CGSize, cuco::default_hash_function<Key>>,
                                   cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>>;

  auto set = cuco::static_set{
    num_keys, desired_load_factor, cuco::empty_key<Key>{key_sentinel}, {}, probe{}};

  test_unique_sequence(set, num_keys);
}
