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

#include <cuco/static_set.cuh>

#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>

#include <catch2/catch_template_test_macros.hpp>

template <typename Set>
void test_unique_sequence(Set& set, std::size_t num_keys)
{
  using Key = typename Set::key_type;

  thrust::device_vector<Key> d_keys(num_keys);
  thrust::sequence(d_keys.begin(), d_keys.end());
  auto keys_begin = d_keys.begin();

  SECTION("Non-inserted keys should not be contained.")
  {
    REQUIRE(set.size() == 0);

    auto keys_end = set.retrieve_all(keys_begin);
    REQUIRE(std::distance(keys_begin, keys_end) == 0);
  }

  set.insert(keys_begin, keys_begin + num_keys);
  REQUIRE(set.size() == num_keys);

  SECTION("All inserted key/value pairs should be contained.")
  {
    thrust::device_vector<Key> d_res(num_keys);
    auto d_res_end = set.retrieve_all(d_res.begin());
    thrust::sort(d_res.begin(), d_res_end);
    REQUIRE(cuco::test::equal(
      d_res.begin(), d_res_end, thrust::counting_iterator<Key>(0), cuda::std::equal_to<Key>{}));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_set::retrieve_all tests",
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

  constexpr std::size_t gold_capacity = [&]() {
    if constexpr (cuco::is_double_hashing<probe>::value) {
      return (CGSize == 1) ? 409   // 409 x 1 x 2
                           : 422;  // 211 x 2 x 2
    } else {
      return 400;
    }
  }();

  auto set = cuco::static_set{num_keys, desired_load_factor, cuco::empty_key<Key>{-1}, {}, probe{}};

  REQUIRE(set.capacity() == gold_capacity);

  test_unique_sequence(set, num_keys);
}
