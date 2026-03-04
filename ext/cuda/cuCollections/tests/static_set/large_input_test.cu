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
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

#include <catch2/catch_template_test_macros.hpp>

template <typename Set>
void test_unique_sequence(Set& set, bool* res_begin, std::size_t num_keys)
{
  using Key = typename Set::key_type;

  auto const keys_begin = thrust::counting_iterator<Key>(0);
  auto const keys_end   = thrust::counting_iterator<Key>(num_keys);

  SECTION("Non-inserted keys should not be contained.")
  {
    REQUIRE(set.size() == 0);

    set.contains(keys_begin, keys_end, res_begin);
    REQUIRE(cuco::test::none_of(res_begin, res_begin + num_keys, cuda::std::identity{}));
  }

  set.insert(keys_begin, keys_end);
  REQUIRE(set.size() == num_keys);

  SECTION("All inserted key/value pairs should be contained.")
  {
    set.contains(keys_begin, keys_end, res_begin);
    REQUIRE(cuco::test::all_of(res_begin, res_begin + num_keys, cuda::std::identity{}));
  }

  SECTION("All inserted key/value pairs can be retrieved.")
  {
    auto output_keys = thrust::device_vector<Key>(num_keys);

    auto const keys_end = set.retrieve_all(output_keys.begin());
    REQUIRE(static_cast<std::size_t>(std::distance(output_keys.begin(), keys_end)) == num_keys);

    thrust::sort(output_keys.begin(), keys_end);

    REQUIRE(cuco::test::equal(output_keys.begin(),
                              output_keys.end(),
                              thrust::counting_iterator<Key>(0),
                              cuda::std::equal_to<Key>{}));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "cuco::static_set large input test",
  "",
  ((typename Key, cuco::test::probe_sequence Probe, int CGSize), Key, Probe, CGSize),
  (int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, cuco::test::probe_sequence::double_hashing, 2))
{
  constexpr std::size_t num_keys{1'200'000'000};

  using extent_type = cuco::extent<std::size_t>;
  using probe       = cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>;

  try {
    auto set = cuco::static_set{num_keys * 2, cuco::empty_key<Key>{-1}, {}, probe{}};

    thrust::device_vector<bool> d_contained(num_keys);
    test_unique_sequence(set, d_contained.data().get(), num_keys);
  } catch (cuco::cuda_error&) {
    SKIP("Out of memory");
  } catch (std::bad_alloc&) {
    SKIP("Out of memory");
  }
}
