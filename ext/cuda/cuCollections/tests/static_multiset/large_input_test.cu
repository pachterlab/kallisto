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

#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/discard_iterator.h>
#include <thrust/sort.h>

#include <catch2/catch_template_test_macros.hpp>

#include <cstdint>
#include <iterator>

template <typename Set>
void test_unique_sequence(Set& set, typename Set::value_type* res_begin, std::size_t num_keys)
{
  using Key = typename Set::key_type;

  auto const keys_begin = thrust::counting_iterator<Key>(0);
  auto const keys_end   = keys_begin + num_keys;

  set.insert(keys_begin, keys_end);
  REQUIRE(set.size() == num_keys);

  SECTION("All inserted keys can be retrieved.")
  {
    auto const [_, res_end] =
      set.retrieve(keys_begin, keys_end, thrust::make_discard_iterator(), res_begin);
    REQUIRE(static_cast<std::size_t>(std::distance(res_begin, res_end)) == num_keys);

    thrust::sort(thrust::device, res_begin, res_end);

    REQUIRE(cuco::test::equal(res_begin, res_end, keys_begin, cuda::std::equal_to<Key>{}));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "cuco::static_multiset large input test",
  "",
  ((typename Key, cuco::test::probe_sequence Probe, int CGSize), Key, Probe, CGSize),
  (int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, cuco::test::probe_sequence::double_hashing, 2))
{
  constexpr std::size_t num_keys{1'200'000'000};

  using extent_type = cuco::extent<std::size_t>;
  using probe       = cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>;

  try {
    auto set = cuco::static_multiset{num_keys * 2, cuco::empty_key<Key>{-1}, {}, probe{}};

    thrust::device_vector<Key> d_retrieved(num_keys);
    test_unique_sequence(set, d_retrieved.data().get(), num_keys);
  } catch (cuco::cuda_error&) {
    SKIP("Out of memory");
  } catch (std::bad_alloc&) {
    SKIP("Out of memory");
  }
}
