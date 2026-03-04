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

#include <cuco/bloom_filter.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <exception>

using size_type = int32_t;

template <typename Filter>
void test_unique_sequence(Filter& filter, size_type num_keys)
{
  using Key = typename Filter::key_type;

  // Generate keys
  thrust::device_vector<Key> keys(num_keys);
  thrust::sequence(thrust::device, keys.begin(), keys.end());

  thrust::device_vector<bool> contained(num_keys, false);

  auto is_even =
    cuda::proclaim_return_type<bool>([] __device__(auto const& i) { return i % 2 == 0; });

  SECTION("Non-inserted keys should not be contained.")
  {
    filter.contains(keys.begin(), keys.end(), contained.begin());
    REQUIRE(cuco::test::none_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  SECTION("All inserted keys should be contained.")
  {
    filter.add(keys.begin(), keys.end());
    filter.contains(keys.begin(), keys.end(), contained.begin());
    REQUIRE(cuco::test::all_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  SECTION("After clearing the filter no keys should be contained.")
  {
    filter.clear();
    filter.contains(keys.begin(), keys.end(), contained.begin());
    REQUIRE(cuco::test::none_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  SECTION("All conditionally inserted keys should be contained")
  {
    filter.add_if(keys.begin(), keys.end(), thrust::counting_iterator<std::size_t>(0), is_even);
    filter.contains_if(keys.begin(),
                       keys.end(),
                       thrust::counting_iterator<std::size_t>(0),
                       is_even,
                       contained.begin());
    REQUIRE(cuco::test::equal(
      contained.begin(),
      contained.end(),
      thrust::counting_iterator<std::size_t>(0),
      cuda::proclaim_return_type<bool>([] __device__(auto const& idx_contained, auto const& idx) {
        return ((idx % 2) == 0) == idx_contained;
      })));
  }

  // TODO test FPR but how?
}

TEMPLATE_TEST_CASE_SIG(
  "bloom_filter default policy tests",
  "",
  ((class Key, class Policy), Key, Policy),
  (int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint32_t, 1>),
  (int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint32_t, 8>),
  (int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint64_t, 1>),
  (int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint64_t, 8>))
{
  using filter_type =
    cuco::bloom_filter<Key, cuco::extent<size_t>, cuda::thread_scope_device, Policy>;
  constexpr size_type num_keys{400};

  uint32_t pattern_bits = Policy::words_per_block + GENERATE(0, 1, 2, 3, 4);

  // some parameter combinations might be invalid so we skip them
  try {
    [[maybe_unused]] auto policy = Policy{pattern_bits};
  } catch (std::exception const& e) {
    SKIP(e.what());
  }

  auto filter = filter_type{1000, {}, {pattern_bits}};

  test_unique_sequence(filter, num_keys);
}

TEMPLATE_TEST_CASE_SIG("bloom_filter arrow policy tests",
                       "",
                       ((class Key, class Policy), Key, Policy),
                       (int32_t, cuco::arrow_filter_policy<int32_t>),
                       (uint64_t, cuco::arrow_filter_policy<uint64_t>),
                       (float, cuco::arrow_filter_policy<float>))
{
  using filter_type =
    cuco::bloom_filter<Key, cuco::extent<size_t>, cuda::thread_scope_device, Policy>;
  constexpr size_type num_keys{400};

  auto filter = filter_type{1000};

  test_unique_sequence(filter, num_keys);
}
