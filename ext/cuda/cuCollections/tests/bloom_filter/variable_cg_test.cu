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

#include <cuco/bloom_filter.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/fill.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cstdint>
#include <exception>

using size_type = int32_t;

template <int32_t AddCGSize, int32_t ContainsCGSize, typename Filter>
void test_variable_cg_size(Filter& filter, size_type num_keys)
{
  constexpr int32_t block_size = 128;
  constexpr int32_t grid_size  = 128;

  using Key = typename Filter::key_type;

  auto ref = filter.ref();

  // Generate keys
  thrust::device_vector<Key> keys(num_keys);
  thrust::sequence(thrust::device, keys.begin(), keys.end());

  thrust::device_vector<bool> contained(num_keys, false);

  auto const always_true = thrust::constant_iterator<bool>{true};

  SECTION("Check if fallback kernels work for varying combinations of CG sizes.")
  {
    cuco::detail::bloom_filter_ns::add_if_n<AddCGSize, block_size>
      <<<grid_size, block_size>>>(keys.begin(), num_keys, always_true, cuda::std::identity{}, ref);
    cuco::detail::bloom_filter_ns::contains_if_n<ContainsCGSize, block_size>
      <<<grid_size, block_size>>>(
        keys.begin(), num_keys, always_true, cuda::std::identity{}, contained.begin(), ref);
    REQUIRE(cuco::test::all_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  filter.clear();
  thrust::fill(contained.begin(), contained.end(), false);  // reset output vector

  SECTION("Check if adaptive add kernel works with fallback contains kernel.")
  {
    cuco::detail::bloom_filter_ns::add<block_size>
      <<<grid_size, block_size>>>(keys.begin(), num_keys, ref);
    cuco::detail::bloom_filter_ns::contains_if_n<ContainsCGSize, block_size>
      <<<grid_size, block_size>>>(
        keys.begin(), num_keys, always_true, cuda::std::identity{}, contained.begin(), ref);
    REQUIRE(cuco::test::all_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  // TODO adaptive vs. adaptive and fallback add vs. adaptive contains (requires #673)
}

TEMPLATE_TEST_CASE_SIG(
  "bloom_filter variable CG size tests",
  "",
  ((int32_t AddCGSize, int32_t ContainsCGSize, class Key, class Policy),
   AddCGSize,
   ContainsCGSize,
   Key,
   Policy),
  (1, 4, int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint32_t, 1>),
  (1, 4, int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint32_t, 8>),
  (1, 4, int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint64_t, 1>),
  (1, 4, int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint64_t, 8>),
  (4, 1, int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint32_t, 1>),
  (4, 1, int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint32_t, 8>),
  (4, 1, int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint64_t, 1>),
  (4, 1, int32_t, cuco::default_filter_policy<cuco::xxhash_64<int32_t>, uint64_t, 8>))
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

  test_variable_cg_size<AddCGSize, ContainsCGSize>(filter, num_keys);
}
