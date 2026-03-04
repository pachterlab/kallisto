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

#include <cuco/hash_functions.cuh>
#include <cuco/hyperloglog.cuh>

#include <thrust/device_vector.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>

TEMPLATE_TEST_CASE_SIG("hyperloglog: unique sequence",
                       "",
                       ((typename T, typename Hash), T, Hash),
                       (int32_t, cuco::xxhash_64<int32_t>),
                       (int64_t, cuco::xxhash_64<int64_t>),
                       (__int128_t, cuco::xxhash_64<__int128_t>))
{
  auto num_items_pow2 = GENERATE(25, 26, 28);
  auto hll_precision  = GENERATE(8, 10, 12, 13, 18, 20);
  auto sketch_size_kb = 4 * (1ull << hll_precision) / 1024;
  INFO("hll_precision=" << hll_precision);
  INFO("sketch_size_kb=" << sketch_size_kb);
  INFO("num_items=2^" << num_items_pow2);
  auto num_items = 1ull << num_items_pow2;

  // This factor determines the error threshold for passing the test
  double constexpr tolerance_factor = 2.5;
  // RSD for a given precision is given by the following formula
  double const relative_standard_deviation =
    1.04 / std::sqrt(static_cast<double>(1ull << hll_precision));

  thrust::device_vector<T> items(num_items);

  // Generate `num_items` distinct items
  thrust::sequence(items.begin(), items.end(), 0);

  // Initialize the estimator
  cuco::hyperloglog<T, cuda::thread_scope_device, Hash> estimator{
    cuco::sketch_size_kb(sketch_size_kb)};

  REQUIRE(estimator.estimate() == 0);

  // Add all items to the estimator
  estimator.add(items.begin(), items.end());

  auto const estimate = estimator.estimate();

  // Adding the same items again should not affect the result
  estimator.add(items.begin(), items.begin() + num_items / 2);
  REQUIRE(estimator.estimate() == estimate);

  // Clearing the estimator should reset the estimate
  estimator.clear();
  REQUIRE(estimator.estimate() == 0);

  double const relative_error =
    std::abs((static_cast<double>(estimate) / static_cast<double>(num_items)) - 1.0);

  // Check if the error is acceptable
  REQUIRE(relative_error < tolerance_factor * relative_standard_deviation);
}
