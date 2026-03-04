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

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>  // std::memcpy

/**
 * @file spark_parity_test.cu
 * @brief Unit test to ensure parity with Spark's HLL implementation
 *
 * The following unit tests mimic Spark's unit tests which can be found here:
 * https://github.com/apache/spark/blob/d10dbaa31a44878df5c7e144f111e18261346531/sql/catalyst/src/test/scala/org/apache/spark/sql/catalyst/expressions/aggregate/HyperLogLogPlusPlusSuite.scala
 *
 */

// TODO implement this test once add_if is available
// TEST_CASE("hyperloglog: Spark parity: add nulls", "")

TEST_CASE("hyperloglog: Spark parity: deterministic cardinality estimation", "")
{
  using T              = int;
  using estimator_type = cuco::hyperloglog<T, cuda::thread_scope_device, cuco::xxhash_64<T>>;

  constexpr size_t repeats = 10;
  // This factor determines the error threshold for passing the test
  constexpr double tolerance_factor = 3.0;
  auto num_items          = GENERATE(100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000);
  auto standard_deviation = GENERATE(0.1, 0.05, 0.025, 0.01, 0.001);

  auto expected_hll_precision = std::max(
    static_cast<int32_t>(4),
    static_cast<int32_t>(std::ceil(2.0 * std::log(1.106 / standard_deviation) / std::log(2.0))));
  auto expected_sketch_bytes = 4 * (1ull << expected_hll_precision);

  INFO("num_items" << num_items);
  INFO("standard_deviation=" << standard_deviation);
  INFO("expected_hll_precision=" << expected_hll_precision);
  INFO("expected_sketch_bytes=" << expected_sketch_bytes);

  auto sd = cuco::standard_deviation(standard_deviation);
  auto sb = cuco::sketch_size_kb(expected_sketch_bytes / 1024.0);

  // Validate sketch size calculation
  REQUIRE(estimator_type::sketch_bytes(sd) >= 64);
  REQUIRE(estimator_type::sketch_bytes(sd) == expected_sketch_bytes);
  REQUIRE(estimator_type::sketch_bytes(sd) == estimator_type::sketch_bytes(sb));

  auto items_begin =
    thrust::make_transform_iterator(thrust::make_counting_iterator<size_t>(0),
                                    cuda::proclaim_return_type<T>([repeats] __device__(auto i) {
                                      return static_cast<T>(i / repeats);
                                    }));

  estimator_type estimator{sd};

  REQUIRE(estimator.estimate() == 0);

  // Add all items to the estimator
  estimator.add(items_begin, items_begin + num_items);

  auto const estimate = estimator.estimate();

  double const relative_error =
    std::abs((static_cast<double>(estimate) / static_cast<double>(num_items / repeats)) - 1.0);
  // RSD for a given precision is given by the following formula
  double const expected_standard_deviation =
    1.04 / std::sqrt(static_cast<double>(1ull << expected_hll_precision));

  // Check if the error is acceptable
  REQUIRE(relative_error < expected_standard_deviation * tolerance_factor);
}

// the following test is omitted since we refrain from doing randomized unit tests in cuco
// TEST_CASE("hyperloglog: Spark parity: random cardinality estimation", "")

TEST_CASE("hyperloglog: Spark parity: merging HLL instances", "")
{
  using T              = int;
  using estimator_type = cuco::hyperloglog<T, cuda::thread_scope_device, cuco::xxhash_64<T>>;

  auto num_items          = 1000000;
  auto standard_deviation = cuco::standard_deviation(0.05);

  auto items_begin = thrust::make_counting_iterator<T>(0);

  // count lower half of input
  estimator_type lower{standard_deviation};
  lower.add(items_begin, items_begin + num_items / 2);

  // count upper half of input
  estimator_type upper{standard_deviation};
  upper.add(items_begin + num_items / 2, items_begin + num_items);

  // merge upper into lower so lower has seen the entire input
  lower.merge(upper);

  auto reversed_items_begin = thrust::make_transform_iterator(
    items_begin, cuda::proclaim_return_type<T>([num_items] __device__(auto i) {
      return static_cast<T>(num_items - i);
    }));

  // count the entire input vector but in reversed order
  estimator_type entire{standard_deviation};
  entire.add(reversed_items_begin, reversed_items_begin + num_items);

  auto const entire_sketch = entire.sketch();
  auto const lower_sketch  = lower.sketch();

  // check if sketches are bitwise identical
  REQUIRE(cuco::test::equal(entire_sketch.data(),
                            entire_sketch.data() + entire_sketch.size(),
                            lower_sketch.data(),
                            cuda::std::equal_to{}));
}

/*
The following unit tests fail since xxhash_64 does not deduplicate different bit patterns for NaN
values and +-0.0. They are thus counted as distinct items.

TEST_CASE("hyperloglog: Spark parity: add 0.0 and -0.0", "")
{
  using T = double;
  using estimator_type =
    cuco::hyperloglog<T, cuda::thread_scope_device, cuco::xxhash_64<T>>;

  auto standard_deviation = cuco::standard_deviation(0.05);

  auto items = thrust::device_vector<T>({0.0, -0.0});

  estimator_type estimator{standard_deviation};
  estimator.add(items.begin(), items.end());

  REQUIRE(estimator.estimate() == 1);
}

TEST_CASE("hyperloglog: Spark parity: add NaN", "")
{
  using T = double;
  using estimator_type =
    cuco::hyperloglog<T, cuda::thread_scope_device, cuco::xxhash_64<T>>;

  auto standard_deviation = cuco::standard_deviation(0.05);

  // Define the special bit pattern for the NaN.
  uint64_t nan_bits = 0x7ff1234512345678ULL;
  double special_nan;
  std::memcpy(&special_nan, &nan_bits, sizeof(special_nan));

  auto items = thrust::device_vector<T>({0.0, special_nan});

  estimator_type estimator{standard_deviation};
  estimator.add(items.begin(), items.end());

  REQUIRE(estimator.estimate() == 1);
}
*/
