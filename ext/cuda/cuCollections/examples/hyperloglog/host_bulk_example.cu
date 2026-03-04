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
#include <cuco/hyperloglog.cuh>

#include <thrust/device_vector.h>
#include <thrust/sequence.h>

#include <cmath>
#include <cstddef>
#include <iostream>

/**
 * @file host_bulk_example.cu
 * @brief Demonstrates usage of `cuco::hyperloglog` "bulk" host APIs.
 */
int main(void)
{
  using T                         = int;
  constexpr std::size_t num_items = 1ull << 28;  // 1GB

  thrust::device_vector<T> items(num_items);

  // Generate `num_items` distinct items
  thrust::sequence(items.begin(), items.end(), 0);

  // We define the desired standard deviation of the approximation error
  // 0.0122197 is the default value and corresponds to a 32KB sketch size
  auto const sd = cuco::standard_deviation{0.0122197};

  // Initialize the estimator
  cuco::hyperloglog<T> estimator{sd};

  // Add all items to the estimator
  estimator.add(items.begin(), items.end());

  // Adding the same items again will not affect the result
  estimator.add(items.begin(), items.begin() + num_items / 2);

  // Calculate the cardinality estimate
  std::size_t const estimated_cardinality = estimator.estimate();

  std::cout << "True cardinality: " << num_items
            << "\nEstimated cardinality: " << estimated_cardinality << "\nError: "
            << std::abs(
                 static_cast<double>(estimated_cardinality) / static_cast<double>(num_items) - 1.0)
            << std::endl;

  return 0;
}