/*
 * Copyright (c) 2022-2025, NVIDIA CORPORATION.
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

#include <cuco/static_set.cuh>

#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/logical.h>
#include <thrust/sequence.h>

#include <iostream>
#include <limits>

/**
 * @file host_bulk_example.cu
 * @brief Demonstrates usage of the static_set "bulk" host APIs.
 *
 * The bulk APIs are only invocable from the host and are used for doing operations like `insert` or
 * `contains` on a set of keys.
 *
 */
int main(void)
{
  using Key = int;

  // Empty slots are represented by reserved "sentinel" values. These values should be selected such
  // that they never occur in your input data.
  Key constexpr empty_key_sentinel = -1;

  // Number of keys to be inserted
  std::size_t constexpr num_keys = 50'000;

  // Compute capacity based on a 50% load factor
  auto constexpr load_factor = 0.5;
  std::size_t const capacity = std::ceil(num_keys / load_factor);

  // Constructs a set with at least `capacity` slots using -1 as the empty keys sentinel.
  cuco::static_set<Key> set{capacity, cuco::empty_key{empty_key_sentinel}};

  // Create a sequence of keys {0, 1, 2, .., i}
  thrust::device_vector<Key> keys(num_keys);
  thrust::sequence(keys.begin(), keys.end(), 0);

  // Inserts all keys into the hash set
  set.insert(keys.begin(), keys.end());

  // Storage for result
  thrust::device_vector<bool> found(num_keys);

  // Check if all keys are contained in the set
  set.contains(keys.begin(), keys.end(), found.begin());

  // Verify that all keys have been found
  bool const all_keys_found = thrust::all_of(found.begin(), found.end(), cuda::std::identity{});

  if (all_keys_found) { std::cout << "Success! Found all keys.\n"; }

  return 0;
}
