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

#include <cuco/static_multiset.cuh>

#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/logical.h>
#include <thrust/sequence.h>

#include <iostream>
#include <limits>

/**
 * @file host_bulk_example.cu
 * @brief Demonstrates usage of the static_multiset "bulk" host APIs.
 *
 * The bulk APIs are only invocable from the host and are used for doing operations like `insert` or
 * `retrieve` on a multiset of keys.
 *
 */
int main(void)
{
  using key_type = int;

  // Empty slots are represented by reserved "sentinel" values. These values should be selected such
  // that they never occur in your input data.
  key_type constexpr empty_key_sentinel = -1;

  // Number of keys to be inserted
  std::size_t constexpr num_keys = 50'000;

  // Compute capacity based on a 50% load factor
  auto constexpr load_factor = 0.5;
  std::size_t const capacity = std::ceil(num_keys / load_factor);

  // Constructs a set with at least `capacity` slots using -1 as the empty keys sentinel.
  cuco::static_multiset<key_type> multiset{capacity, cuco::empty_key{empty_key_sentinel}};

  // Create a sequence of keys {0, 1, 2, .., i}
  // We're going to insert each key twice so we only need 'num_keys / 2' distinct keys.
  thrust::device_vector<key_type> keys(num_keys / 2);
  thrust::sequence(keys.begin(), keys.end(), 0);

  // Inserts all keys into the hash set
  multiset.insert(keys.begin(), keys.end());
  // Insert the same set of keys again, so each distinct key should occur twice in the multiset
  multiset.insert(keys.begin(), keys.end());

  // Counts the occurrences of matching keys contained in the multiset.
  std::size_t const counted_output_size = multiset.count(keys.begin(), keys.end());

  // Storage for result
  thrust::device_vector<key_type> output_probes(counted_output_size);
  thrust::device_vector<key_type> output_matches(counted_output_size);

  // Retrieve all matching keys
  auto const [output_probes_end, _] =
    multiset.retrieve(keys.begin(), keys.end(), output_probes.begin(), output_matches.begin());
  std::size_t const retrieved_output_size = output_probes_end - output_probes.begin();

  if ((retrieved_output_size == counted_output_size) and (retrieved_output_size == num_keys)) {
    std::cout << "Success! Found all keys.\n";
  } else {
    std::cout << "Fail! Something went wrong.\n";
  }

  return 0;
}
