/*
 * Copyright (c) 2020-2025, NVIDIA CORPORATION.
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

#include <cuco/static_map.cuh>

#include <thrust/device_vector.h>
#include <thrust/equal.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>

/**
 * @file host_bulk_example.cu
 * @brief Demonstrates usage of the static_map "bulk" host APIs.
 *
 * The bulk APIs are only invocable from the host and are used for doing operations like insert or
 * find on a set of keys.
 *
 */

int main(void)
{
  using Key   = int;
  using Value = int;

  // Empty slots are represented by reserved "sentinel" values. These values should be selected such
  // that they never occur in your input data.
  Key constexpr empty_key_sentinel     = -1;
  Value constexpr empty_value_sentinel = -1;

  // Number of key/value pairs to be inserted
  std::size_t constexpr num_keys = 50'000;

  // Compute capacity based on a 50% load factor
  auto constexpr load_factor = 0.5;
  std::size_t const capacity = std::ceil(num_keys / load_factor);

  // Constructs a map with "capacity" slots using -1 and -1 as the empty key/value sentinels.
  auto map = cuco::static_map{
    capacity, cuco::empty_key{empty_key_sentinel}, cuco::empty_value{empty_value_sentinel}};

  // Create a sequence of keys and values
  thrust::device_vector<Key> insert_keys(num_keys);
  thrust::sequence(insert_keys.begin(), insert_keys.end(), 0);
  thrust::device_vector<Value> insert_values(num_keys);
  thrust::sequence(insert_values.begin(), insert_values.end(), 0);
  // Combine keys and values into pairs {{0,0}, {1,1}, ... {i,i}}
  auto pairs = thrust::make_transform_iterator(
    thrust::counting_iterator<std::size_t>{0},
    cuda::proclaim_return_type<cuco::pair<Key, Value>>(
      [keys = insert_keys.begin(), values = insert_values.begin()] __device__(auto i) {
        return cuco::pair<Key, Value>{keys[i], values[i]};
      }));

  // Inserts all pairs into the map
  map.insert(pairs, pairs + num_keys);

  // Storage for found values
  thrust::device_vector<Value> found_values(num_keys);

  // Finds all keys {0, 1, 2, ...} and stores associated values into `found_values`
  // If a key `keys_to_find[i]` doesn't exist, `found_values[i] == empty_value_sentinel`
  map.find(insert_keys.begin(), insert_keys.end(), found_values.begin());

  // Verify that all the found values match the inserted values
  bool const all_values_match =
    thrust::equal(found_values.begin(), found_values.end(), insert_values.begin());

  if (all_values_match) { std::cout << "Success! Found all values.\n"; }

  return 0;
}
