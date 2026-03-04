/*
 * Copyright (c) 2020-2024, NVIDIA CORPORATION.
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

#include <cuda/std/functional>
#include <cuda/std/tuple>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sequence.h>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>

/**
 * @file device_ref_example.cu
 * @brief Demonstrates usage of the device side APIs for individual operations like insert/find.
 *
 * Individual operations like a single insert or find can be performed in device code via the
 * "static_map_ref" types.
 *
 * @note This example is for demonstration purposes only. It is not intended to show the most
 * performant way to do the example algorithm.
 *
 */

/**
 * @brief Inserts keys that pass the specified predicated into the map.
 *
 * @tparam Map Type of the map device reference
 * @tparam KeyIter Input iterator whose value_type convertible to Map::key_type
 * @tparam ValueIter Input iterator whose value_type is convertible to Map::mapped_type
 * @tparam Predicate Unary predicate
 *
 * @param[in] map_ref Reference of the map into which inserts will be performed
 * @param[in] key_begin The beginning of the range of keys to insert
 * @param[in] value_begin The beginning of the range of values associated with each key to insert
 * @param[in] num_keys The total number of keys and values
 * @param[in] pred Unary predicate applied to each key. Only keys that pass the predicated will be
 * inserted.
 * @param[out] num_inserted The total number of keys successfully inserted
 */
template <typename Map, typename KeyIter, typename ValueIter, typename Predicate>
__global__ void filtered_insert(Map map_ref,
                                KeyIter key_begin,
                                ValueIter value_begin,
                                std::size_t num_keys,
                                Predicate pred,
                                int* num_inserted)
{
  auto tid = threadIdx.x + blockIdx.x * blockDim.x;

  std::size_t counter = 0;
  while (tid < num_keys) {
    // Only insert keys that pass the predicate
    if (pred(key_begin[tid])) {
      // Map::insert returns `true` if it is the first time the given key was
      // inserted and `false` if the key already existed
      if (map_ref.insert(cuco::pair{key_begin[tid], value_begin[tid]})) {
        ++counter;  // Count number of successfully inserted keys
      }
    }
    tid += gridDim.x * blockDim.x;
  }

  // Update global count of inserted keys
  atomicAdd(num_inserted, counter);
}

/**
 * @brief For keys that have a match in the map, increments their corresponding value by one.
 *
 * @tparam Map Type of the map device reference
 * @tparam KeyIter Input iterator whose value_type convertible to Map::key_type
 *
 * @param map_ref Reference of the map into which queries will be performed
 * @param key_begin The beginning of the range of keys to query
 * @param num_keys The total number of keys
 */
template <typename Map, typename KeyIter>
__global__ void increment_values(Map map_ref, KeyIter key_begin, std::size_t num_keys)
{
  auto tid = threadIdx.x + blockIdx.x * blockDim.x;
  while (tid < num_keys) {
    // If the key exists in the map, find returns an iterator to the specified key. Otherwise it
    // returns map.end()
    auto found = map_ref.find(key_begin[tid]);
    if (found != map_ref.end()) {
      // If the key exists, atomically increment the associated value
      auto ref =
        cuda::atomic_ref<typename Map::mapped_type, cuda::thread_scope_device>{found->second};
      ref.fetch_add(1, cuda::memory_order_relaxed);
    }
    tid += gridDim.x * blockDim.x;
  }
}

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

  // Create a sequence of keys and values {{0,0}, {1,1}, ... {i,i}}
  thrust::device_vector<Key> insert_keys(num_keys);
  thrust::sequence(insert_keys.begin(), insert_keys.end(), 0);
  thrust::device_vector<Value> insert_values(num_keys);
  thrust::sequence(insert_values.begin(), insert_values.end(), 0);

  // Compute capacity based on a 50% load factor
  auto constexpr load_factor = 0.5;
  std::size_t const capacity = std::ceil(num_keys / load_factor);

  // Constructs a map with "capacity" slots using -1 and -1 as the empty key/value sentinels.
  auto map = cuco::static_map{capacity,
                              cuco::empty_key{empty_key_sentinel},
                              cuco::empty_value{empty_value_sentinel},
                              cuda::std::equal_to<Key>{},
                              cuco::linear_probing<1, cuco::default_hash_function<Key>>{}};

  // Get a non-owning, mutable reference of the map that allows inserts to pass by value into the
  // kernel
  auto insert_ref = map.ref(cuco::insert);

  // Predicate will only insert even keys
  auto is_even = [] __device__(auto key) { return (key % 2) == 0; };

  // Allocate storage for count of number of inserted keys
  thrust::device_vector<int> num_inserted(1);

  auto constexpr block_size = 256;
  auto const grid_size      = (num_keys + block_size - 1) / block_size;
  filtered_insert<<<grid_size, block_size>>>(insert_ref,
                                             insert_keys.begin(),
                                             insert_values.begin(),
                                             num_keys,
                                             is_even,
                                             num_inserted.data().get());

  std::cout << "Number of keys inserted: " << num_inserted[0] << std::endl;

  // Get a non-owning reference of the map that allows find operations to pass by value into the
  // kernel
  auto find_ref = map.ref(cuco::find);

  increment_values<<<grid_size, block_size>>>(find_ref, insert_keys.begin(), num_keys);

  // Retrieve contents of all the non-empty slots in the map
  thrust::device_vector<Key> contained_keys(num_inserted[0]);
  thrust::device_vector<Value> contained_values(num_inserted[0]);
  map.retrieve_all(contained_keys.begin(), contained_values.begin());

  auto tuple_iter =
    thrust::make_zip_iterator(cuda::std::tuple{contained_keys.begin(), contained_values.begin()});
  // Iterate over all slot contents and verify that `slot.key + 1 == slot.value` is always true.
  auto result = thrust::all_of(
    thrust::device, tuple_iter, tuple_iter + num_inserted[0], [] __device__(auto const& tuple) {
      return cuda::std::get<0>(tuple) + 1 == cuda::std::get<1>(tuple);
    });

  if (result) { std::cout << "Success! Target values are properly incremented.\n"; }

  return 0;
}
