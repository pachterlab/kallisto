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

#include <cuco/static_map.cuh>

#include <cuda/std/array>
#include <cuda/std/limits>

#include <cooperative_groups.h>

#include <cstdio>
#include <iostream>

/**
 * @file preallocated_memory_example.cu
 * @brief Demonstrates usage of static_map with pre-allocated device memory.
 *
 * This example shows how to use a static_map with device memory that is allocated
 * at compile time. This can be useful for cases where you want to avoid dynamic memory
 * allocation or when working with memory-constrained environments.
 *
 * @note This example is for demonstration purposes only. It is not intended to show the most
 * performant way to do the example algorithm.
 */

// Basic types
using Value = uint32_t;
using Key   = int;

// Sentinel values for empty slots
Key constexpr empty_key_sentinel     = -1;
Value constexpr empty_value_sentinel = cuda::std::numeric_limits<Value>::min();

// Map configuration
std::size_t constexpr capacity = 100'000;

// Type aliases for cleaner code
using probing_scheme_type = cuco::linear_probing<1, cuco::default_hash_function<Key>>;
using storage_type        = cuco::bucket_storage<cuco::pair<Key, Value>, 1>;
using storage_ref_type    = typename storage_type::ref_type;
using map_ref_type        = cuco::static_map_ref<Key,
                                                 Value,
                                                 cuda::thread_scope_device,
                                                 cuda::std::equal_to<Key>,
                                                 probing_scheme_type,
                                                 storage_ref_type>;

// Pre-allocated device memory
__device__ auto constexpr valid_extent =
  cuco::make_valid_extent<probing_scheme_type, cuco::storage<1>>(cuco::extent<size_t, capacity>{});

__device__ cuda::std::array<typename storage_ref_type::value_type, valid_extent.value()>
  storage_array;
__device__ int found_count = 0;  // Track number of items found (for verification)

// Helper function to create map reference from pre-allocated storage
__device__ map_ref_type create_map()
{
  storage_ref_type storage_ref{valid_extent.value(), storage_array.data()};
  return map_ref_type(cuco::empty_key{empty_key_sentinel},
                      cuco::empty_value{empty_value_sentinel},
                      cuda::std::equal_to<Key>{},
                      probing_scheme_type{},
                      cuco::cuda_thread_scope<cuda::thread_scope_device>{},
                      storage_ref);
}

/**
 * @brief Initialize the pre-allocated storage
 */
__global__ void init_kernel()
{
  auto map         = create_map();
  auto const block = cooperative_groups::this_thread_block();

  // Initialize storage using all threads in the block
  map.initialize(block);
  block.sync();
}

/**
 * @brief Insert key-value pairs (key -> key*2)
 */
__global__ void insert_kernel()
{
  auto map        = create_map();
  auto insert_ref = map.rebind_operators(cuco::insert);

  // Each thread inserts one key-value pair
  auto key = threadIdx.x + blockIdx.x;
  insert_ref.insert(cuco::pair(key, key * 2));
}

/**
 * @brief Find and verify inserted key-value pairs
 */
__global__ void find_kernel()
{
  auto map      = create_map();
  auto find_ref = map.rebind_operators(cuco::find);

  // Each thread looks up one key
  auto key    = threadIdx.x + blockIdx.x;
  auto result = find_ref.find(key);

  // Count successful finds (no printf output)
  if (result != find_ref.end()) {
    auto ref = cuda::atomic_ref<int, cuda::thread_scope_device>{found_count};
    ref.fetch_add(1, cuda::memory_order_relaxed);
  }
}

int main()
{
  // Step 1: Initialize the pre-allocated storage
  init_kernel<<<1, 128>>>();
  CUCO_CUDA_TRY(cudaDeviceSynchronize());

  // Step 2: Insert some key-value pairs
  insert_kernel<<<2, 32>>>();
  CUCO_CUDA_TRY(cudaDeviceSynchronize());

  // Step 3: Find and verify the inserted pairs
  find_kernel<<<2, 32>>>();

  // Check results - expect to find all 64 keys (2 blocks * 32 threads)
  int host_found_count;
  CUCO_CUDA_TRY(cudaMemcpyFromSymbol(&host_found_count, found_count, sizeof(int)));
  int expected_count = 2 * 32;  // Total number of keys inserted and queried

  if (host_found_count == expected_count) {
    std::cout << "Success: all keys are found" << std::endl;
  } else {
    std::cout << "Fail" << std::endl;
  }

  return 0;
}