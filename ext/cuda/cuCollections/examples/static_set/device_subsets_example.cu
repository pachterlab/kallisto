/*
 * Copyright (c) 2023-2025, NVIDIA CORPORATION.
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

#include <cuco/static_set_ref.cuh>

#include <cuda/std/array>
#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/reduce.h>
#include <thrust/scan.h>

#include <cooperative_groups.h>

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <numeric>

/**
 * @file device_subsets_example.cu
 * @brief Demonstrates how to use one bulk set storage to create multiple subsets and perform
 * individual operations via device-side ref APIs.
 *
 * To optimize memory usage, especially when dealing with expensive data allocation and multiple
 * hashsets, a practical solution involves employing a single bulk storage for generating subsets.
 * This eliminates the need for separate memory allocation and deallocation for each container. This
 * can be achieved by using the lightweight non-owning ref type.
 *
 * @note This example is for demonstration purposes only. It is not intended to show the most
 * performant way to do the example algorithm.
 */

auto constexpr cg_size     = 8;   ///< A CUDA Cooperative Group of 8 threads to handle each subset
auto constexpr bucket_size = 1;   ///< Number of concurrent slots handled by each thread
auto constexpr N           = 10;  ///< Number of elements to insert and query

using key_type = int;  ///< Key type
using probing_scheme_type =
  cuco::linear_probing<cg_size,
                       cuco::default_hash_function<key_type>>;  ///< Type controls CG granularity
                                                                ///< and probing scheme (linear
                                                                ///< probing v.s. double hashing)
/// Type of bulk allocation storage
using storage_type = cuco::bucket_storage<key_type, bucket_size>;
/// Lightweight non-owning storage ref type
using storage_ref_type = typename storage_type::ref_type;
using ref_type         = cuco::static_set_ref<key_type,
                                              cuda::thread_scope_device,
                                              cuda::std::equal_to<key_type>,
                                              probing_scheme_type,
                                              storage_ref_type>;  ///< Set ref type

/// Sample data to insert and query
__device__ constexpr cuda::std::array<key_type, N> data = {1, 3, 5, 7, 9, 11, 13, 15, 17, 19};
/// Empty slots are represented by reserved "sentinel" values. These values should be selected such
/// that they never occur in your input data.
key_type constexpr empty_key_sentinel = -1;

/**
 * @brief Inserts sample data into subsets by using cooperative group
 *
 * Each Cooperative Group creates its own subset and inserts `N` sample data.
 *
 * @param set_refs Pointer to the array of subset objects
 */
__global__ void insert(ref_type* set_refs)
{
  namespace cg = cooperative_groups;

  auto const tile = cg::tiled_partition<cg_size, cg::thread_block>(cg::this_thread_block());
  // Get subset (or CG) index
  auto const idx = (blockDim.x * blockIdx.x + threadIdx.x) / cg_size;

  auto raw_set_ref    = *(set_refs + idx);
  auto insert_set_ref = raw_set_ref.rebind_operators(cuco::insert);

  // Insert `N` elemtns into the set with CG insert
  for (int i = 0; i < N; i++) {
    insert_set_ref.insert(tile, data[i]);
  }
}

/**
 * @brief All inserted data can be found
 *
 * Each Cooperative Group reconstructs its own subset ref based on the storage parameters and
 * verifies all inserted data can be found.
 *
 * @param set_refs Pointer to the array of subset objects
 */
__global__ void find(ref_type* set_refs)
{
  namespace cg = cooperative_groups;

  auto const tile = cg::tiled_partition<cg_size, cg::thread_block>(cg::this_thread_block());
  auto const idx  = (blockDim.x * blockIdx.x + threadIdx.x) / cg_size;

  auto raw_set_ref  = *(set_refs + idx);
  auto find_set_ref = raw_set_ref.rebind_operators(cuco::find);

  // Result denoting if any of the inserted data is not found
  __shared__ int result;
  if (threadIdx.x == 0) { result = 0; }
  __syncthreads();

  for (int i = 0; i < N; i++) {
    // Query the set with inserted data
    auto const found = find_set_ref.find(tile, data[i]);
    // Record if the inserted data has been found
    atomicOr(&result, *found != data[i]);
  }
  __syncthreads();

  if (threadIdx.x == 0) {
    // If the result is still 0, all inserted data are found.
    if (result == 0) { printf("Success! Found all inserted elements.\n"); }
  }
}

int main()
{
  // Number of subsets to be created
  auto constexpr num = 16;
  // Each subset may have a different requested size
  auto constexpr subset_sizes =
    std::array<std::size_t, num>{20, 20, 20, 20, 30, 30, 30, 30, 40, 40, 40, 40, 50, 50, 50, 50};

  auto valid_sizes = std::vector<std::size_t>();
  valid_sizes.reserve(num);

  for (size_t i = 0; i < num; ++i) {
    valid_sizes.emplace_back(static_cast<std::size_t>(
      cuco::make_valid_extent<probing_scheme_type, storage_ref_type>(subset_sizes[i])));
  }

  std::vector<std::size_t> offsets(num + 1, 0);

  // prefix sum to compute offsets and total number of slots
  std::size_t current_sum = 0;
  for (std::size_t i = 0; i < valid_sizes.size(); ++i) {
    current_sum += valid_sizes[i];
    offsets[i + 1] = current_sum;
  }

  // total number of slots is located at the back of the offsets array
  auto const total_num_slots = offsets.back();

  // Create a single bulk storage used by all subsets
  auto set_storage = storage_type{
    total_num_slots, cuco::cuda_allocator<key_type>{}, cuda::stream_ref{cudaStream_t{nullptr}}};
  // Initializes the storage with the given sentinel
  set_storage.initialize(empty_key_sentinel);

  std::vector<ref_type> set_refs;

  // create subsets
  for (std::size_t i = 0; i < num; ++i) {
    storage_ref_type storage_ref{valid_sizes[i], set_storage.data() + offsets[i]};
    set_refs.emplace_back(
      ref_type{cuco::empty_key<key_type>{empty_key_sentinel}, {}, {}, {}, storage_ref});
  }

  thrust::device_vector<ref_type> d_set_refs(set_refs);

  // Insert sample data
  insert<<<1, 128>>>(d_set_refs.data().get());
  // Find all inserted data
  find<<<1, 128>>>(d_set_refs.data().get());

  return 0;
}
