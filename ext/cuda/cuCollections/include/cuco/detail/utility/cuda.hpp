/*
 * Copyright (c) 2021-2025, NVIDIA CORPORATION.
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
 */

#pragma once

#include <cuco/detail/error.hpp>
#include <cuco/detail/utility/cuda.cuh>

#include <cstdint>

namespace cuco {
namespace detail {

constexpr std::int32_t default_block_size() noexcept { return 128; }
/// Default stride
constexpr std::int32_t default_stride() noexcept { return 1; }

/**
 * @brief Computes the desired 1D grid size with the given parameters
 *
 * @param num Number of elements to handle in the kernel
 * @param cg_size Number of threads per CUDA Cooperative Group
 * @param stride Number of elements to be handled by each thread
 * @param block_size Number of threads in each thread block
 *
 * @return The resulting grid size
 */
constexpr auto grid_size(index_type num,
                         std::int32_t cg_size    = 1,
                         std::int32_t stride     = default_stride(),
                         std::int32_t block_size = default_block_size()) noexcept
{
  return (cg_size * num + stride * block_size - 1) / (stride * block_size);
}

/**
 * @brief Computes the ideal 1D grid size with the given parameters
 *
 * @tparam Kernel Kernel type
 *
 * @param block_size Number of threads in each thread block
 * @param kernel CUDA kernel to launch
 * @param dynamic_shm_size Dynamic shared memory size
 *
 * @return The grid size that delivers the highest occupancy
 */
template <typename Kernel>
constexpr auto max_occupancy_grid_size(std::int32_t block_size,
                                       Kernel kernel,
                                       std::size_t dynamic_shm_size = 0)
{
  int device = 0;
  CUCO_CUDA_TRY(cudaGetDevice(&device));

  int num_multiprocessors = -1;
  CUCO_CUDA_TRY(
    cudaDeviceGetAttribute(&num_multiprocessors, cudaDevAttrMultiProcessorCount, device));

  int max_active_blocks_per_multiprocessor{};
  CUCO_CUDA_TRY(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
    &max_active_blocks_per_multiprocessor, kernel, block_size, dynamic_shm_size));

  return max_active_blocks_per_multiprocessor * num_multiprocessors;
}

}  // namespace detail
}  // namespace cuco
