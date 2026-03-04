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
 */

#pragma once

#include <cuda/std/cstdint>

#include <cooperative_groups.h>

#if defined(CUCO_DISABLE_KERNEL_VISIBILITY_WARNING_SUPPRESSION)
#define CUCO_SUPPRESS_KERNEL_WARNINGS
#elif defined(__NVCC__) && (defined(__GNUC__) || defined(__clang__))
// handle when nvcc is the CUDA compiler and gcc or clang is host
#define CUCO_SUPPRESS_KERNEL_WARNINGS _Pragma("nv_diag_suppress 1407")
_Pragma("GCC diagnostic ignored \"-Wattributes\"")
#elif defined(__clang__)
// handle when clang is the CUDA compiler
#define CUCO_SUPPRESS_KERNEL_WARNINGS _Pragma("clang diagnostic ignored \"-Wattributes\"")
#elif defined(__NVCOMPILER)
#define CUCO_SUPPRESS_KERNEL_WARNINGS #pragma diag_suppress attribute_requires_external_linkage
#endif

#ifndef CUCO_KERNEL
#define CUCO_KERNEL __attribute__((visibility("hidden"))) __global__
#endif
namespace cuco {
namespace detail {

using index_type = cuda::std::int64_t;  ///< CUDA thread index type

/// Default block size
/// CUDA warp size
[[nodiscard]] __device__ constexpr cuda::std::int32_t warp_size() noexcept { return 32; }

/**
 * @brief Returns the global thread index in a 1D scalar grid
 *
 * @return The global thread index
 */
[[nodiscard]] __device__ inline index_type global_thread_id() noexcept
{
  return index_type{threadIdx.x} + index_type{blockDim.x} * index_type{blockIdx.x};
}

/**
 * @brief Returns the grid stride of a 1D grid
 *
 * @return The grid stride
 */
[[nodiscard]] __device__ inline index_type grid_stride() noexcept
{
  return index_type{gridDim.x} * index_type{blockDim.x};
}

/**
 * @brief Constexpr helper to extract the size of a Cooperative Group.
 *
 * @tparam Tile The Cooperative Group type
 */
template <typename Tile>
struct tile_size;

/**
 * @brief Specialization of `cuco::detail::tile_size` for 'cooperative_groups::thread_block_tile'.
 *
 * @tparam CGSize The Cooperative Group size
 * @tparam ParentCG The Cooperative Group the tile has been created from
 */
template <uint32_t CGSize, class ParentCG>
struct tile_size<cooperative_groups::thread_block_tile<CGSize, ParentCG>> {
  static constexpr uint32_t value = CGSize;  ///< Size of the `thread_block_tile`
};

template <typename Tile>
__device__ constexpr uint32_t tile_size_v = tile_size<Tile>::value;

}  // namespace detail
}  // namespace cuco
