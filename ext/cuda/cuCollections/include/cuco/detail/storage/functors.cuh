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
#pragma once

namespace cuco::detail {
/**
 * @brief Functor for initializing device memory with a given value
 *
 * @tparam SizeType Type used for indexing
 * @tparam T Type of value being initialized
 */
template <typename SizeType, typename T>
struct initialize_functor {
  T* const _d_ptr;  ///< Pointer to device memory
  T const _key;     ///< Value to initialize memory with

  /**
   * @brief Constructs functor for initializing device memory
   *
   * @param d_ptr Pointer to device memory to initialize
   * @param key Value to initialize memory with
   */
  __host__ __device__ initialize_functor(T* d_ptr, T key) noexcept : _d_ptr{d_ptr}, _key{key} {}

  /**
   * @brief Device function to initialize memory at given index
   *
   * @param idx Index into device memory
   */
  __device__ __forceinline__ void operator()(SizeType idx) const noexcept { _d_ptr[idx] = _key; }
};
}  // namespace cuco::detail
