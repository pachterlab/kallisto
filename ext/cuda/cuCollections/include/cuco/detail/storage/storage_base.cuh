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

#pragma once

#include <cuda/stream_ref>

#include <cstddef>

namespace cuco {
namespace detail {
/**
 * @brief Custom deleter for unique pointer.
 *
 * @tparam SizeType Type of device storage size
 * @tparam Allocator Type of allocator used for device storage
 */
template <typename SizeType, typename Allocator>
struct custom_deleter {
  using pointer = typename Allocator::value_type*;  ///< Value pointer type

  /**
   * @brief Constructor of custom deleter.
   *
   * @param size Number of values to deallocate
   * @param allocator Allocator used for deallocating device storage
   * @param stream Stream to use for deallocation
   */
  explicit constexpr custom_deleter(SizeType size, Allocator& allocator, cuda::stream_ref stream)
    : size_{size}, allocator_{allocator}, stream_{stream}
  {
  }

  /**
   * @brief Operator for deallocation
   *
   * @param ptr Pointer to the first value for deallocation
   */
  void operator()(pointer ptr) { allocator_.deallocate(ptr, size_, stream_); }

  SizeType size_;            ///< Number of values to delete
  Allocator& allocator_;     ///< Allocator used deallocating values
  cuda::stream_ref stream_;  ///< Stream used for deallocation
};

/**
 * @brief Base class of open addressing storage.
 *
 * This class should not be used directly.
 *
 * @tparam Extent Type of extent denoting storage capacity
 */
template <typename Extent>
class storage_base {
 public:
  using extent_type = Extent;                            ///< Storage extent type
  using size_type   = typename extent_type::value_type;  ///< Storage size type

  /**
   * @brief Constructor of base storage.
   *
   * @param size Number of elements to (de)allocate
   */
  __host__ __device__ explicit constexpr storage_base(Extent size) : extent_{size} {}

  /**
   * @brief Gets the total number of elements in the current storage.
   *
   * @return The total number of elements
   */
  [[nodiscard]] __host__ __device__ constexpr size_type capacity() const noexcept
  {
    return static_cast<size_type>(extent_);
  }

  /**
   * @brief Gets the extent of the current storage.
   *
   * @return The extent.
   */
  [[nodiscard]] __host__ __device__ constexpr extent_type extent() const noexcept
  {
    return extent_;
  }

 protected:
  extent_type extent_;  ///< Total number of elements
};

}  // namespace detail
}  // namespace cuco
