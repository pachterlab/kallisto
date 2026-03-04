/*
 * Copyright (c) 2025 NVIDIA CORPORATION.
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

#include <cuco/detail/roaring_bitmap/roaring_bitmap_storage.cuh>
#include <cuco/roaring_bitmap_ref.cuh>
#include <cuco/utility/allocator.hpp>

#include <cuda/std/cstddef>
#include <cuda/stream_ref>

namespace cuco::experimental {

/**
 * @brief GPU-accelerated container that owns a serialized Roaring bitmap.
 *
 * The `roaring_bitmap` provides host-side bulk membership queries over a bitmap stored in the
 * [Roaring bitmap format specification](https://github.com/RoaringBitmap/RoaringFormatSpec).
 * The serialized bytes are copied to device-accessible storage upon construction, and queries are
 * executed on the GPU.
 *
 * In addition to bulk host APIs such as `contains`/`contains_async`, this container exposes a
 * non-owning reference object via `ref()` that can be used for device-side per-thread queries.
 *
 * @tparam T Key type. Must be `cuda::std::uint32_t` or `cuda::std::uint64_t`.
 * @tparam Allocator Allocator type used to manage device-accessible storage for the serialized
 *                   bytes.
 */
template <class T, class Allocator = cuco::cuda_allocator<cuda::std::byte>>
class roaring_bitmap {
 public:
  using value_type     = T;                                             ///< Key type
  using storage_type   = detail::roaring_bitmap_storage<T, Allocator>;  ///< Storage implementation
  using allocator_type = typename storage_type::allocator_type;         ///< Allocator type
  using ref_type       = roaring_bitmap_ref<value_type>;  ///< Non-owning reference type

  /**
   * @brief Constructs a `roaring_bitmap` by copying the serialized bytes to device-accessible
   *        storage.
   *
   * @param bitmap Pointer to the beginning of the serialized bitmap in host memory
   * @param alloc Allocator used to allocate device-accessible storage
   * @param stream CUDA stream used for device memory operations during construction
   */
  roaring_bitmap(cuda::std::byte const* bitmap,
                 Allocator const& alloc  = {},
                 cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Copy constructor
   *
   * @param other The roaring_bitmap to copy from
   */
  roaring_bitmap(roaring_bitmap const& other) = default;

  /**
   * @brief Move constructor
   *
   * @param other The roaring_bitmap to move from
   */
  roaring_bitmap(roaring_bitmap&& other) = default;

  /**
   * @brief Copy assignment operator
   *
   * @param other The roaring_bitmap to copy from
   * @return Reference to this roaring_bitmap
   */
  roaring_bitmap& operator=(roaring_bitmap const& other) = default;

  /**
   * @brief Move assignment operator
   *
   * @param other The roaring_bitmap to move from
   * @return Reference to this roaring_bitmap
   */
  roaring_bitmap& operator=(roaring_bitmap&& other) = default;

  ~roaring_bitmap() = default;  ///< Destructor

  /**
   * @brief Bulk membership query for keys in `[first, last)`.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   *       `contains_async`.
   *
   * @tparam InputIt  Device-accessible random access input iterator of keys convertible to `T`
   * @tparam OutputIt Device-accessible random access output iterator whose `value_type` is
   * constructible from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last  End of the sequence of keys
   * @param contained Output iterator where results are written; `true` iff the corresponding key
   *                  is present in the bitmap
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt, class OutputIt>
  void contains(InputIt first,
                InputIt last,
                OutputIt contained,
                cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Asynchronously performs a bulk membership query for keys in `[first, last)`.
   *
   * @tparam InputIt  Device-accessible random access input iterator of keys convertible to `T`
   * @tparam OutputIt Device-accessible random access output iterator to `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last  End of the sequence of keys
   * @param contained Output iterator where results are written; `true` iff the corresponding key
   *                  is present in the bitmap
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt, class OutputIt>
  void contains_async(InputIt first,
                      InputIt last,
                      OutputIt contained,
                      cuda::stream_ref stream = cuda::stream_ref{
                        cudaStream_t{nullptr}}) const noexcept;

  /**
   * @brief Number of keys stored in the bitmap.
   *
   * @return Count of keys in the bitmap
   */
  [[nodiscard]] cuda::std::size_t size() const noexcept;

  /**
   * @brief Checks whether the bitmap contains no keys.
   *
   * @return `true` iff `size() == 0`
   */
  [[nodiscard]] bool empty() const noexcept;

  /**
   * @brief Returns a pointer to the beginning of the serialized bitmap bytes in device-accessible
   *        storage.
   *
   * @return Pointer to the serialized storage
   */
  [[nodiscard]] cuda::std::byte const* data() const noexcept;

  /**
   * @brief Size in bytes of the serialized bitmap storage.
   *
   * @return Number of bytes occupied by the serialized bitmap
   */
  [[nodiscard]] cuda::std::size_t size_bytes() const noexcept;

  /**
   * @brief Returns the allocator used to manage device-accessible storage.
   *
   * @return Allocator instance
   */
  [[nodiscard]] allocator_type allocator() const noexcept;

  /**
   * @brief Returns a non-owning reference to the underlying bitmap suitable for device-side use.
   *
   * The returned reference type provides device functions such as `contains(T)` for per-thread
   * membership testing.
   *
   * @return Non-owning reference to the underlying bitmap
   */
  [[nodiscard]] ref_type ref() const noexcept;

 private:
  storage_type storage_;  ///< Storage type
};

}  // namespace cuco::experimental

#include <cuco/detail/roaring_bitmap/roaring_bitmap.inl>