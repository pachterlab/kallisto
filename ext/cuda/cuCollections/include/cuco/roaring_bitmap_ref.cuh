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

#include <cuco/detail/roaring_bitmap/roaring_bitmap_impl.cuh>

#include <cuda/std/cstddef>
#include <cuda/stream_ref>

namespace cuco::experimental {

/**
 * @brief Non-owning reference to a Roaring bitmap stored in its serialized format.
 *
 * A `roaring_bitmap_ref` provides device and host APIs to query membership against a bitmap that
 * is laid out according to the [Roaring bitmap format
 * specification](https://github.com/RoaringBitmap/RoaringFormatSpec). The object does not own the
 * underlying storage; it simply provides algorithms over the referenced bytes.
 *
 * @note The reference reads directly from the serialized representation without deserializing.
 *       It supports 32-bit and 64-bit key types. For 32-bit bitmaps the layout follows the
 *       "Standard 32-bit Roaring Bitmap" format; for 64-bit bitmaps, the "portable" format is
 * supported.
 *
 * @tparam T Key type stored in the bitmap. Must be `cuda::std::uint32_t` or `cuda::std::uint64_t`.
 */
template <class T>
class roaring_bitmap_ref {
  using impl_type = detail::roaring_bitmap_impl<T>;

 public:
  using value_type       = T;                                     ///< Key type stored in the bitmap
  using storage_ref_type = typename impl_type::storage_ref_type;  ///< Implementation storage ref

  /**
   * @brief Constructs a non-owning reference from an implementation-specific storage reference.
   *
   * @param storage_ref Reference to the underlying serialized bitmap storage
   */
  __host__ __device__ roaring_bitmap_ref(storage_ref_type const& storage_ref);

  /**
   * @brief Constructs a device-side reference from a raw pointer to a 32-bit Roaring bitmap.
   *
   * @note This constructor is only available when `T == cuda::std::uint32_t` and can be used in
   *       device code to create a lightweight view over device-resident serialized bytes.
   *
   * @param bitmap Pointer to the beginning of the serialized bitmap in device memory
   */
  template <class U = T,
            class   = cuda::std::enable_if_t<cuda::std::is_same_v<U, cuda::std::uint32_t>>>
  __device__ roaring_bitmap_ref(cuda::std::byte const* bitmap);

  /**
   * @brief Bulk membership query for keys in `[first, last)`.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   *       `contains_async`.
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
  __host__ void contains(InputIt first,
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
  __host__ void contains_async(InputIt first,
                               InputIt last,
                               OutputIt contained,
                               cuda::stream_ref stream = cuda::stream_ref{
                                 cudaStream_t{nullptr}}) const noexcept;

  /**
   * @brief Device-side membership query for a single key.
   *
   * @param value Key to test for membership
   *
   * @return `true` iff `value` is contained in the bitmap
   */
  __device__ bool contains(T value) const;

  /**
   * @brief Number of keys stored in the bitmap.
   *
   * @return Count of keys in the bitmap
   */
  [[nodiscard]] __host__ __device__ cuda::std::size_t size() const noexcept;

  /**
   * @brief Checks whether the bitmap contains no keys.
   *
   * @return `true` iff `size() == 0`
   */
  [[nodiscard]] __host__ __device__ bool empty() const noexcept;

  /**
   * @brief Returns a pointer to the beginning of the serialized bitmap bytes.
   *
   * @return Pointer to the serialized storage
   */
  [[nodiscard]] __host__ __device__ cuda::std::byte const* data() const noexcept;

  /**
   * @brief Size in bytes of the serialized bitmap storage.
   *
   * @return Number of bytes occupied by the serialized bitmap
   */
  [[nodiscard]] __host__ __device__ cuda::std::size_t size_bytes() const noexcept;

 private:
  impl_type impl_;
};

}  // namespace cuco::experimental

#include <cuco/detail/roaring_bitmap/roaring_bitmap_ref.inl>