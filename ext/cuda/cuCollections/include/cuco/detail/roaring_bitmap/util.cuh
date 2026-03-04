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

#include <cuco/detail/error.hpp>
#include <cuco/utility/traits.hpp>

#include <cuda/std/cstddef>
#include <cuda/std/cstdint>
#include <cuda/std/iterator>
#include <cuda/std/memory>

#include <nv/target>
#include <vector>

namespace cuco::experimental::detail {

template <class T>
__host__ __device__ __forceinline__ T aligned_load(cuda::std::byte const* ptr)
{
  return *reinterpret_cast<T const*>(cuda::std::assume_aligned<alignof(T)>(ptr));
}

template <class T>
__host__ __device__ __forceinline__ T misaligned_load(cuda::std::byte const* ptr)
{
  T value;
  cuda::std::memcpy(&value, ptr, sizeof(T));
  return value;
}

__host__ __device__ __forceinline__ bool check_bit(cuda::std::byte const* bitmap,
                                                   cuda::std::uint32_t index)
{
  // check if the bit at index is set
  return static_cast<cuda::std::uint8_t>(bitmap[index / 8]) &
         (cuda::std::uint8_t(1) << (index % 8));
}

template <class T>
struct roaring_bitmap_metadata {
  static_assert(cuco::dependent_false<T>, "T must be either uint32_t or uint64_t");
};

/**
 * @brief Metadata structure for 32-bit roaring bitmap
 *
 * Contains metadata information for a 32-bit roaring bitmap including size, container information,
 * and validity status.
 */
template <>
struct roaring_bitmap_metadata<cuda::std::uint32_t> {
  /// Maximum number of elements in an array container before converting to bitmap
  static constexpr cuda::std::uint32_t max_array_container_card = 4096;

  /// Total size of the bitmap in bytes
  cuda::std::size_t size_bytes = 0;
  /// Number of keys/elements in the bitmap
  cuda::std::size_t num_keys = 0;
  /// Bitmap indicating which containers are run containers
  cuda::std::uint32_t run_container_bitmap = 0;
  /// Offset to key cardinality data
  cuda::std::uint32_t key_cards = 0;
  /// Offset to container offset data
  cuda::std::uint32_t container_offsets = 0;
  /// Number of containers in the bitmap
  cuda::std::int32_t num_containers = 0;
  /// Whether the bitmap contains run containers
  bool has_run = false;
  /// Whether the metadata is valid
  bool valid = false;

  /**
   * @brief Constructs metadata from a serialized bitmap
   *
   * @param bitmap Pointer to the beginning of the serialized bitmap
   */
  __host__ __device__ roaring_bitmap_metadata(cuda::std::byte const* bitmap)
  {
    constexpr cuda::std::uint32_t serial_cookie_no_runcontainer = 12346;
    constexpr cuda::std::uint32_t serial_cookie                 = 12347;
    // constexpr cuda::std::uint32_t frozen_cookie                 = 13766; // not implemented
    constexpr cuda::std::int32_t no_offset_threshold     = 4;
    constexpr cuda::std::int32_t max_containers          = 1 << 16;
    constexpr cuda::std::uint32_t cookie_mask            = 0xFFFF;
    constexpr cuda::std::uint32_t cookie_shift           = 16;
    constexpr cuda::std::uint32_t bitset_container_bytes = 8192;

    cuda::std::byte const* buf = bitmap;

    cuda::std::uint32_t cookie;
    cuda::std::memcpy(&cookie, buf, sizeof(cuda::std::uint32_t));
    buf += sizeof(cuda::std::uint32_t);
    if ((cookie & cookie_mask) != serial_cookie && cookie != serial_cookie_no_runcontainer) {
      valid = false;
      NV_IF_TARGET(
        NV_IS_HOST,
        CUCO_FAIL(
          "Invalid bitmap format: cookie type invalid or not supported");)  // TODO device error
                                                                            // handling
      return;
    }

    if ((cookie & cookie_mask) == serial_cookie)
      // upper 16 bits of cookie are the number of containers - 1
      num_containers = (cookie >> cookie_shift) + 1;
    else {
      // following 4 bytes are the number of containers
      cuda::std::memcpy(&num_containers, buf, sizeof(cuda::std::uint32_t));
      buf += sizeof(cuda::std::uint32_t);
    }
    if (num_containers < 0 or num_containers > max_containers) {
      valid = false;
      NV_IF_TARGET(
        NV_IS_HOST,
        CUCO_FAIL(
          "Invalid bitmap format: num_containers out of range");)  // TODO device error handling
      return;
    }

    has_run = (cookie & cookie_mask) == serial_cookie;
    if (has_run) {
      cuda::std::size_t s  = (num_containers + 7) / 8;  // ceil bytes to store run container bitmap
      run_container_bitmap = cuda::std::distance(bitmap, buf);
      buf += s;
    }

    key_cards = cuda::std::distance(bitmap, buf);
    // if the current address is aligned to 2 bytes, then all containers are aligned to at least 2
    // bytes
    bool const aligned_16 = (reinterpret_cast<cuda::std::uintptr_t>(bitmap + key_cards) %
                             sizeof(cuda::std::uint16_t)) == 0;
    buf += num_containers * 2 * sizeof(cuda::std::uint16_t);

    if ((!has_run) || (num_containers >= no_offset_threshold)) {
      container_offsets = cuda::std::distance(bitmap, buf);
      buf += num_containers * sizeof(cuda::std::uint32_t);
    } else {
      valid = false;
      NV_IF_TARGET(
        NV_IS_HOST,
        CUCO_FAIL("Invalid bitmap format: not implemented");)  // TODO device error handling
      return;
    }

    cuda::std::uint32_t card = 0;
    for (cuda::std::int32_t i = 0; i < num_containers; i++) {
      cuda::std::byte const* card_ptr =
        bitmap + key_cards + (i * 2 + 1) * sizeof(cuda::std::uint16_t);
      if (aligned_16) {
        card = 1u + aligned_load<cuda::std::uint16_t>(card_ptr);
      } else {
        card = 1u + misaligned_load<cuda::std::uint16_t>(card_ptr);
      }
      num_keys += card;
    }

    // find end of roaring bitmap (re-use card from last container)
    cuda::std::byte const* end =
      bitmap + misaligned_load<cuda::std::uint32_t>(
                 bitmap + container_offsets + (num_containers - 1) * sizeof(cuda::std::uint32_t));
    if (has_run and check_bit(bitmap + run_container_bitmap, num_containers - 1)) {
      cuda::std::uint16_t const num_runs = misaligned_load<cuda::std::uint16_t>(end);
      end += sizeof(cuda::std::uint16_t) + num_runs * 2 * sizeof(cuda::std::uint16_t);
    } else {
      if (card <= max_array_container_card) {
        end += card * sizeof(cuda::std::uint16_t);
      } else {
        end += bitset_container_bytes;  // fixed size bitset container
      }
    }

    size_bytes = static_cast<cuda::std::size_t>(cuda::std::distance(bitmap, end));
    valid      = true;
  }
};

/**
 * @brief Metadata structure for 64-bit roaring bitmap
 *
 * Contains metadata information for a 64-bit roaring bitmap including bucket information,
 * size, and validity status.
 */
template <>
struct roaring_bitmap_metadata<cuda::std::uint64_t> {
  /// Number of buckets in the 64-bit bitmap
  cuda::std::size_t num_buckets = 0;
  /// Total size of the bitmap in bytes
  cuda::std::size_t size_bytes = 0;
  /// Number of keys/elements in the bitmap
  cuda::std::size_t num_keys = 0;
  /// Whether the metadata is valid
  bool valid = false;

  /**
   * @brief Metadata for individual buckets in a 64-bit roaring bitmap
   *
   * Each bucket contains a 32-bit roaring bitmap with its own metadata.
   */
  struct bucket_metadata {
    /// Byte offset of this bucket in the serialized data
    cuda::std::size_t byte_offset;
    /// Key associated with this bucket (upper 32 bits)
    cuda::std::uint32_t key;
    /// Metadata for the 32-bit roaring bitmap in this bucket
    roaring_bitmap_metadata<cuda::std::uint32_t> metadata;

    /**
     * @brief Constructs bucket metadata
     *
     * @param offset Byte offset of the bucket
     * @param k Key associated with the bucket
     * @param meta Metadata for the bucket's roaring bitmap
     */
    bucket_metadata(cuda::std::size_t offset,
                    cuda::std::uint32_t k,
                    roaring_bitmap_metadata<cuda::std::uint32_t> const& meta)
      : byte_offset{offset}, key{k}, metadata{meta}
    {
    }
  };

  /**
   * @brief Constructs metadata from a serialized 64-bit bitmap with bucket metadata
   *
   * @param bitmap Pointer to the beginning of the serialized bitmap
   * @param bucket_metadata Vector to store metadata for each bucket
   */
  __host__ roaring_bitmap_metadata(cuda::std::byte const* bitmap,
                                   std::vector<bucket_metadata>& bucket_metadata)
  {
    cuda::std::size_t byte_offset     = 0;
    cuda::std::byte const* bitmap_ptr = bitmap;
    cuda::std::memcpy(&num_buckets, bitmap_ptr, sizeof(cuda::std::uint64_t));
    byte_offset += sizeof(cuda::std::uint64_t);  // skip num_buckets

    bucket_metadata.clear();
    bucket_metadata.reserve(num_buckets);

    for (cuda::std::size_t i = 0; i < num_buckets; ++i) {
      cuda::std::uint32_t bucket_key;
      cuda::std::memcpy(&bucket_key, bitmap_ptr + byte_offset, sizeof(cuda::std::uint32_t));
      byte_offset += sizeof(cuda::std::uint32_t);  // skip bucket key
      roaring_bitmap_metadata<cuda::std::uint32_t> bucket_meta{bitmap_ptr + byte_offset};
      if (!bucket_meta.valid) {
        valid = false;
        return;
      }
      bucket_metadata.emplace_back(byte_offset, bucket_key, bucket_meta);
      num_keys += bucket_meta.num_keys;
      byte_offset += bucket_meta.size_bytes;  // skip bucket
    }
    size_bytes = byte_offset;
    valid      = true;
  }

  /**
   * @brief Constructs metadata from a serialized 64-bit bitmap
   *
   * @param bitmap Pointer to the beginning of the serialized bitmap
   */
  __host__ __device__ roaring_bitmap_metadata(cuda::std::byte const* bitmap)
  {
    cuda::std::size_t byte_offset     = 0;
    cuda::std::byte const* bitmap_ptr = bitmap;
    cuda::std::memcpy(&num_buckets, bitmap_ptr, sizeof(cuda::std::uint64_t));
    byte_offset += sizeof(cuda::std::uint64_t);  // skip num_buckets

    for (cuda::std::size_t i = 0; i < num_buckets; ++i) {
      byte_offset += sizeof(cuda::std::uint32_t);  // skip bucket key
      roaring_bitmap_metadata<cuda::std::uint32_t> bucket_meta{bitmap_ptr + byte_offset};
      if (!bucket_meta.valid) {
        valid = false;
        return;
      }
      num_keys += bucket_meta.num_keys;
      byte_offset += bucket_meta.size_bytes;  // skip bucket
    }
    size_bytes = byte_offset;
    valid      = true;
  }
};
}  // namespace cuco::experimental::detail