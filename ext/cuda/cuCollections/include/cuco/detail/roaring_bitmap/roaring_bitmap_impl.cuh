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
#include <cuco/detail/roaring_bitmap/roaring_bitmap_storage.cuh>
#include <cuco/detail/roaring_bitmap/util.cuh>
#include <cuco/utility/traits.hpp>

#include <cub/device/device_transform.cuh>
#include <cuda/std/cstddef>
#include <cuda/std/cstdint>
#include <cuda/std/functional>
#include <cuda/std/iterator>
#include <cuda/stream_ref>
#include <thrust/iterator/constant_iterator.h>

namespace cuco::experimental::detail {

// primary template
template <class T>
class roaring_bitmap_impl {
  static_assert(cuco::dependent_false<T>, "T must be either uint32_t or uint64_t");
};

template <>
class roaring_bitmap_impl<cuda::std::uint32_t> {
 public:
  using storage_ref_type = roaring_bitmap_storage_ref<cuda::std::uint32_t>;

  static constexpr cuda::std::int32_t binary_search_threshold = 8;  // TODO determine optimal value

  __host__ __device__ roaring_bitmap_impl(storage_ref_type const& storage_ref)
    : storage_ref_{storage_ref},
      offsets_aligned_{(reinterpret_cast<cuda::std::uintptr_t>(
                         storage_ref_.data() + storage_ref_.metadata().container_offsets)) %
                         sizeof(cuda::std::uint32_t) ==
                       0},
      aligned_16_{(reinterpret_cast<cuda::std::uintptr_t>(storage_ref_.data() +
                                                          storage_ref_.metadata().key_cards)) %
                    sizeof(cuda::std::uint16_t) ==
                  0}  // if base address of key_cards is aligned, then all containers are aligned
  {
  }

  template <class InputIt, class OutputIt>
  __host__ void contains(InputIt first,
                         InputIt last,
                         OutputIt contained,
                         cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const
  {
    this->contains_async(first, last, contained, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
    stream.sync();
#else
    stream.wait();
#endif
  }

  template <class InputIt, class OutputIt>
  __host__ void contains_async(InputIt first,
                               InputIt last,
                               OutputIt contained,
                               cuda::stream_ref stream = cuda::stream_ref{
                                 cudaStream_t{nullptr}}) const noexcept
  {
    if (this->empty()) {
      cub::DeviceTransform::Transform(
        thrust::constant_iterator<bool>(false),
        contained,
        cuda::std::distance(first, last),
        cuda::proclaim_return_type<bool>([] __device__(auto /* dummy */) { return false; }),
        stream.get());
    } else {
      cub::DeviceTransform::Transform(
        first,
        contained,
        cuda::std::distance(first, last),
        cuda::proclaim_return_type<bool>(
          [*this] __device__(auto key) { return this->contains(key); }),
        stream.get());
    }
  }

  __device__ bool contains(cuda::std::uint32_t value) const
  {
    if (storage_ref_.metadata().num_keys == 0) { return false; }

    if (aligned_16_) {
      return this->dispatch_contains<true>(value);
    } else {
      return this->dispatch_contains<false>(value);
    }
  }

  template <bool Aligned>
  __device__ bool dispatch_contains(cuda::std::uint32_t value) const
  {
    cuda::std::uint16_t const upper = value >> 16;
    cuda::std::uint16_t const lower = value & 0xFFFF;
    cuda::std::uint16_t key;

    if (storage_ref_.metadata().num_containers < binary_search_threshold) {
// linear search
#pragma unroll
      for (cuda::std::int32_t i = 0; i < storage_ref_.metadata().num_containers; i++) {
        cuda::std::byte const* key_ptr =
          storage_ref_.key_cards() + (i * 2) * sizeof(cuda::std::uint16_t);
        if constexpr (Aligned) {
          key = aligned_load<cuda::std::uint16_t>(key_ptr);
        } else {
          key = misaligned_load<cuda::std::uint16_t>(key_ptr);
        }
        if (key == upper) { return this->contains_container<Aligned>(lower, i); }
        if (key > upper) { return false; }
      }
    } else {
      // binary search
      cuda::std::uint32_t left  = 0;
      cuda::std::uint32_t right = storage_ref_.metadata().num_containers;
      while (left < right) {
        cuda::std::uint32_t mid = left + (right - left) / 2;
        cuda::std::byte const* key_ptr =
          storage_ref_.key_cards() + (mid * 2) * sizeof(cuda::std::uint16_t);
        if constexpr (Aligned) {
          key = aligned_load<cuda::std::uint16_t>(key_ptr);
        } else {
          key = misaligned_load<cuda::std::uint16_t>(key_ptr);
        }

        if (key == upper) {
          return this->contains_container<Aligned>(lower, mid);
        } else if (key < upper) {
          left = mid + 1;
        } else {
          right = mid;
        }
      }
    }
    return false;
  }

  [[nodiscard]] __host__ __device__ cuda::std::size_t size() const noexcept
  {
    return storage_ref_.metadata().num_keys;
  }

  [[nodiscard]] __host__ __device__ bool empty() const noexcept { return this->size() == 0; }

  [[nodiscard]] __host__ __device__ cuda::std::byte const* data() const noexcept
  {
    return storage_ref_.data();
  }

  [[nodiscard]] __host__ __device__ cuda::std::size_t size_bytes() const noexcept
  {
    return storage_ref_.metadata().size_bytes;
  }

  template <bool Aligned>
  __device__ bool contains_container(cuda::std::uint16_t lower, cuda::std::uint32_t index) const
  {
    cuda::std::uint32_t offset;
    cuda::std::byte const* offset_ptr =
      storage_ref_.container_offsets() + index * sizeof(cuda::std::uint32_t);
    if (offsets_aligned_) {
      offset = aligned_load<cuda::std::uint32_t>(offset_ptr);
    } else {
      offset = misaligned_load<cuda::std::uint32_t>(offset_ptr);
    }
    cuda::std::byte const* container = storage_ref_.data() + offset;
    if (storage_ref_.metadata().has_run and check_bit(storage_ref_.run_container_bitmap(), index)) {
      return this->contains_run_container<Aligned>(container, lower);
    } else {
      cuda::std::uint32_t card;
      cuda::std::byte const* card_ptr =
        storage_ref_.key_cards() + (index * 2 + 1) * sizeof(cuda::std::uint16_t);
      if constexpr (Aligned) {
        card = 1u + aligned_load<cuda::std::uint16_t>(card_ptr);
      } else {
        card = 1u + misaligned_load<cuda::std::uint16_t>(card_ptr);
      }
      if (card <= storage_ref_type::metadata_type::max_array_container_card) {
        return this->contains_array_container<Aligned>(container, lower, card);
      } else {
        return this->contains_bitset_container(container, lower);
      }
    }
  }

  template <bool Aligned>
  __device__ bool contains_array_container(cuda::std::byte const* container,
                                           cuda::std::uint16_t lower,
                                           cuda::std::uint32_t card) const
  {
    cuda::std::uint16_t elem;
    // Use linear search for small arrays, binary search for larger ones
    if (card < binary_search_threshold) {
      for (cuda::std::uint32_t i = 0; i < card; i++) {
        cuda::std::byte const* elem_ptr = container + i * sizeof(cuda::std::uint16_t);
        if constexpr (Aligned) {
          elem = aligned_load<cuda::std::uint16_t>(elem_ptr);
        } else {
          elem = misaligned_load<cuda::std::uint16_t>(elem_ptr);
        }
        if (elem == lower) { return true; }
      }
      return false;
    } else {
      cuda::std::uint32_t left  = 0;
      cuda::std::uint32_t right = card;

      while (left < right) {
        cuda::std::uint32_t mid         = left + (right - left) / 2;
        cuda::std::byte const* elem_ptr = container + mid * sizeof(cuda::std::uint16_t);
        if constexpr (Aligned) {
          elem = aligned_load<cuda::std::uint16_t>(elem_ptr);
        } else {
          elem = misaligned_load<cuda::std::uint16_t>(elem_ptr);
        }
        if (elem == lower) {
          return true;
        } else if (elem < lower) {
          left = mid + 1;
        } else {
          right = mid;
        }
      }
      return false;
    }
  }

  __device__ bool contains_bitset_container(cuda::std::byte const* container,
                                            cuda::std::uint16_t lower) const
  {
    return check_bit(container, lower);
  }

  template <bool Aligned>
  __device__ bool contains_run_container(cuda::std::byte const* container,
                                         cuda::std::uint16_t lower) const
  {
    // TODO implement binary search
    cuda::std::uint16_t num_runs;
    if constexpr (Aligned) {
      num_runs = aligned_load<cuda::std::uint16_t>(container);
    } else {
      num_runs = misaligned_load<cuda::std::uint16_t>(container);
    }

    cuda::std::uint16_t start;
    cuda::std::uint32_t end;

    for (cuda::std::uint32_t i = 0; i < num_runs; i++) {
      // the first 16 bits of the run container denotes the number of runs
      // followed by the sequence of runs as (start, end) U16 pairs
      cuda::std::byte const* start_ptr = container + (i * 2 + 1) * sizeof(cuda::std::uint16_t);
      // TODO load start+end in one instruction
      if constexpr (Aligned) {
        start = aligned_load<cuda::std::uint16_t>(start_ptr);
        end   = static_cast<cuda::std::uint32_t>(start) +
              aligned_load<cuda::std::uint16_t>(start_ptr + sizeof(cuda::std::uint16_t));
      } else {
        start = misaligned_load<cuda::std::uint16_t>(start_ptr);
        end   = static_cast<cuda::std::uint32_t>(start) +
              misaligned_load<cuda::std::uint16_t>(start_ptr + sizeof(cuda::std::uint16_t));
      }
      if (start <= lower && end >= lower) { return true; }
      if (start > lower) { break; }
    }
    return false;
  }

  storage_ref_type storage_ref_;
  bool offsets_aligned_;
  bool aligned_16_;
};

template <>
class roaring_bitmap_impl<cuda::std::uint64_t> {
 public:
  using bucket_type      = roaring_bitmap_impl<cuda::std::uint32_t>;
  using storage_ref_type = roaring_bitmap_storage_ref<cuda::std::uint64_t>;

  __host__ __device__ roaring_bitmap_impl(storage_ref_type const& storage_ref)
    : storage_ref_{storage_ref}
  {
  }

  template <class InputIt, class OutputIt>
  __host__ void contains(InputIt first,
                         InputIt last,
                         OutputIt contained,
                         cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const
  {
    this->contains_async(first, last, contained, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
    stream.sync();
#else
    stream.wait();
#endif
  }

  template <class InputIt, class OutputIt>
  __host__ void contains_async(InputIt first,
                               InputIt last,
                               OutputIt contained,
                               cuda::stream_ref stream = cuda::stream_ref{
                                 cudaStream_t{nullptr}}) const noexcept
  {
    if (this->empty()) {
      cub::DeviceTransform::Transform(
        thrust::constant_iterator<bool>(false),
        contained,
        cuda::std::distance(first, last),
        cuda::proclaim_return_type<bool>([] __device__(auto /* dummy */) { return false; }),
        stream.get());
    } else {
      cub::DeviceTransform::Transform(
        first,
        contained,
        cuda::std::distance(first, last),
        cuda::proclaim_return_type<bool>(
          [*this] __device__(auto key) { return this->contains(key); }),
        stream.get());
    }
  }

  __device__ bool contains(cuda::std::uint64_t value) const
  {
    cuda::std::uint32_t bucket_key   = value >> 32;
    cuda::std::uint32_t bucket_value = value & 0xFFFFFFFF;

    // binary search in storage_ref_.buckets()
    cuda::std::uint32_t left  = 0;
    cuda::std::uint32_t right = storage_ref_.metadata().num_buckets;
    while (left < right) {
      cuda::std::uint32_t mid = left + (right - left) / 2;
      if (storage_ref_.buckets()[mid].first == bucket_key) {
        return bucket_type{storage_ref_.buckets()[mid].second}.contains(
          bucket_value);  // TODO is constructing the ref in-place a bad idea?
      } else if (storage_ref_.buckets()[mid].first < bucket_key) {
        left = mid + 1;
      } else {
        right = mid;
      }
    }
    return false;
  }

  [[nodiscard]] __host__ __device__ cuda::std::size_t size() const noexcept
  {
    return storage_ref_.metadata().num_keys;
  }

  [[nodiscard]] __host__ __device__ bool empty() const noexcept { return this->size() == 0; }

  [[nodiscard]] __host__ __device__ cuda::std::byte const* data() const noexcept
  {
    return storage_ref_.data();
  }

  [[nodiscard]] __host__ __device__ cuda::std::size_t size_bytes() const noexcept
  {
    return storage_ref_.metadata().size_bytes;
  }

  storage_ref_type storage_ref_;
};

}  // namespace cuco::experimental::detail