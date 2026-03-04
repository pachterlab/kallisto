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

#include <cuco/detail/storage/functors.cuh>
#include <cuco/detail/utility/cuda.hpp>
#include <cuco/extent.cuh>

#include <cub/device/device_for.cuh>
#include <cuda/std/array>
#include <cuda/std/bit>
#include <cuda/stream_ref>

#include <cassert>

namespace cuco {

template <typename T, int BucketSize, typename Extent>
__host__ __device__ constexpr bucket_storage_ref<T, BucketSize, Extent>::bucket_storage_ref(
  Extent size, value_type* slots) noexcept
  : extent_{size}, slots_{slots}
{
  static_assert(cuda::std::has_single_bit(alignment), "Alignment must be a power of 2");
  assert(reinterpret_cast<cuda::std::uintptr_t>(slots) % alignment == 0 &&
         "Storage must be properly aligned");
}

template <typename T, int BucketSize, typename Extent>
__host__ __device__ constexpr bucket_storage_ref<T, BucketSize, Extent>::iterator
bucket_storage_ref<T, BucketSize, Extent>::end() noexcept
{
  return iterator{reinterpret_cast<value_type*>(this->data() + this->capacity())};
}

template <typename T, int BucketSize, typename Extent>
__host__ __device__ constexpr bucket_storage_ref<T, BucketSize, Extent>::iterator
bucket_storage_ref<T, BucketSize, Extent>::end() const noexcept
{
  return iterator{reinterpret_cast<value_type*>(this->data() + this->capacity())};
}

template <typename T, int BucketSize, typename Extent>
__host__ __device__ constexpr bucket_storage_ref<T, BucketSize, Extent>::value_type*
bucket_storage_ref<T, BucketSize, Extent>::data() noexcept
{
  return slots_;
}

template <typename T, int BucketSize, typename Extent>
__host__ __device__ constexpr bucket_storage_ref<T, BucketSize, Extent>::value_type*
bucket_storage_ref<T, BucketSize, Extent>::data() const noexcept
{
  return slots_;
}

template <typename T, int BucketSize, typename Extent>
__device__ constexpr bucket_storage_ref<T, BucketSize, Extent>::bucket_type
bucket_storage_ref<T, BucketSize, Extent>::operator[](size_type index) const noexcept
{
  return *reinterpret_cast<bucket_type*>(this->data() + index);
}

template <typename T, int BucketSize, typename Extent>
__host__ __device__ constexpr typename bucket_storage_ref<T, BucketSize, Extent>::size_type
bucket_storage_ref<T, BucketSize, Extent>::num_buckets() const noexcept
{
  return static_cast<size_type>(extent_) / bucket_size;
}

template <typename T, int BucketSize, typename Extent>
__host__ __device__ constexpr typename bucket_storage_ref<T, BucketSize, Extent>::size_type
bucket_storage_ref<T, BucketSize, Extent>::capacity() const noexcept
{
  return static_cast<size_type>(extent_);
}

template <typename T, int BucketSize, typename Extent>
__host__ __device__ constexpr typename bucket_storage_ref<T, BucketSize, Extent>::extent_type
bucket_storage_ref<T, BucketSize, Extent>::extent() const noexcept
{
  return extent_;
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
constexpr bucket_storage<T, BucketSize, Extent, Allocator>::bucket_storage(
  Extent size, Allocator const& allocator, cuda::stream_ref stream)
  : extent_{size},
    allocator_{allocator},
    slots_{allocator_.allocate(capacity(), stream),
           slot_deleter_type{capacity(), allocator_, stream}}
{
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
constexpr bucket_storage<T, BucketSize, Extent, Allocator>::value_type*
bucket_storage<T, BucketSize, Extent, Allocator>::data() const noexcept
{
  return slots_.get();
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
constexpr bucket_storage<T, BucketSize, Extent, Allocator>::allocator_type
bucket_storage<T, BucketSize, Extent, Allocator>::allocator() const noexcept
{
  return allocator_;
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
constexpr bucket_storage<T, BucketSize, Extent, Allocator>::ref_type
bucket_storage<T, BucketSize, Extent, Allocator>::ref() const noexcept
{
  return ref_type{this->extent(), this->data()};
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
void bucket_storage<T, BucketSize, Extent, Allocator>::initialize(value_type key,
                                                                  cuda::stream_ref stream)
{
  this->initialize_async(key, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
void bucket_storage<T, BucketSize, Extent, Allocator>::initialize_async(value_type key,
                                                                        cuda::stream_ref stream)
{
  if (this->capacity() == 0) { return; }

  auto ftor = cuco::detail::initialize_functor<size_type, T>{this->data(), key};
  CUCO_CUDA_TRY(cub::DeviceFor::Bulk(this->capacity(), ftor, stream.get()));
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
__host__ __device__ constexpr typename bucket_storage<T, BucketSize, Extent, Allocator>::size_type
bucket_storage<T, BucketSize, Extent, Allocator>::num_buckets() const noexcept
{
  return static_cast<size_type>(extent_) / bucket_size;
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
__host__ __device__ constexpr typename bucket_storage<T, BucketSize, Extent, Allocator>::size_type
bucket_storage<T, BucketSize, Extent, Allocator>::capacity() const noexcept
{
  return static_cast<size_type>(extent_);
}

template <typename T, int BucketSize, typename Extent, typename Allocator>
__host__ __device__ constexpr typename bucket_storage<T, BucketSize, Extent, Allocator>::extent_type
bucket_storage<T, BucketSize, Extent, Allocator>::extent() const noexcept
{
  return extent_;
}

}  // namespace cuco
