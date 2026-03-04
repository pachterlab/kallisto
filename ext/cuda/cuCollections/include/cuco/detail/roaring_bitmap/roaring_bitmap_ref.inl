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
#include <cuda/std/type_traits>
#include <cuda/stream_ref>

namespace cuco::experimental {

template <class T>
__host__ __device__ roaring_bitmap_ref<T>::roaring_bitmap_ref(storage_ref_type const& storage_ref)
  : impl_{storage_ref}
{
}

template <class T>
template <class U /* = T */,
          class /* = cuda::std::enable_if_t<cuda::std::is_same_v<U, cuda::std::uint32_t>> */>
__device__ roaring_bitmap_ref<T>::roaring_bitmap_ref(cuda::std::byte const* bitmap) : impl_{bitmap}
{
}

template <class T>
template <class InputIt, class OutputIt>
__host__ void roaring_bitmap_ref<T>::contains(InputIt first,
                                              InputIt last,
                                              OutputIt output,
                                              cuda::stream_ref stream) const
{
  impl_.contains(first, last, output, stream);
}

template <class T>
template <class InputIt, class OutputIt>
__host__ void roaring_bitmap_ref<T>::contains_async(InputIt first,
                                                    InputIt last,
                                                    OutputIt output,
                                                    cuda::stream_ref stream) const noexcept
{
  impl_.contains_async(first, last, output, stream);
}

template <class T>
__device__ bool roaring_bitmap_ref<T>::contains(T value) const
{
  return impl_.contains(value);
}

template <class T>
__host__ __device__ cuda::std::size_t roaring_bitmap_ref<T>::size() const noexcept
{
  return impl_.size();
}

template <class T>
__host__ __device__ bool roaring_bitmap_ref<T>::empty() const noexcept
{
  return impl_.empty();
}

template <class T>
__host__ __device__ cuda::std::byte const* roaring_bitmap_ref<T>::data() const noexcept
{
  return impl_.data();
}

template <class T>
__host__ __device__ cuda::std::size_t roaring_bitmap_ref<T>::size_bytes() const noexcept
{
  return impl_.size_bytes();
}

}  // namespace cuco::experimental