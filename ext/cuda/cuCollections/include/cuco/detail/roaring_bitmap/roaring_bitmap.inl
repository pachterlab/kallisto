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

#include <cuda/std/cstddef>
#include <cuda/stream_ref>

namespace cuco::experimental {

template <class T, class Allocator>
roaring_bitmap<T, Allocator>::roaring_bitmap(cuda::std::byte const* bitmap,
                                             Allocator const& alloc,
                                             cuda::stream_ref stream)
  : storage_{bitmap, alloc, stream}
{
}

template <class T, class Allocator>
template <class InputIt, class OutputIt>
void roaring_bitmap<T, Allocator>::contains(InputIt first,
                                            InputIt last,
                                            OutputIt output,
                                            cuda::stream_ref stream) const
{
  ref_type{storage_.ref()}.contains(first, last, output, stream);
}

template <class T, class Allocator>
template <class InputIt, class OutputIt>
void roaring_bitmap<T, Allocator>::contains_async(InputIt first,
                                                  InputIt last,
                                                  OutputIt output,
                                                  cuda::stream_ref stream) const noexcept
{
  ref_type{storage_.ref()}.contains_async(first, last, output, stream);
}

template <class T, class Allocator>
cuda::std::size_t roaring_bitmap<T, Allocator>::size() const noexcept
{
  return ref_type{storage_.ref()}.size();
}

template <class T, class Allocator>
bool roaring_bitmap<T, Allocator>::empty() const noexcept
{
  return ref_type{storage_.ref()}.empty();
}

template <class T, class Allocator>
cuda::std::byte const* roaring_bitmap<T, Allocator>::data() const noexcept
{
  return ref_type{storage_.ref()}.data();
}

template <class T, class Allocator>
cuda::std::size_t roaring_bitmap<T, Allocator>::size_bytes() const noexcept
{
  return ref_type{storage_.ref()}.size_bytes();
}

template <class T, class Allocator>
typename roaring_bitmap<T, Allocator>::allocator_type roaring_bitmap<T, Allocator>::allocator()
  const noexcept
{
  return storage_.allocator();
}

template <class T, class Allocator>
typename roaring_bitmap<T, Allocator>::ref_type roaring_bitmap<T, Allocator>::ref() const noexcept
{
  return ref_type{storage_.ref()};
}
}  // namespace cuco::experimental