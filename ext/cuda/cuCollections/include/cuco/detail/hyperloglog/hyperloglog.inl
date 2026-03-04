/*
 * Copyright (c) 2024, NVIDIA CORPORATION.
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

namespace cuco {

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr hyperloglog<T, Scope, Hash, Allocator>::hyperloglog(cuco::sketch_size_kb sketch_size_kb,
                                                              Hash const& hash,
                                                              Allocator const& alloc,
                                                              cuda::stream_ref stream)
  : allocator_{alloc},
    sketch_{allocator_.allocate(sketch_bytes(sketch_size_kb) / sizeof(register_type), stream),
            detail::custom_deleter{
              sketch_bytes(sketch_size_kb) / sizeof(register_type), allocator_, stream}},
    ref_{cuda::std::span{reinterpret_cast<cuda::std::byte*>(sketch_.get()),
                         sketch_bytes(sketch_size_kb)},
         hash}
{
  this->clear_async(stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr hyperloglog<T, Scope, Hash, Allocator>::hyperloglog(
  cuco::standard_deviation standard_deviation,
  Hash const& hash,
  Allocator const& alloc,
  cuda::stream_ref stream)
  : allocator_{alloc},
    sketch_{allocator_.allocate(sketch_bytes(standard_deviation) / sizeof(register_type), stream),
            detail::custom_deleter{
              sketch_bytes(standard_deviation) / sizeof(register_type), allocator_, stream}},
    ref_{cuda::std::span{reinterpret_cast<cuda::std::byte*>(sketch_.get()),
                         sketch_bytes(standard_deviation)},
         hash}
{
  this->clear_async(stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr void hyperloglog<T, Scope, Hash, Allocator>::clear_async(cuda::stream_ref stream) noexcept
{
  ref_.clear_async(stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr void hyperloglog<T, Scope, Hash, Allocator>::clear(cuda::stream_ref stream)
{
  ref_.clear(stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
template <class InputIt>
constexpr void hyperloglog<T, Scope, Hash, Allocator>::add_async(InputIt first,
                                                                 InputIt last,
                                                                 cuda::stream_ref stream)
{
  ref_.add_async(first, last, stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
template <class InputIt>
constexpr void hyperloglog<T, Scope, Hash, Allocator>::add(InputIt first,
                                                           InputIt last,
                                                           cuda::stream_ref stream)
{
  ref_.add(first, last, stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
template <cuda::thread_scope OtherScope, class OtherAllocator>
constexpr void hyperloglog<T, Scope, Hash, Allocator>::merge_async(
  hyperloglog<T, OtherScope, Hash, OtherAllocator> const& other, cuda::stream_ref stream)
{
  ref_.merge_async(other.ref_, stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
template <cuda::thread_scope OtherScope, class OtherAllocator>
constexpr void hyperloglog<T, Scope, Hash, Allocator>::merge(
  hyperloglog<T, OtherScope, Hash, OtherAllocator> const& other, cuda::stream_ref stream)
{
  ref_.merge(other.ref_, stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
template <cuda::thread_scope OtherScope>
constexpr void hyperloglog<T, Scope, Hash, Allocator>::merge_async(
  ref_type<OtherScope> const& other_ref, cuda::stream_ref stream)
{
  ref_.merge_async(other_ref, stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
template <cuda::thread_scope OtherScope>
constexpr void hyperloglog<T, Scope, Hash, Allocator>::merge(ref_type<OtherScope> const& other_ref,
                                                             cuda::stream_ref stream)
{
  ref_.merge(other_ref, stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr std::size_t hyperloglog<T, Scope, Hash, Allocator>::estimate(
  cuda::stream_ref stream) const
{
  return ref_.estimate(stream);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr typename hyperloglog<T, Scope, Hash, Allocator>::template ref_type<>
hyperloglog<T, Scope, Hash, Allocator>::ref() const noexcept
{
  return {this->sketch(), this->hash_function()};
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr auto hyperloglog<T, Scope, Hash, Allocator>::hash_function() const noexcept
{
  return ref_.hash_function();
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr cuda::std::span<cuda::std::byte> hyperloglog<T, Scope, Hash, Allocator>::sketch()
  const noexcept
{
  return ref_.sketch();
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr size_t hyperloglog<T, Scope, Hash, Allocator>::sketch_bytes() const noexcept
{
  return ref_.sketch_bytes();
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr size_t hyperloglog<T, Scope, Hash, Allocator>::sketch_bytes(
  cuco::sketch_size_kb sketch_size_kb) noexcept
{
  return ref_type<>::sketch_bytes(sketch_size_kb);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr size_t hyperloglog<T, Scope, Hash, Allocator>::sketch_bytes(
  cuco::standard_deviation standard_deviation) noexcept
{
  return ref_type<>::sketch_bytes(standard_deviation);
}

template <class T, cuda::thread_scope Scope, class Hash, class Allocator>
constexpr size_t hyperloglog<T, Scope, Hash, Allocator>::sketch_alignment() noexcept
{
  return ref_type<>::sketch_alignment();
}

}  // namespace cuco