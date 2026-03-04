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

template <class T, cuda::thread_scope Scope, class Hash>
__host__ __device__ constexpr hyperloglog_ref<T, Scope, Hash>::hyperloglog_ref(
  cuda::std::span<cuda::std::byte> sketch_span, Hash const& hash)
  : impl_{sketch_span, hash}
{
}

template <class T, cuda::thread_scope Scope, class Hash>
template <class CG>
__device__ constexpr void hyperloglog_ref<T, Scope, Hash>::clear(CG group) noexcept
{
  impl_.clear(group);
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ constexpr void hyperloglog_ref<T, Scope, Hash>::clear_async(
  cuda::stream_ref stream) noexcept
{
  impl_.clear_async(stream);
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ constexpr void hyperloglog_ref<T, Scope, Hash>::clear(cuda::stream_ref stream)
{
  impl_.clear(stream);
}

template <class T, cuda::thread_scope Scope, class Hash>
__device__ constexpr void hyperloglog_ref<T, Scope, Hash>::add(T const& item) noexcept
{
  impl_.add(item);
}

template <class T, cuda::thread_scope Scope, class Hash>
template <class InputIt>
__host__ constexpr void hyperloglog_ref<T, Scope, Hash>::add_async(InputIt first,
                                                                   InputIt last,
                                                                   cuda::stream_ref stream)
{
  impl_.add_async(first, last, stream);
}

template <class T, cuda::thread_scope Scope, class Hash>
template <class InputIt>
__host__ constexpr void hyperloglog_ref<T, Scope, Hash>::add(InputIt first,
                                                             InputIt last,
                                                             cuda::stream_ref stream)
{
  impl_.add(first, last, stream);
}

template <class T, cuda::thread_scope Scope, class Hash>
template <class CG, cuda::thread_scope OtherScope>
__device__ constexpr void hyperloglog_ref<T, Scope, Hash>::merge(
  CG group, hyperloglog_ref<T, OtherScope, Hash> const& other)
{
  impl_.merge(group, other.impl_);
}

template <class T, cuda::thread_scope Scope, class Hash>
template <cuda::thread_scope OtherScope>
__host__ constexpr void hyperloglog_ref<T, Scope, Hash>::merge_async(
  hyperloglog_ref<T, OtherScope, Hash> const& other, cuda::stream_ref stream)
{
  impl_.merge_async(other.impl_, stream);
}

template <class T, cuda::thread_scope Scope, class Hash>
template <cuda::thread_scope OtherScope>
__host__ constexpr void hyperloglog_ref<T, Scope, Hash>::merge(
  hyperloglog_ref<T, OtherScope, Hash> const& other, cuda::stream_ref stream)
{
  impl_.merge(other.impl_, stream);
}

template <class T, cuda::thread_scope Scope, class Hash>
__device__ std::size_t hyperloglog_ref<T, Scope, Hash>::estimate(
  cooperative_groups::thread_block const& group) const noexcept
{
  return impl_.estimate(group);
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ constexpr std::size_t hyperloglog_ref<T, Scope, Hash>::estimate(
  cuda::stream_ref stream) const
{
  return impl_.estimate(stream);
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ __device__ constexpr auto hyperloglog_ref<T, Scope, Hash>::hash_function() const noexcept
{
  return impl_.hash_function();
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ __device__ constexpr cuda::std::span<cuda::std::byte>
hyperloglog_ref<T, Scope, Hash>::sketch() const noexcept
{
  return impl_.sketch();
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ __device__ constexpr std::size_t hyperloglog_ref<T, Scope, Hash>::sketch_bytes()
  const noexcept
{
  return impl_.sketch_bytes();
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ __device__ constexpr std::size_t hyperloglog_ref<T, Scope, Hash>::sketch_bytes(
  cuco::sketch_size_kb sketch_size_kb) noexcept
{
  return impl_type::sketch_bytes(sketch_size_kb);
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ __device__ constexpr std::size_t hyperloglog_ref<T, Scope, Hash>::sketch_bytes(
  cuco::standard_deviation standard_deviation) noexcept
{
  return impl_type::sketch_bytes(standard_deviation);
}

template <class T, cuda::thread_scope Scope, class Hash>
__host__ __device__ constexpr std::size_t
hyperloglog_ref<T, Scope, Hash>::sketch_alignment() noexcept
{
  return impl_type::sketch_alignment();
}

}  // namespace cuco