/*
 * Copyright (c) 2024-2025, NVIDIA CORPORATION.
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

#include <cuco/utility/cuda_thread_scope.cuh>

#include <cuda/atomic>
#include <cuda/stream_ref>

namespace cuco {

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
__host__ __device__ constexpr bloom_filter_ref<Key, Extent, Scope, Policy>::bloom_filter_ref(
  filter_block_type* data, Extent num_blocks, cuda_thread_scope<Scope>, Policy const& policy)
  : impl_{data, num_blocks, {}, policy}
{
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
__host__ __device__ constexpr bloom_filter_ref<Key, Extent, Scope, Policy>::bloom_filter_ref(
  word_type* data, Extent num_blocks, cuda_thread_scope<Scope>, Policy const& policy)
  : impl_{data, num_blocks, {}, policy}
{
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class CG>
__device__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::clear(CG group)
{
  impl_.clear(group);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::clear(cuda::stream_ref stream)
{
  impl_.clear(stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::clear_async(
  cuda::stream_ref stream)
{
  impl_.clear_async(stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class ProbeKey>
__device__ void bloom_filter_ref<Key, Extent, Scope, Policy>::add(ProbeKey const& key)
{
  impl_.add(key);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class CG, class ProbeKey>
__device__ void bloom_filter_ref<Key, Extent, Scope, Policy>::add(CG group, ProbeKey const& key)
{
  impl_.add(group, key);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class CG, class InputIt>
__device__ void bloom_filter_ref<Key, Extent, Scope, Policy>::add(CG group,
                                                                  InputIt first,
                                                                  InputIt last)
{
  impl_.add(group, first, last);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class InputIt>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::add(InputIt first,
                                                                          InputIt last,
                                                                          cuda::stream_ref stream)
{
  impl_.add(first, last, stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class InputIt>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::add_async(
  InputIt first, InputIt last, cuda::stream_ref stream)
{
  impl_.add_async(first, last, stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class InputIt, class StencilIt, class Predicate>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::add_if(
  InputIt first, InputIt last, StencilIt stencil, Predicate pred, cuda::stream_ref stream)
{
  impl_.add_if(first, last, stencil, pred, stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class InputIt, class StencilIt, class Predicate>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::add_if_async(
  InputIt first, InputIt last, StencilIt stencil, Predicate pred, cuda::stream_ref stream) noexcept
{
  impl_.add_if_async(first, last, stencil, pred, stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class ProbeKey>
[[nodiscard]] __device__ bool bloom_filter_ref<Key, Extent, Scope, Policy>::contains(
  ProbeKey const& key) const
{
  return impl_.contains(key);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class CG, class ProbeKey>
[[nodiscard]] __device__ bool bloom_filter_ref<Key, Extent, Scope, Policy>::contains(
  CG group, ProbeKey const& key) const
{
  return impl_.contains(group, key);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class InputIt, class OutputIt>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::contains(
  InputIt first, InputIt last, OutputIt output_begin, cuda::stream_ref stream) const
{
  impl_.contains(first, last, output_begin, stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class InputIt, class OutputIt>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::contains_async(
  InputIt first, InputIt last, OutputIt output_begin, cuda::stream_ref stream) const noexcept
{
  impl_.contains_async(first, last, output_begin, stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class InputIt, class StencilIt, class Predicate, class OutputIt>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::contains_if(
  InputIt first,
  InputIt last,
  StencilIt stencil,
  Predicate pred,
  OutputIt output_begin,
  cuda::stream_ref stream) const
{
  impl_.contains_if(first, last, stencil, pred, output_begin, stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
template <class InputIt, class StencilIt, class Predicate, class OutputIt>
__host__ constexpr void bloom_filter_ref<Key, Extent, Scope, Policy>::contains_if_async(
  InputIt first,
  InputIt last,
  StencilIt stencil,
  Predicate pred,
  OutputIt output_begin,
  cuda::stream_ref stream) const noexcept
{
  impl_.contains_if_async(first, last, stencil, pred, output_begin, stream);
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
[[nodiscard]] __host__ __device__ constexpr
  typename bloom_filter_ref<Key, Extent, Scope, Policy>::word_type*
  bloom_filter_ref<Key, Extent, Scope, Policy>::data() noexcept
{
  return impl_.data();
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
[[nodiscard]] __host__ __device__ constexpr
  typename bloom_filter_ref<Key, Extent, Scope, Policy>::word_type const*
  bloom_filter_ref<Key, Extent, Scope, Policy>::data() const noexcept
{
  return impl_.data();
}

template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
[[nodiscard]] __host__ __device__ constexpr
  typename bloom_filter_ref<Key, Extent, Scope, Policy>::extent_type
  bloom_filter_ref<Key, Extent, Scope, Policy>::block_extent() const noexcept
{
  return impl_.block_extent();
}

}  // namespace cuco