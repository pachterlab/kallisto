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

#include <cuco/detail/bitwise_compare.cuh>
#include <cuco/detail/utility/cuda.hpp>
#include <cuco/detail/utils.hpp>
#include <cuco/operator.hpp>
#include <cuco/static_multiset_ref.cuh>

#include <cstddef>

namespace cuco {

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  static_multiset(Extent capacity,
                  empty_key<Key> empty_key_sentinel,
                  KeyEqual const& pred,
                  ProbingScheme const& probing_scheme,
                  cuda_thread_scope<Scope>,
                  Storage,
                  Allocator const& alloc,
                  cuda::stream_ref stream)
  : impl_{std::make_unique<impl_type>(
      capacity, empty_key_sentinel, pred, probing_scheme, alloc, stream)}
{
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  static_multiset(Extent n,
                  double desired_load_factor,
                  empty_key<Key> empty_key_sentinel,
                  KeyEqual const& pred,
                  ProbingScheme const& probing_scheme,
                  cuda_thread_scope<Scope>,
                  Storage,
                  Allocator const& alloc,
                  cuda::stream_ref stream)
  : impl_{std::make_unique<impl_type>(
      n, desired_load_factor, empty_key_sentinel, pred, probing_scheme, alloc, stream)}
{
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  static_multiset(Extent capacity,
                  empty_key<Key> empty_key_sentinel,
                  erased_key<Key> erased_key_sentinel,
                  KeyEqual const& pred,
                  ProbingScheme const& probing_scheme,
                  cuda_thread_scope<Scope>,
                  Storage,
                  Allocator const& alloc,
                  cuda::stream_ref stream)
  : impl_{std::make_unique<impl_type>(
      capacity, empty_key_sentinel, erased_key_sentinel, pred, probing_scheme, alloc, stream)}
{
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::clear(
  cuda::stream_ref stream)
{
  impl_->clear(stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::clear_async(
  cuda::stream_ref stream) noexcept
{
  impl_->clear_async(stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::insert(
  InputIt first, InputIt last, cuda::stream_ref stream)
{
  this->insert_async(first, last, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::insert_async(
  InputIt first, InputIt last, cuda::stream_ref stream) noexcept
{
  impl_->insert_async(first, last, ref(op::insert), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate>
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::insert_if(
  InputIt first, InputIt last, StencilIt stencil, Predicate pred, cuda::stream_ref stream)
{
  return impl_->insert_if(first, last, stencil, pred, ref(op::insert), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  insert_if_async(InputIt first,
                  InputIt last,
                  StencilIt stencil,
                  Predicate pred,
                  cuda::stream_ref stream) noexcept
{
  impl_->insert_if_async(first, last, stencil, pred, ref(op::insert), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::contains(
  InputIt first, InputIt last, OutputIt output_begin, cuda::stream_ref stream) const
{
  this->contains_async(first, last, output_begin, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  contains_async(InputIt first,
                 InputIt last,
                 OutputIt output_begin,
                 cuda::stream_ref stream) const noexcept
{
  impl_->contains_async(first, last, output_begin, ref(op::contains), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::contains_if(
  InputIt first,
  InputIt last,
  StencilIt stencil,
  Predicate pred,
  OutputIt output_begin,
  cuda::stream_ref stream) const
{
  this->contains_if_async(first, last, stencil, pred, output_begin, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  contains_if_async(InputIt first,
                    InputIt last,
                    StencilIt stencil,
                    Predicate pred,
                    OutputIt output_begin,
                    cuda::stream_ref stream) const noexcept
{
  impl_->contains_if_async(first, last, stencil, pred, output_begin, ref(op::contains), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::find(
  InputIt first, InputIt last, OutputIt output_begin, cuda::stream_ref stream) const
{
  find_async(first, last, output_begin, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::find_async(
  InputIt first, InputIt last, OutputIt output_begin, cuda::stream_ref stream) const
{
  impl_->find_async(first, last, output_begin, ref(op::find), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::find_if(
  InputIt first,
  InputIt last,
  StencilIt stencil,
  Predicate pred,
  OutputIt output_begin,
  cuda::stream_ref stream) const
{
  this->find_if_async(first, last, stencil, pred, output_begin, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  find_if_async(InputIt first,
                InputIt last,
                StencilIt stencil,
                Predicate pred,
                OutputIt output_begin,
                cuda::stream_ref stream) const
{
  impl_->find_if_async(first, last, stencil, pred, output_begin, ref(op::find), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename CallbackOp>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::for_each(
  CallbackOp&& callback_op, cuda::stream_ref stream) const
{
  impl_->for_each_async(std::forward<CallbackOp>(callback_op), stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename CallbackOp>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  for_each_async(CallbackOp&& callback_op, cuda::stream_ref stream) const
{
  impl_->for_each_async(std::forward<CallbackOp>(callback_op), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename CallbackOp>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::for_each(
  InputIt first, InputIt last, CallbackOp&& callback_op, cuda::stream_ref stream) const
{
  impl_->for_each_async(
    first, last, std::forward<CallbackOp>(callback_op), ref(op::for_each), stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename CallbackOp>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  for_each_async(InputIt first,
                 InputIt last,
                 CallbackOp&& callback_op,
                 cuda::stream_ref stream) const noexcept
{
  impl_->for_each_async(
    first, last, std::forward<CallbackOp>(callback_op), ref(op::for_each), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt>
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::count(
  InputIt first, InputIt last, cuda::stream_ref stream) const
{
  return impl_->count(first, last, ref(op::count), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename ProbeKeyEqual, typename ProbeHash>
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::count(
  InputIt first,
  InputIt last,
  ProbeKeyEqual const& probe_key_equal,
  ProbeHash const& probe_hash,
  cuda::stream_ref stream) const
{
  return impl_->count(
    first,
    last,
    ref(op::count).rebind_key_eq(probe_key_equal).rebind_hash_function(probe_hash),
    stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename ProbeKeyEqual, typename ProbeHash>
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::count_outer(
  InputIt first,
  InputIt last,
  ProbeKeyEqual const& probe_key_equal,
  ProbeHash const& probe_hash,
  cuda::stream_ref stream) const
{
  return impl_->count_outer(
    first,
    last,
    ref(op::count).rebind_key_eq(probe_key_equal).rebind_hash_function(probe_hash),
    stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename ProbeKeyEqual, typename ProbeHash, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::count_each(
  InputIt first,
  InputIt last,
  ProbeKeyEqual const& probe_key_equal,
  ProbeHash const& probe_hash,
  OutputIt output_begin,
  cuda::stream_ref stream) const
{
  impl_->count_each(first,
                    last,
                    output_begin,
                    ref(op::count).rebind_key_eq(probe_key_equal).rebind_hash_function(probe_hash),
                    stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename ProbeKeyEqual, typename ProbeHash, typename OutputIt>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  count_each_outer(InputIt first,
                   InputIt last,
                   ProbeKeyEqual const& probe_key_equal,
                   ProbeHash const& probe_hash,
                   OutputIt output_begin,
                   cuda::stream_ref stream) const
{
  impl_->count_each_outer(
    first,
    last,
    output_begin,
    ref(op::count).rebind_key_eq(probe_key_equal).rebind_hash_function(probe_hash),
    stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <class InputProbeIt, class OutputProbeIt, class OutputMatchIt>
std::pair<OutputProbeIt, OutputMatchIt>
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::retrieve(
  InputProbeIt first,
  InputProbeIt last,
  OutputProbeIt output_probe,
  OutputMatchIt output_match,
  cuda::stream_ref stream) const
{
  return impl_->retrieve(first, last, output_probe, output_match, this->ref(op::retrieve), stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <class InputProbeIt,
          class ProbeEqual,
          class ProbeHash,
          class OutputProbeIt,
          class OutputMatchIt>
std::pair<OutputProbeIt, OutputMatchIt>
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::retrieve(
  InputProbeIt first,
  InputProbeIt last,
  ProbeEqual const& probe_equal,
  ProbeHash const& probe_hash,
  OutputProbeIt output_probe,
  OutputMatchIt output_match,
  cuda::stream_ref stream) const
{
  auto const probe_ref =
    this->ref(op::retrieve).rebind_key_eq(probe_equal).rebind_hash_function(probe_hash);
  return impl_->retrieve(first, last, output_probe, output_match, probe_ref, stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <class InputProbeIt,
          class ProbeEqual,
          class ProbeHash,
          class OutputProbeIt,
          class OutputMatchIt>
std::pair<OutputProbeIt, OutputMatchIt>
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::retrieve_outer(
  InputProbeIt first,
  InputProbeIt last,
  ProbeEqual const& probe_equal,
  ProbeHash const& probe_hash,
  OutputProbeIt output_probe,
  OutputMatchIt output_match,
  cuda::stream_ref stream) const
{
  auto const probe_ref =
    this->ref(op::retrieve).rebind_key_eq(probe_equal).rebind_hash_function(probe_hash);
  return impl_->retrieve_outer(first, last, output_probe, output_match, probe_ref, stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename OutputIt>
OutputIt
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::retrieve_all(
  OutputIt output_begin, cuda::stream_ref stream) const
{
  return impl_->retrieve_all(output_begin, stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::rehash(
  cuda::stream_ref stream)
{
  impl_->rehash(*this, stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::rehash(
  size_type capacity, cuda::stream_ref stream)
{
  auto const extent = make_valid_extent<static_multiset>(capacity);
  impl_->rehash(extent, *this, stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::rehash_async(
  cuda::stream_ref stream)
{
  impl_->rehash_async(*this, stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::rehash_async(
  size_type capacity, cuda::stream_ref stream)
{
  auto const extent = make_valid_extent<static_multiset>(capacity);
  impl_->rehash_async(extent, *this, stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size(
  cuda::stream_ref stream) const
{
  return impl_->size(stream);
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr auto
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::capacity()
  const noexcept
{
  return impl_->capacity();
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
__host__ auto
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::data() const
  -> value_type*
{
  return impl_->data();
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::key_type
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  empty_key_sentinel() const noexcept
{
  return impl_->empty_key_sentinel();
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::key_type
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  erased_key_sentinel() const noexcept
{
  return impl_->erased_key_sentinel();
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  key_equal
  static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::key_eq()
    const noexcept
{
  return impl_->key_eq();
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::hasher
static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::hash_function()
  const noexcept
{
  return impl_->hash_function();
}

template <class Key,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename... Operators>
auto static_multiset<Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::ref(
  Operators...) const noexcept
{
  static_assert(sizeof...(Operators), "No operators specified");
  return cuco::detail::bitwise_compare(this->empty_key_sentinel(), this->erased_key_sentinel())
           ? ref_type<Operators...>{cuco::empty_key<key_type>(this->empty_key_sentinel()),
                                    impl_->key_eq(),
                                    impl_->probing_scheme(),
                                    cuda_thread_scope<Scope>{},
                                    impl_->storage_ref()}
           : ref_type<Operators...>{cuco::empty_key<key_type>(this->empty_key_sentinel()),
                                    cuco::erased_key<key_type>(this->erased_key_sentinel()),
                                    impl_->key_eq(),
                                    impl_->probing_scheme(),
                                    cuda_thread_scope<Scope>{},
                                    impl_->storage_ref()};
}
}  // namespace cuco
