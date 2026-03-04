/*
 * Copyright (c) 2021-2025, NVIDIA CORPORATION.
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

#include <cuco/detail/storage/counter_storage.cuh>
#include <cuco/detail/utility/cuda.hpp>
#include <cuco/detail/utils.cuh>

#include <cuda/std/tuple>
#include <thrust/count.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/transform_iterator.h>

#include <iterator>

namespace cuco {
namespace experimental {
template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  static_multimap(Extent capacity,
                  empty_key<Key> empty_key_sentinel,
                  empty_value<T> empty_value_sentinel,
                  KeyEqual const& pred,
                  ProbingScheme const& probing_scheme,
                  cuda_thread_scope<Scope>,
                  Storage,
                  Allocator const& alloc,
                  cuda::stream_ref stream)
  : impl_{std::make_unique<impl_type>(capacity,
                                      cuco::pair{empty_key_sentinel, empty_value_sentinel},
                                      pred,
                                      probing_scheme,
                                      alloc,
                                      stream)},
    empty_value_sentinel_{empty_value_sentinel}
{
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  static_multimap(Extent n,
                  double desired_load_factor,
                  empty_key<Key> empty_key_sentinel,
                  empty_value<T> empty_value_sentinel,
                  KeyEqual const& pred,
                  ProbingScheme const& probing_scheme,
                  cuda_thread_scope<Scope>,
                  Storage,
                  Allocator const& alloc,
                  cuda::stream_ref stream)
  : impl_{std::make_unique<impl_type>(n,
                                      desired_load_factor,
                                      cuco::pair{empty_key_sentinel, empty_value_sentinel},
                                      pred,
                                      probing_scheme,
                                      alloc,
                                      stream)},
    empty_value_sentinel_{empty_value_sentinel}
{
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  static_multimap(Extent capacity,
                  empty_key<Key> empty_key_sentinel,
                  empty_value<T> empty_value_sentinel,
                  erased_key<Key> erased_key_sentinel,
                  KeyEqual const& pred,
                  ProbingScheme const& probing_scheme,
                  cuda_thread_scope<Scope>,
                  Storage,
                  Allocator const& alloc,
                  cuda::stream_ref stream)
  : impl_{std::make_unique<impl_type>(capacity,
                                      cuco::pair{empty_key_sentinel, empty_value_sentinel},
                                      erased_key_sentinel,
                                      pred,
                                      probing_scheme,
                                      alloc,
                                      stream)},
    empty_value_sentinel_{empty_value_sentinel}
{
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::clear(
  cuda::stream_ref stream)
{
  impl_->clear(stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  clear_async(cuda::stream_ref stream) noexcept
{
  impl_->clear_async(stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::insert(
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  insert_async(InputIt first, InputIt last, cuda::stream_ref stream) noexcept
{
  impl_->insert_async(first, last, ref(op::insert), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate>
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::insert_if(
  InputIt first, InputIt last, StencilIt stencil, Predicate pred, cuda::stream_ref stream)
{
  return impl_->insert_if(first, last, stencil, pred, ref(op::insert), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  insert_if_async(InputIt first,
                  InputIt last,
                  StencilIt stencil,
                  Predicate pred,
                  cuda::stream_ref stream) noexcept
{
  impl_->insert_if_async(first, last, stencil, pred, ref(op::insert), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::contains(
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  contains_async(InputIt first,
                 InputIt last,
                 OutputIt output_begin,
                 cuda::stream_ref stream) const noexcept
{
  impl_->contains_async(first, last, output_begin, ref(op::contains), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  contains_if(InputIt first,
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::find(
  InputIt first, InputIt last, OutputIt output_begin, cuda::stream_ref stream) const
{
  this->find_async(first, last, output_begin, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  find_async(InputIt first, InputIt last, OutputIt output_begin, cuda::stream_ref stream) const
{
  impl_->find_async(first, last, output_begin, ref(op::find), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::find_if(
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename CallbackOp>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::for_each(
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename CallbackOp>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  for_each_async(CallbackOp&& callback_op, cuda::stream_ref stream) const
{
  impl_->for_each_async(std::forward<CallbackOp>(callback_op), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename CallbackOp>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::for_each(
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename CallbackOp>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  for_each_async(InputIt first,
                 InputIt last,
                 CallbackOp&& callback_op,
                 cuda::stream_ref stream) const noexcept
{
  impl_->for_each_async(
    first, last, std::forward<CallbackOp>(callback_op), ref(op::for_each), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt>
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::count(
  InputIt first, InputIt last, cuda::stream_ref stream) const
{
  return impl_->count(first, last, ref(op::count), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename ProbeEqual, typename ProbeHash>
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::count(
  InputIt first,
  InputIt last,
  ProbeEqual const& probe_equal,
  ProbeHash const& probe_hash,
  cuda::stream_ref stream) const
{
  return impl_->count(first,
                      last,
                      ref(op::count).rebind_key_eq(probe_equal).rebind_hash_function(probe_hash),
                      stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputIt, typename OutputProbeIt, typename OutputMatchIt>
std::pair<OutputProbeIt, OutputMatchIt>
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::retrieve(
  InputIt first,
  InputIt last,
  OutputProbeIt output_probe,
  OutputMatchIt output_match,
  cuda::stream_ref stream) const
{
  return impl_->retrieve(first, last, output_probe, output_match, this->ref(op::retrieve), stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename InputProbeIt,
          class ProbeEqual,
          class ProbeHash,
          typename OutputProbeIt,
          typename OutputMatchIt>
std::pair<OutputProbeIt, OutputMatchIt>
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::retrieve(
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
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename KeyOut, typename ValueOut>
std::pair<KeyOut, ValueOut>
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::retrieve_all(
  KeyOut keys_out, ValueOut values_out, cuda::stream_ref stream) const
{
  auto const zipped_out_begin = thrust::make_zip_iterator(cuda::std::tuple{keys_out, values_out});
  auto const zipped_out_end   = impl_->retrieve_all(zipped_out_begin, stream);
  auto const num              = std::distance(zipped_out_begin, zipped_out_end);

  return std::make_pair(keys_out + num, values_out + num);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::rehash(
  cuda::stream_ref stream)
{
  impl_->rehash(*this, stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::rehash(
  size_type capacity, cuda::stream_ref stream)
{
  auto const extent = make_valid_extent<static_multimap>(capacity);
  impl_->rehash(extent, *this, stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  rehash_async(cuda::stream_ref stream)
{
  impl_->rehash_async(*this, stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
void static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  rehash_async(size_type capacity, cuda::stream_ref stream)
{
  auto const extent = make_valid_extent<static_multimap>(capacity);
  impl_->rehash_async(extent, *this, stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size_type
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::size(
  cuda::stream_ref stream) const
{
  return impl_->size(stream);
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr auto
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::capacity()
  const noexcept
{
  return impl_->capacity();
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
__host__ auto
static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::data() const
  -> value_type*
{
  return impl_->data();
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  key_type
  static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
    empty_key_sentinel() const noexcept
{
  return impl_->empty_key_sentinel();
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  mapped_type
  static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
    empty_value_sentinel() const noexcept
{
  return empty_value_sentinel_;
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  key_type
  static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
    erased_key_sentinel() const noexcept
{
  return impl_->erased_key_sentinel();
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  key_equal
  static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::key_eq()
    const noexcept
{
  return impl_->key_eq();
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
constexpr static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  hasher
  static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
    hash_function() const noexcept
{
  return impl_->hash_function();
}

template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
template <typename... Operators>
auto static_multimap<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::ref(
  Operators...) const noexcept
{
  static_assert(sizeof...(Operators), "No operators specified");
  return cuco::detail::bitwise_compare(this->empty_key_sentinel(), this->erased_key_sentinel())
           ? ref_type<Operators...>{cuco::empty_key<key_type>(this->empty_key_sentinel()),
                                    cuco::empty_value<mapped_type>(this->empty_value_sentinel()),
                                    impl_->key_eq(),
                                    impl_->probing_scheme(),
                                    cuda_thread_scope<Scope>{},
                                    impl_->storage_ref()}
           : ref_type<Operators...>{cuco::empty_key<key_type>(this->empty_key_sentinel()),
                                    cuco::empty_value<mapped_type>(this->empty_value_sentinel()),
                                    cuco::erased_key<key_type>(this->erased_key_sentinel()),
                                    impl_->key_eq(),
                                    impl_->probing_scheme(),
                                    cuda_thread_scope<Scope>{},
                                    impl_->storage_ref()};
}
}  // namespace experimental

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::static_multimap(
  std::size_t capacity,
  empty_key<Key> empty_key_sentinel,
  empty_value<Value> empty_value_sentinel,
  cudaStream_t stream,
  Allocator const& alloc)
  : capacity_{cuco::detail::get_valid_capacity<cg_size(), vector_width(), uses_vector_load()>(
      capacity)},
    empty_key_sentinel_{empty_key_sentinel.value},
    empty_value_sentinel_{empty_value_sentinel.value},
    allocator_{alloc},
    delete_slots_{allocator_, capacity_, cuda::stream_ref{stream}},
    slots_{allocator_.allocate(capacity_, cuda::stream_ref{stream}), delete_slots_}
{
  auto constexpr block_size = 128;
  auto constexpr stride     = 4;
  auto const grid_size      = (get_capacity() + stride * block_size - 1) / (stride * block_size);

  detail::initialize<atomic_key_type, atomic_mapped_type><<<grid_size, block_size, 0, stream>>>(
    slots_.get(), empty_key_sentinel_, empty_value_sentinel_, get_capacity());
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt>
void static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::insert(InputIt first,
                                                                          InputIt last,
                                                                          cudaStream_t stream)
{
  auto const num_keys = cuco::detail::distance(first, last);
  if (num_keys == 0) { return; }

  auto constexpr block_size = 128;
  auto constexpr stride     = 1;
  auto const grid_size = (cg_size() * num_keys + stride * block_size - 1) / (stride * block_size);
  auto view            = get_device_mutable_view();

  detail::insert<block_size, cg_size()>
    <<<grid_size, block_size, 0, stream>>>(first, num_keys, view);
  CUCO_CUDA_TRY(cudaStreamSynchronize(stream));
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename StencilIt, typename Predicate>
void static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::insert_if(
  InputIt first, InputIt last, StencilIt stencil, Predicate pred, cudaStream_t stream)
{
  auto const num_keys = cuco::detail::distance(first, last);
  if (num_keys == 0) { return; }

  auto constexpr block_size = 128;
  auto constexpr stride     = 1;
  auto const grid_size = (cg_size() * num_keys + stride * block_size - 1) / (stride * block_size);
  auto view            = get_device_mutable_view();

  detail::insert_if_n<block_size, cg_size()>
    <<<grid_size, block_size, 0, stream>>>(first, stencil, num_keys, view, pred);
  CUCO_CUDA_TRY(cudaStreamSynchronize(stream));
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename OutputIt, typename KeyEqual>
void static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::contains(
  InputIt first, InputIt last, OutputIt output_begin, KeyEqual key_equal, cudaStream_t stream) const
{
  auto const num_keys = cuco::detail::distance(first, last);
  if (num_keys == 0) { return; }

  auto constexpr is_pair_contains = false;
  auto constexpr block_size       = 128;
  auto constexpr stride           = 1;
  auto const grid_size = (cg_size() * num_keys + stride * block_size - 1) / (stride * block_size);
  auto view            = get_device_view();

  detail::contains<is_pair_contains, block_size, cg_size()>
    <<<grid_size, block_size, 0, stream>>>(first, num_keys, output_begin, view, key_equal);
  CUCO_CUDA_TRY(cudaStreamSynchronize(stream));
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename OutputIt, typename PairEqual>
void static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::pair_contains(
  InputIt first, InputIt last, OutputIt output_begin, PairEqual pair_equal, cudaStream_t stream)
  const
{
  auto const num_pairs = cuco::detail::distance(first, last);
  if (num_pairs == 0) { return; }

  auto constexpr is_pair_contains = true;
  auto constexpr block_size       = 128;
  auto constexpr stride           = 1;
  auto const grid_size = (cg_size() * num_pairs + stride * block_size - 1) / (stride * block_size);
  auto view            = get_device_view();

  detail::contains<is_pair_contains, block_size, cg_size()>
    <<<grid_size, block_size, 0, stream>>>(first, num_pairs, output_begin, view, pair_equal);
  CUCO_CUDA_TRY(cudaStreamSynchronize(stream));
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename KeyEqual>
std::size_t static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::count(
  InputIt first, InputIt last, cudaStream_t stream, KeyEqual key_equal) const
{
  auto const num_keys = cuco::detail::distance(first, last);
  if (num_keys == 0) { return 0; }

  auto constexpr is_outer   = false;
  auto constexpr block_size = 128;
  auto constexpr stride     = 1;

  auto view            = get_device_view();
  auto const grid_size = (cg_size() * num_keys + stride * block_size - 1) / (stride * block_size);

  auto counter = detail::counter_storage<size_type, Scope, allocator_type>{allocator_, stream};
  counter.reset(stream);

  detail::count<block_size, cg_size(), is_outer>
    <<<grid_size, block_size, 0, stream>>>(first, num_keys, counter.data(), view, key_equal);

  return counter.load_to_host(stream);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename KeyEqual>
std::size_t static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::count_outer(
  InputIt first, InputIt last, cudaStream_t stream, KeyEqual key_equal) const
{
  auto const num_keys = cuco::detail::distance(first, last);
  if (num_keys == 0) { return 0; }

  auto constexpr is_outer   = true;
  auto constexpr block_size = 128;
  auto constexpr stride     = 1;

  auto view            = get_device_view();
  auto const grid_size = (cg_size() * num_keys + stride * block_size - 1) / (stride * block_size);

  auto counter = detail::counter_storage<size_type, Scope, allocator_type>{allocator_, stream};
  counter.reset(stream);

  detail::count<block_size, cg_size(), is_outer>
    <<<grid_size, block_size, 0, stream>>>(first, num_keys, counter.data(), view, key_equal);

  return counter.load_to_host(stream);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename PairEqual>
std::size_t static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::pair_count(
  InputIt first, InputIt last, PairEqual pair_equal, cudaStream_t stream) const
{
  auto const num_pairs = cuco::detail::distance(first, last);
  if (num_pairs == 0) { return 0; }

  auto constexpr is_outer   = false;
  auto constexpr block_size = 128;
  auto constexpr stride     = 1;

  auto view            = get_device_view();
  auto const grid_size = (cg_size() * num_pairs + stride * block_size - 1) / (stride * block_size);

  auto counter = detail::counter_storage<size_type, Scope, allocator_type>{allocator_, stream};
  counter.reset(stream);

  detail::pair_count<block_size, cg_size(), is_outer>
    <<<grid_size, block_size, 0, stream>>>(first, num_pairs, counter.data(), view, pair_equal);

  return counter.load_to_host(stream);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename PairEqual>
std::size_t static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::pair_count_outer(
  InputIt first, InputIt last, PairEqual pair_equal, cudaStream_t stream) const
{
  auto const num_pairs = cuco::detail::distance(first, last);
  if (num_pairs == 0) { return 0; }

  auto constexpr is_outer   = true;
  auto constexpr block_size = 128;
  auto constexpr stride     = 1;

  auto view            = get_device_view();
  auto const grid_size = (cg_size() * num_pairs + stride * block_size - 1) / (stride * block_size);

  auto counter = detail::counter_storage<size_type, Scope, allocator_type>{allocator_, stream};
  counter.reset(stream);

  detail::pair_count<block_size, cg_size(), is_outer>
    <<<grid_size, block_size, 0, stream>>>(first, num_pairs, counter.data(), view, pair_equal);

  return counter.load_to_host(stream);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename OutputIt, typename KeyEqual>
OutputIt static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::retrieve(
  InputIt first, InputIt last, OutputIt output_begin, cudaStream_t stream, KeyEqual key_equal) const
{
  auto const num_keys = cuco::detail::distance(first, last);
  if (num_keys == 0) { return output_begin; }

  // Using per-warp buffer for vector loads and per-CG buffer for scalar loads
  constexpr auto buffer_size = uses_vector_load() ? (warp_size() * 3u) : (cg_size() * 3u);
  constexpr auto is_outer    = false;

  auto view                   = get_device_view();
  auto const flushing_cg_size = [&]() {
    if constexpr (uses_vector_load()) { return warp_size(); }
    return cg_size();
  }();

  auto const grid_size = detail::grid_size(num_keys, cg_size());

  auto counter = detail::counter_storage<size_type, Scope, allocator_type>{allocator_, stream};
  counter.reset(stream);

  detail::retrieve<detail::default_block_size(), flushing_cg_size, cg_size(), buffer_size, is_outer>
    <<<grid_size, detail::default_block_size(), 0, stream>>>(
      first, num_keys, output_begin, counter.data(), view, key_equal);

  return output_begin + counter.load_to_host(stream);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename OutputIt, typename KeyEqual>
OutputIt static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::retrieve_outer(
  InputIt first, InputIt last, OutputIt output_begin, cudaStream_t stream, KeyEqual key_equal) const
{
  auto const num_keys = cuco::detail::distance(first, last);
  if (num_keys == 0) { return output_begin; }

  // Using per-warp buffer for vector loads and per-CG buffer for scalar loads
  constexpr auto buffer_size = uses_vector_load() ? (warp_size() * 3u) : (cg_size() * 3u);
  constexpr auto is_outer    = true;

  auto view                   = get_device_view();
  auto const flushing_cg_size = [&]() {
    if constexpr (uses_vector_load()) { return warp_size(); }
    return cg_size();
  }();

  auto const grid_size = detail::grid_size(num_keys, cg_size());

  auto counter = detail::counter_storage<size_type, Scope, allocator_type>{allocator_, stream};
  counter.reset(stream);

  detail::retrieve<detail::default_block_size(), flushing_cg_size, cg_size(), buffer_size, is_outer>
    <<<grid_size, detail::default_block_size(), 0, stream>>>(
      first, num_keys, output_begin, counter.data(), view, key_equal);

  return output_begin + counter.load_to_host(stream);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename OutputIt1, typename OutputIt2, typename PairEqual>
std::pair<OutputIt1, OutputIt2>
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::pair_retrieve(
  InputIt first,
  InputIt last,
  OutputIt1 probe_output_begin,
  OutputIt2 contained_output_begin,
  PairEqual pair_equal,
  cudaStream_t stream) const
{
  auto const num_pairs = cuco::detail::distance(first, last);
  if (num_pairs == 0) { return std::make_pair(probe_output_begin, contained_output_begin); }

  // Using per-warp buffer for vector loads and per-CG buffer for scalar loads
  constexpr auto buffer_size = uses_vector_load() ? (warp_size() * 3u) : (cg_size() * 3u);
  constexpr auto block_size  = 128;
  constexpr auto is_outer    = false;
  constexpr auto stride      = 1;

  auto view                   = get_device_view();
  auto const flushing_cg_size = [&]() {
    if constexpr (uses_vector_load()) { return warp_size(); }
    return cg_size();
  }();
  auto const grid_size = (cg_size() * num_pairs + stride * block_size - 1) / (stride * block_size);

  auto counter = detail::counter_storage<size_type, Scope, allocator_type>{allocator_, stream};
  counter.reset(stream);

  detail::pair_retrieve<block_size, flushing_cg_size, cg_size(), buffer_size, is_outer>
    <<<grid_size, block_size, 0, stream>>>(first,
                                           num_pairs,
                                           probe_output_begin,
                                           contained_output_begin,
                                           counter.data(),
                                           view,
                                           pair_equal);

  auto const h_count = counter.load_to_host(stream);
  return {probe_output_begin + h_count, contained_output_begin + h_count};
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename InputIt, typename OutputIt1, typename OutputIt2, typename PairEqual>
std::pair<OutputIt1, OutputIt2>
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::pair_retrieve_outer(
  InputIt first,
  InputIt last,
  OutputIt1 probe_output_begin,
  OutputIt2 contained_output_begin,
  PairEqual pair_equal,
  cudaStream_t stream) const
{
  auto const num_pairs = cuco::detail::distance(first, last);
  if (num_pairs == 0) { return std::make_pair(probe_output_begin, contained_output_begin); }

  // Using per-warp buffer for vector loads and per-CG buffer for scalar loads
  constexpr auto buffer_size = uses_vector_load() ? (warp_size() * 3u) : (cg_size() * 3u);
  constexpr auto block_size  = 128;
  constexpr auto is_outer    = true;
  constexpr auto stride      = 1;

  auto view                   = get_device_view();
  auto const flushing_cg_size = [&]() {
    if constexpr (uses_vector_load()) { return warp_size(); }
    return cg_size();
  }();
  auto const grid_size = (cg_size() * num_pairs + stride * block_size - 1) / (stride * block_size);

  auto counter = detail::counter_storage<size_type, Scope, allocator_type>{allocator_, stream};
  counter.reset(stream);

  detail::pair_retrieve<block_size, flushing_cg_size, cg_size(), buffer_size, is_outer>
    <<<grid_size, block_size, 0, stream>>>(first,
                                           num_pairs,
                                           probe_output_begin,
                                           contained_output_begin,
                                           counter.data(),
                                           view,
                                           pair_equal);

  auto const h_count = counter.load_to_host(stream);
  return {probe_output_begin + h_count, contained_output_begin + h_count};
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename ParentCG>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_mutable_view::insert(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
  value_type const& insert_pair) noexcept
{
  impl_.template insert<uses_vector_load()>(g, insert_pair);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename CG>
__device__ __forceinline__ static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::make_copy(
  CG g, pair_atomic_type* const memory_to_use, device_view source_device_view) noexcept
{
#if defined(CUCO_HAS_CUDA_BARRIER)
  __shared__ cuda::barrier<cuda::thread_scope::thread_scope_block> barrier;
  if (g.thread_rank() == 0) { init(&barrier, g.size()); }
  g.sync();

  cuda::memcpy_async(g,
                     memory_to_use,
                     source_device_view.get_slots(),
                     sizeof(pair_atomic_type) * source_device_view.get_capacity(),
                     barrier);

  barrier.arrive_and_wait();
#else
  pair_atomic_type const* const slots_ptr = source_device_view.get_slots();
  for (std::size_t i = g.thread_rank(); i < source_device_view.get_capacity(); i += g.size()) {
    new (&memory_to_use[i].first)
      atomic_key_type{slots_ptr[i].first.load(cuda::memory_order_relaxed)};
    new (&memory_to_use[i].second)
      atomic_mapped_type{slots_ptr[i].second.load(cuda::memory_order_relaxed)};
  }
  g.sync();
#endif

  return device_view(memory_to_use,
                     source_device_view.get_capacity(),
                     source_device_view.get_empty_key_sentinel(),
                     source_device_view.get_empty_value_sentinel());
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename CG, typename atomicT, typename OutputIt>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::flush_output_buffer(
  CG g,
  uint32_t const num_outputs,
  value_type* output_buffer,
  atomicT* num_matches,
  OutputIt output_begin) noexcept
{
  impl_.flush_output_buffer(g, num_outputs, output_buffer, num_matches, output_begin);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename CG, typename atomicT, typename OutputIt1, typename OutputIt2>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::flush_output_buffer(
  CG g,
  uint32_t const num_outputs,
  value_type* probe_output_buffer,
  value_type* contained_output_buffer,
  atomicT* num_matches,
  OutputIt1 probe_output_begin,
  OutputIt2 contained_output_begin) noexcept
{
  impl_.flush_output_buffer(g,
                            num_outputs,
                            probe_output_buffer,
                            contained_output_buffer,
                            num_matches,
                            probe_output_begin,
                            contained_output_begin);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename ProbeKey, typename KeyEqual, typename ParentCG>
__device__ __forceinline__ bool
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::contains(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
  ProbeKey const& k,
  KeyEqual key_equal) const noexcept
{
  constexpr bool is_pair_contains = false;
  return impl_.template contains<is_pair_contains, uses_vector_load()>(g, k, key_equal);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename ProbePair, typename PairEqual, typename ParentCG>
__device__ __forceinline__ bool
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::pair_contains(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
  ProbePair const& p,
  PairEqual pair_equal) const noexcept
{
  constexpr bool is_pair_contains = true;
  return impl_.template contains<is_pair_contains, uses_vector_load()>(g, p, pair_equal);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename KeyEqual, typename ParentCG>
__device__ __forceinline__ std::size_t
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::count(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
  Key const& k,
  KeyEqual key_equal) noexcept
{
  constexpr bool is_outer = false;
  return impl_.count<uses_vector_load(), is_outer>(g, k, key_equal);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename KeyEqual, typename ParentCG>
__device__ __forceinline__ std::size_t
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::count_outer(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
  Key const& k,
  KeyEqual key_equal) noexcept
{
  constexpr bool is_outer = true;
  return impl_.count<uses_vector_load(), is_outer>(g, k, key_equal);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename PairEqual, typename ParentCG>
__device__ __forceinline__ std::size_t
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::pair_count(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
  value_type const& pair,
  PairEqual pair_equal) noexcept
{
  constexpr bool is_outer = false;
  return impl_.pair_count<uses_vector_load(), is_outer>(g, pair, pair_equal);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename PairEqual, typename ParentCG>
__device__ __forceinline__ std::size_t
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::pair_count_outer(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
  value_type const& pair,
  PairEqual pair_equal) noexcept
{
  constexpr bool is_outer = true;
  return impl_.pair_count<uses_vector_load(), is_outer>(g, pair, pair_equal);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <uint32_t buffer_size,
          typename FlushingCG,
          typename atomicT,
          typename OutputIt,
          typename KeyEqual,
          typename ParentCG>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::retrieve(
  FlushingCG flushing_cg,
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
  Key const& k,
  uint32_t* flushing_cg_counter,
  value_type* output_buffer,
  atomicT* num_matches,
  OutputIt output_begin,
  KeyEqual key_equal) noexcept
{
  constexpr bool is_outer = false;
  if constexpr (uses_vector_load()) {
    impl_.template retrieve<buffer_size, is_outer>(flushing_cg,
                                                   probing_cg,
                                                   k,
                                                   flushing_cg_counter,
                                                   output_buffer,
                                                   num_matches,
                                                   output_begin,
                                                   key_equal);
  } else  // In the case of scalar load, flushing CG is the same as probing CG
  {
    impl_.template retrieve<buffer_size, is_outer>(
      probing_cg, k, flushing_cg_counter, output_buffer, num_matches, output_begin, key_equal);
  }
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <uint32_t buffer_size,
          typename FlushingCG,
          typename atomicT,
          typename OutputIt,
          typename KeyEqual,
          typename ParentCG>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::retrieve_outer(
  FlushingCG flushing_cg,
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
  Key const& k,
  uint32_t* flushing_cg_counter,
  value_type* output_buffer,
  atomicT* num_matches,
  OutputIt output_begin,
  KeyEqual key_equal) noexcept
{
  constexpr bool is_outer = true;
  if constexpr (uses_vector_load()) {
    impl_.template retrieve<buffer_size, is_outer>(flushing_cg,
                                                   probing_cg,
                                                   k,
                                                   flushing_cg_counter,
                                                   output_buffer,
                                                   num_matches,
                                                   output_begin,
                                                   key_equal);
  } else  // In the case of scalar load, flushing CG is the same as probing CG
  {
    impl_.template retrieve<buffer_size, is_outer>(
      probing_cg, k, flushing_cg_counter, output_buffer, num_matches, output_begin, key_equal);
  }
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename OutputIt1,
          typename OutputIt2,
          typename OutputIt3,
          typename OutputIt4,
          typename PairEqual,
          typename ParentCG>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::pair_retrieve(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
  value_type const& pair,
  OutputIt1 probe_key_begin,
  OutputIt2 probe_val_begin,
  OutputIt3 contained_key_begin,
  OutputIt4 contained_val_begin,
  PairEqual pair_equal) noexcept
{
  constexpr bool is_outer = false;
  impl_.pair_retrieve<is_outer, uses_vector_load()>(probing_cg,
                                                    pair,
                                                    probe_key_begin,
                                                    probe_val_begin,
                                                    contained_key_begin,
                                                    contained_val_begin,
                                                    pair_equal);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <uint32_t buffer_size,
          typename FlushingCG,
          typename atomicT,
          typename OutputIt1,
          typename OutputIt2,
          typename PairEqual,
          typename ParentCG>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::pair_retrieve(
  FlushingCG flushing_cg,
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
  value_type const& pair,
  uint32_t* flushing_cg_counter,
  value_type* probe_output_buffer,
  value_type* contained_output_buffer,
  atomicT* num_matches,
  OutputIt1 probe_output_begin,
  OutputIt2 contained_output_begin,
  PairEqual pair_equal) noexcept
{
  constexpr bool is_outer = false;
  if constexpr (uses_vector_load()) {
    impl_.pair_retrieve<buffer_size, is_outer>(flushing_cg,
                                               probing_cg,
                                               pair,
                                               flushing_cg_counter,
                                               probe_output_buffer,
                                               contained_output_buffer,
                                               num_matches,
                                               probe_output_begin,
                                               contained_output_begin,
                                               pair_equal);
  } else  // In the case of scalar load, flushing CG is the same as probing CG
  {
    impl_.pair_retrieve<buffer_size, is_outer>(probing_cg,
                                               pair,
                                               flushing_cg_counter,
                                               probe_output_buffer,
                                               contained_output_buffer,
                                               num_matches,
                                               probe_output_begin,
                                               contained_output_begin,
                                               pair_equal);
  }
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <typename OutputIt1,
          typename OutputIt2,
          typename OutputIt3,
          typename OutputIt4,
          typename PairEqual,
          typename ParentCG>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::pair_retrieve_outer(
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
  value_type const& pair,
  OutputIt1 probe_key_begin,
  OutputIt2 probe_val_begin,
  OutputIt3 contained_key_begin,
  OutputIt4 contained_val_begin,
  PairEqual pair_equal) noexcept
{
  constexpr bool is_outer = true;
  impl_.pair_retrieve<is_outer, uses_vector_load()>(probing_cg,
                                                    pair,
                                                    probe_key_begin,
                                                    probe_val_begin,
                                                    contained_key_begin,
                                                    contained_val_begin,
                                                    pair_equal);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
template <uint32_t buffer_size,
          typename FlushingCG,
          typename atomicT,
          typename OutputIt1,
          typename OutputIt2,
          typename PairEqual,
          typename ParentCG>
__device__ __forceinline__ void
static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view::pair_retrieve_outer(
  FlushingCG flushing_cg,
  cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
  value_type const& pair,
  uint32_t* flushing_cg_counter,
  value_type* probe_output_buffer,
  value_type* contained_output_buffer,
  atomicT* num_matches,
  OutputIt1 probe_output_begin,
  OutputIt2 contained_output_begin,
  PairEqual pair_equal) noexcept
{
  constexpr bool is_outer = true;
  if constexpr (uses_vector_load()) {
    impl_.pair_retrieve<buffer_size, is_outer>(flushing_cg,
                                               probing_cg,
                                               pair,
                                               flushing_cg_counter,
                                               probe_output_buffer,
                                               contained_output_buffer,
                                               num_matches,
                                               probe_output_begin,
                                               contained_output_begin,
                                               pair_equal);
  } else  // In the case of scalar load, flushing CG is the same as probing CG
  {
    impl_.pair_retrieve<buffer_size, is_outer>(probing_cg,
                                               pair,
                                               flushing_cg_counter,
                                               probe_output_buffer,
                                               contained_output_buffer,
                                               num_matches,
                                               probe_output_begin,
                                               contained_output_begin,
                                               pair_equal);
  }
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
std::size_t static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::get_size(
  cudaStream_t stream) const noexcept
{
  auto begin  = thrust::make_transform_iterator(raw_slots(), detail::slot_to_tuple<Key, Value>{});
  auto filled = cuco::detail::slot_is_filled<Key>{get_empty_key_sentinel()};

  return thrust::count_if(thrust::cuda::par.on(stream), begin, begin + get_capacity(), filled);
}

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
float static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::get_load_factor(
  cudaStream_t stream) const noexcept
{
  auto size = get_size(stream);
  return static_cast<float>(size) / static_cast<float>(capacity_);
}

}  // namespace cuco
