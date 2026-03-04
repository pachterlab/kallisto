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

#include <cuco/operator.hpp>

#include <cuda/atomic>
#include <cuda/std/type_traits>
#include <cuda/std/utility>

#include <cooperative_groups.h>

namespace cuco {

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_multiset_ref<
  Key,
  Scope,
  KeyEqual,
  ProbingScheme,
  StorageRef,
  Operators...>::static_multiset_ref(cuco::empty_key<Key> empty_key_sentinel,
                                     KeyEqual const& predicate,
                                     ProbingScheme const& probing_scheme,
                                     cuda_thread_scope<Scope>,
                                     StorageRef storage_ref) noexcept
  : impl_{empty_key_sentinel, predicate, probing_scheme, storage_ref}
{
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_multiset_ref<
  Key,
  Scope,
  KeyEqual,
  ProbingScheme,
  StorageRef,
  Operators...>::static_multiset_ref(cuco::empty_key<Key> empty_key_sentinel,
                                     cuco::erased_key<Key> erased_key_sentinel,
                                     KeyEqual const& predicate,
                                     ProbingScheme const& probing_scheme,
                                     cuda_thread_scope<Scope>,
                                     StorageRef storage_ref) noexcept
  : impl_{empty_key_sentinel, erased_key_sentinel, predicate, probing_scheme, storage_ref}
{
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename... OtherOperators>
__host__ __device__ constexpr static_multiset_ref<Key,
                                                  Scope,
                                                  KeyEqual,
                                                  ProbingScheme,
                                                  StorageRef,
                                                  Operators...>::
  static_multiset_ref(
    static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, OtherOperators...>&&
      other) noexcept
  : impl_{std::move(other.impl_)}
{
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_multiset_ref<Key,
                                                  Scope,
                                                  KeyEqual,
                                                  ProbingScheme,
                                                  StorageRef,
                                                  Operators...>::key_equal
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::key_eq()
  const noexcept
{
  return this->impl_.key_eq();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_multiset_ref<Key,
                                                  Scope,
                                                  KeyEqual,
                                                  ProbingScheme,
                                                  StorageRef,
                                                  Operators...>::hasher
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::hash_function()
  const noexcept
{
  return impl_.hash_function();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr auto
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::capacity()
  const noexcept
{
  return impl_.capacity();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr auto
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::storage_ref()
  const noexcept
{
  return this->impl_.storage_ref();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr auto
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::probing_scheme()
  const noexcept
{
  return this->impl_.probing_scheme();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_multiset_ref<Key,
                                                  Scope,
                                                  KeyEqual,
                                                  ProbingScheme,
                                                  StorageRef,
                                                  Operators...>::extent_type
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::extent()
  const noexcept
{
  return impl_.extent();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_multiset_ref<Key,
                                                  Scope,
                                                  KeyEqual,
                                                  ProbingScheme,
                                                  StorageRef,
                                                  Operators...>::extent_type
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::bucket_extent()
  const noexcept
{
  return this->extent();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr Key
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::
  empty_key_sentinel() const noexcept
{
  return impl_.empty_key_sentinel();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr Key
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::
  erased_key_sentinel() const noexcept
{
  return impl_.erased_key_sentinel();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_multiset_ref<Key,
                                                  Scope,
                                                  KeyEqual,
                                                  ProbingScheme,
                                                  StorageRef,
                                                  Operators...>::const_iterator
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::end()
  const noexcept
{
  return this->impl_.end();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_multiset_ref<Key,
                                                  Scope,
                                                  KeyEqual,
                                                  ProbingScheme,
                                                  StorageRef,
                                                  Operators...>::iterator
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::end() noexcept
{
  return this->impl_.end();
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename... NewOperators>
__host__ __device__ constexpr auto
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::
  rebind_operators(NewOperators...) const noexcept
{
  return static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, NewOperators...>{
    cuco::empty_key<Key>{this->empty_key_sentinel()},
    this->key_eq(),
    this->probing_scheme(),
    {},
    this->storage_ref()};
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename NewKeyEqual>
__host__ __device__ constexpr auto
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::rebind_key_eq(
  NewKeyEqual const& key_equal) const noexcept
{
  return static_multiset_ref<Key, Scope, NewKeyEqual, ProbingScheme, StorageRef, Operators...>{
    cuco::empty_key<Key>{this->empty_key_sentinel()},
    key_equal,
    this->probing_scheme(),
    {},
    this->storage_ref()};
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename NewHash>
__host__ __device__ constexpr auto
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::
  rebind_hash_function(NewHash const& hash) const
{
  auto const probing_scheme = this->probing_scheme().rebind_hash_function(hash);
  return static_multiset_ref<Key,
                             Scope,
                             KeyEqual,
                             cuda::std::decay_t<decltype(probing_scheme)>,
                             StorageRef,
                             Operators...>{cuco::empty_key<Key>{this->empty_key_sentinel()},
                                           this->key_eq(),
                                           probing_scheme,
                                           {},
                                           this->storage_ref()};
}
template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename CG, cuda::thread_scope NewScope>
__device__ constexpr auto
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::make_copy(
  CG tile, bucket_type* const memory_to_use, cuda_thread_scope<NewScope> scope) const noexcept
{
  auto const storage_ref = this->storage_ref().make_copy(tile, memory_to_use);
  return static_multiset_ref<Key,
                             NewScope,
                             KeyEqual,
                             ProbingScheme,
                             decltype(storage_ref),
                             Operators...>{cuco::empty_key<Key>{this->empty_key_sentinel()},
                                           this->key_eq(),
                                           this->probing_scheme(),
                                           scope,
                                           storage_ref};
}

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename CG>
__device__ constexpr void
static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::initialize(
  CG tile) noexcept
{
  this->storage_ref().initialize(tile, this->empty_key_sentinel());
}

namespace detail {

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::insert_tag,
  static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type =
    static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type   = typename base_type::key_type;
  using value_type = typename base_type::value_type;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Inserts an element.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param value The element to insert
   *
   * @return True if the given element is successfully inserted
   */
  template <typename Value>
  __device__ bool insert(Value value) noexcept
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);
    return ref_.impl_.insert(value);
  }

  /**
   * @brief Inserts an element.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group insert
   * @param value The element to insert
   *
   * @return True if the given element is successfully inserted
   */
  template <typename Value, typename ParentCG>
  __device__ bool insert(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                         Value value) noexcept
  {
    auto& ref_ = static_cast<ref_type&>(*this);
    if (!cuco::detail::bitwise_compare(ref_.erased_key_sentinel(), ref_.empty_key_sentinel())) {
      return ref_.impl_.template insert<true>(group, value);
    } else {
      return ref_.impl_.template insert<false>(group, value);
    }
  }
};

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::contains_tag,
  static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type =
    static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type   = typename base_type::key_type;
  using value_type = typename base_type::value_type;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Indicates whether the probe key `key` was inserted into the container.
   *
   * @tparam ProbeKey Probe key type
   *
   * @param key The key to search for
   *
   * @return A boolean indicating whether the probe key is present
   */
  template <typename ProbeKey>
  [[nodiscard]] __device__ bool contains(ProbeKey key) const noexcept
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.contains(key);
  }

  /**
   * @brief Indicates whether the probe key `key` was inserted into the container.
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group contains
   * @param key The key to search for
   *
   * @return A boolean indicating whether the probe key is present
   */
  template <typename ProbeKey, typename ParentCG>
  [[nodiscard]] __device__ bool contains(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> group, ProbeKey key) const noexcept
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.contains(group, key);
  }
};

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::find_tag,
  static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type =
    static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type       = typename base_type::key_type;
  using value_type     = typename base_type::value_type;
  using iterator       = typename base_type::iterator;
  using const_iterator = typename base_type::const_iterator;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Finds an element in the multiset with key equivalent to the probe key.
   *
   * @note Returns a un-incrementable input iterator to the element whose key is equivalent to
   * `key`. If no such element exists, returns `end()`.
   *
   * @tparam ProbeKey Probe key type
   *
   * @param key The key to search for
   *
   * @return An iterator to the position at which the equivalent key is stored
   */
  template <typename ProbeKey>
  [[nodiscard]] __device__ const_iterator find(ProbeKey key) const noexcept
  {
    // CRTP: cast `this` to the actual ref type
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.find(key);
  }

  /**
   * @brief Finds an element in the multiset with key equivalent to the probe key.
   *
   * @note Returns a un-incrementable input iterator to the element whose key is equivalent to
   * `key`. If no such element exists, returns `end()`.
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform this operation
   * @param key The key to search for
   *
   * @return An iterator to the position at which the equivalent key is stored
   */
  template <typename ProbeKey, typename ParentCG>
  [[nodiscard]] __device__ const_iterator
  find(cooperative_groups::thread_block_tile<cg_size, ParentCG> group, ProbeKey key) const noexcept
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.find(group, key);
  }
};

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::retrieve_tag,
  static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type =
    static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type       = typename base_type::key_type;
  using value_type     = typename base_type::value_type;
  using iterator       = typename base_type::iterator;
  using const_iterator = typename base_type::const_iterator;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Retrieves all the slots corresponding to all keys in the range `[input_probe_begin,
   * input_probe_end)`.
   *
   * If key `k = *(first + i)` exists in the container, copies `k` to `output_probe` and associated
   * slot contents to `output_match`, respectively. The output order is unspecified.
   *
   * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
   * Use `count()` to determine the size of the output range.
   *
   * @tparam BlockSize Size of the thread block this operation is executed in
   * @tparam InputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam AtomicCounter Atomic counter type that follows the same semantics as
   * `cuda::atomic(_ref)`
   *
   * @param block Thread block this operation is executed in
   * @param input_probe_begin Beginning of the input sequence of keys
   * @param input_probe_end End of the input sequence of keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param atomic_counter Counter that is used to determine the next free position in the output
   * sequences
   */
  template <int32_t BlockSize,
            class InputProbeIt,
            class OutputProbeIt,
            class OutputMatchIt,
            class AtomicCounter>
  __device__ void retrieve(cooperative_groups::thread_block const& block,
                           InputProbeIt input_probe_begin,
                           InputProbeIt input_probe_end,
                           OutputProbeIt output_probe,
                           OutputMatchIt output_match,
                           AtomicCounter& atomic_counter) const
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    ref_.impl_.template retrieve<BlockSize>(
      block, input_probe_begin, input_probe_end, output_probe, output_match, atomic_counter);
  }

  /**
   * @brief Retrieves all the slots corresponding to all keys in the range `[input_probe_begin,
   * input_probe_end)`.
   *
   * If key `k = *(first + i)` exists in the container, copies `k` to `output_probe` and associated
   * slot contents to `output_match`, respectively. The output order is unspecified.
   *
   * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
   * Use `count_outer()` to determine the size of the output range.
   *
   * If a key `k` has no matches in the container, then `{key, empty_slot_sentinel}` will be added
   * to the output sequence.
   *
   * @tparam BlockSize Size of the thread block this operation is executed in
   * @tparam InputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam AtomicCounter Atomic counter type that follows the same semantics as
   * `cuda::atomic(_ref)`
   *
   * @param block Thread block this operation is executed in
   * @param input_probe_begin Beginning of the input sequence of keys
   * @param input_probe_end End of the input sequence of keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param atomic_counter Counter that is used to determine the next free position in the output
   * sequences
   */
  template <int32_t BlockSize,
            class InputProbeIt,
            class OutputProbeIt,
            class OutputMatchIt,
            class AtomicCounter>
  __device__ void retrieve_outer(cooperative_groups::thread_block const& block,
                                 InputProbeIt input_probe_begin,
                                 InputProbeIt input_probe_end,
                                 OutputProbeIt output_probe,
                                 OutputMatchIt output_match,
                                 AtomicCounter& atomic_counter) const
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    ref_.impl_.template retrieve_outer<BlockSize>(
      block, input_probe_begin, input_probe_end, output_probe, output_match, atomic_counter);
  }

  /**
   * @brief Retrieves all the slots corresponding to all keys in the range `[input_probe_begin,
   * input_probe_end)` if `pred` of the corresponding stencil returns true.
   *
   * If key `k = *(first + i)` exists in the container and `pred( *(stencil + i) )` returns true,
   * copies `k` to `output_probe` and associated slot content to `output_match`, respectively.
   * The output order is unspecified.
   *
   * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
   * Use `count()` to determine the size of the output range.
   *
   * @tparam BlockSize Size of the thread block this operation is executed in
   * @tparam InputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool`
   * and argument type is convertible from `std::iterator_traits<StencilIt>::value_type`
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam AtomicCounter Atomic counter type that follows the same semantics as
   * `cuda::atomic(_ref)`
   *
   * @param block Thread block this operation is executed in
   * @param input_probe_begin Beginning of the input sequence of keys
   * @param input_probe_end End of the input sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil + n)`
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param atomic_counter Counter that is used to determine the next free position in the output
   * sequences
   */
  template <int BlockSize,
            class InputProbeIt,
            class StencilIt,
            class Predicate,
            class OutputProbeIt,
            class OutputMatchIt,
            class AtomicCounter>
  __device__ void retrieve_if(cooperative_groups::thread_block const& block,
                              InputProbeIt input_probe_begin,
                              InputProbeIt input_probe_end,
                              StencilIt stencil,
                              Predicate pred,
                              OutputProbeIt output_probe,
                              OutputMatchIt output_match,
                              AtomicCounter& atomic_counter) const
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    ref_.impl_.template retrieve_if<BlockSize>(block,
                                               input_probe_begin,
                                               input_probe_end,
                                               stencil,
                                               pred,
                                               output_probe,
                                               output_match,
                                               atomic_counter);
  }
};

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::for_each_tag,
  static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type =
    static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;

  static constexpr auto cg_size = base_type::cg_size;

 public:
  /**
   * @brief Executes a callback on every element in the container with key equivalent to the probe
   * key.
   *
   * @note Passes an un-incrementable input iterator to the element whose key is equivalent to
   * `key` to the callback.
   *
   * @tparam ProbeKey Probe key type
   * @tparam CallbackOp Unary callback functor or device lambda
   *
   * @param key The key to search for
   * @param callback_op Function to call on every element found
   */
  template <class ProbeKey, class CallbackOp>
  __device__ void for_each(ProbeKey key, CallbackOp&& callback_op) const noexcept
  {
    // CRTP: cast `this` to the actual ref type
    auto const& ref_ = static_cast<ref_type const&>(*this);
    ref_.impl_.for_each(key, cuda::std::forward<CallbackOp>(callback_op));
  }

  /**
   * @brief Executes a callback on every element in the container with key equivalent to the probe
   * key.
   *
   * @note Passes an un-incrementable input iterator to the element whose key is equivalent to
   * `key` to the callback.
   *
   * @note This function uses cooperative group semantics, meaning that any thread may call the
   * callback if it finds a matching element. If multiple elements are found within the same group,
   * each thread with a match will call the callback with its associated element.
   *
   * @note Synchronizing `group` within `callback_op` is undefined behavior.
   *
   * @tparam ProbeKey Probe key type
   * @tparam CallbackOp Unary callback functor or device lambda
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform this operation
   * @param key The key to search for
   * @param callback_op Function to call on every element found
   */
  template <class ProbeKey, class CallbackOp, typename ParentCG>
  __device__ void for_each(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                           ProbeKey key,
                           CallbackOp&& callback_op) const noexcept
  {
    // CRTP: cast `this` to the actual ref type
    auto const& ref_ = static_cast<ref_type const&>(*this);
    ref_.impl_.for_each(group, key, cuda::std::forward<CallbackOp>(callback_op));
  }

  /**
   * @brief Executes a callback on every element in the container with key equivalent to the probe
   * key and can additionally perform work that requires synchronizing the Cooperative Group
   * performing this operation.
   *
   * @note Passes an un-incrementable input iterator to the element whose key is equivalent to
   * `key` to the callback.
   *
   * @note This function uses cooperative group semantics, meaning that any thread may call the
   * callback if it finds a matching element. If multiple elements are found within the same group,
   * each thread with a match will call the callback with its associated element.
   *
   * @note Synchronizing `group` within `callback_op` is undefined behavior.
   *
   * @note The `sync_op` function can be used to perform work that requires synchronizing threads in
   * `group` inbetween probing steps, where the number of probing steps performed between
   * synchronization points is capped by `bucket_size * cg_size`. The functor will be called right
   * after the current probing bucket has been traversed.
   *
   * @tparam ProbeKey Probe key type
   * @tparam CallbackOp Unary callback functor or device lambda
   * @tparam SyncOp Functor or device lambda which accepts the current `group` object
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform this operation
   * @param key The key to search for
   * @param callback_op Function to call on every element found
   * @param sync_op Function that is allowed to synchronize `group` inbetween probing buckets
   */
  template <class ProbeKey, class CallbackOp, class SyncOp, typename ParentCG>
  __device__ void for_each(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                           ProbeKey key,
                           CallbackOp&& callback_op,
                           SyncOp&& sync_op) const noexcept
  {
    // CRTP: cast `this` to the actual ref type
    auto const& ref_ = static_cast<ref_type const&>(*this);
    ref_.impl_.for_each(
      group, key, cuda::std::forward<CallbackOp>(callback_op), cuda::std::forward<SyncOp>(sync_op));
  }
};

template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::count_tag,
  static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type =
    static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type   = typename base_type::key_type;
  using value_type = typename base_type::value_type;
  using size_type  = typename base_type::size_type;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Counts the occurrence of a given key contained in multiset
   *
   * @tparam ProbeKey Probe key type
   *
   * @param key The key to count for
   *
   * @return Number of occurrences found by the current thread
   */
  template <typename ProbeKey>
  __device__ size_type count(ProbeKey key) const noexcept
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.count(key);
  }

  /**
   * @brief Counts the occurrence of a given key contained in multiset
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group count
   * @param key The key to count for
   *
   * @return Number of occurrences found by the current thread
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ size_type count(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                             ProbeKey key) const noexcept
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.count(group, key);
  }
};

}  // namespace detail
}  // namespace cuco
