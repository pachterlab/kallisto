/*
 * Copyright (c) 2023-2025, NVIDIA CORPORATION.
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

#include <cuco/detail/bitwise_compare.cuh>
#include <cuco/detail/equal_wrapper.cuh>
#include <cuco/operator.hpp>

#include <cuda/atomic>
#include <cuda/std/iterator>
#include <cuda/std/type_traits>
#include <cuda/std/utility>

#include <cooperative_groups.h>

namespace cuco {

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_map_ref<
  Key,
  T,
  Scope,
  KeyEqual,
  ProbingScheme,
  StorageRef,
  Operators...>::static_map_ref(cuco::empty_key<Key> empty_key_sentinel,
                                cuco::empty_value<T> empty_value_sentinel,
                                KeyEqual const& predicate,
                                ProbingScheme const& probing_scheme,
                                cuda_thread_scope<Scope>,
                                StorageRef storage_ref) noexcept
  : impl_{
      cuco::pair{empty_key_sentinel, empty_value_sentinel}, predicate, probing_scheme, storage_ref}
{
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_map_ref<
  Key,
  T,
  Scope,
  KeyEqual,
  ProbingScheme,
  StorageRef,
  Operators...>::static_map_ref(cuco::empty_key<Key> empty_key_sentinel,
                                cuco::empty_value<T> empty_value_sentinel,
                                cuco::erased_key<Key> erased_key_sentinel,
                                KeyEqual const& predicate,
                                ProbingScheme const& probing_scheme,
                                cuda_thread_scope<Scope>,
                                StorageRef storage_ref) noexcept
  : impl_{cuco::pair{empty_key_sentinel, empty_value_sentinel},
          erased_key_sentinel,
          predicate,
          probing_scheme,
          storage_ref}
{
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename... OtherOperators>
__host__ __device__ constexpr static_map_ref<Key,
                                             T,
                                             Scope,
                                             KeyEqual,
                                             ProbingScheme,
                                             StorageRef,
                                             Operators...>::
  static_map_ref(
    static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, OtherOperators...>&&
      other) noexcept
  : impl_{std::move(other.impl_)}
{
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_map_ref<Key,
                                             T,
                                             Scope,
                                             KeyEqual,
                                             ProbingScheme,
                                             StorageRef,
                                             Operators...>::key_equal
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::key_eq()
  const noexcept
{
  return this->impl_.key_eq();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_map_ref<Key,
                                             T,
                                             Scope,
                                             KeyEqual,
                                             ProbingScheme,
                                             StorageRef,
                                             Operators...>::hasher
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::hash_function()
  const noexcept
{
  return impl_.hash_function();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_map_ref<Key,
                                             T,
                                             Scope,
                                             KeyEqual,
                                             ProbingScheme,
                                             StorageRef,
                                             Operators...>::const_iterator
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::end()
  const noexcept
{
  return this->impl_.end();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_map_ref<Key,
                                             T,
                                             Scope,
                                             KeyEqual,
                                             ProbingScheme,
                                             StorageRef,
                                             Operators...>::iterator
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::end() noexcept
{
  return this->impl_.end();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr auto
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::capacity()
  const noexcept
{
  return impl_.capacity();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr auto
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::storage_ref()
  const noexcept
{
  return this->impl_.storage_ref();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr auto
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::probing_scheme()
  const noexcept
{
  return this->impl_.probing_scheme();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_map_ref<Key,
                                             T,
                                             Scope,
                                             KeyEqual,
                                             ProbingScheme,
                                             StorageRef,
                                             Operators...>::extent_type
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::extent()
  const noexcept
{
  return impl_.extent();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr static_map_ref<Key,
                                             T,
                                             Scope,
                                             KeyEqual,
                                             ProbingScheme,
                                             StorageRef,
                                             Operators...>::extent_type
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::bucket_extent()
  const noexcept
{
  return this->extent();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr Key
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::
  empty_key_sentinel() const noexcept
{
  return impl_.empty_key_sentinel();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr T
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::
  empty_value_sentinel() const noexcept
{
  return impl_.empty_value_sentinel();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
__host__ __device__ constexpr Key
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::
  erased_key_sentinel() const noexcept
{
  return impl_.erased_key_sentinel();
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename... NewOperators>
__host__ __device__ constexpr auto
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::rebind_operators(
  NewOperators...) const noexcept
{
  return static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, NewOperators...>{
    cuco::empty_key<Key>{this->empty_key_sentinel()},
    cuco::empty_value<T>{this->empty_value_sentinel()},
    this->key_eq(),
    this->probing_scheme(),
    {},
    this->storage_ref()};
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename NewKeyEqual>
__host__ __device__ constexpr auto
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::rebind_key_eq(
  NewKeyEqual const& key_equal) const noexcept
{
  return static_map_ref<Key, T, Scope, NewKeyEqual, ProbingScheme, StorageRef, Operators...>{
    cuco::empty_key<Key>{this->empty_key_sentinel()},
    cuco::empty_value<T>{this->empty_value_sentinel()},
    key_equal,
    this->probing_scheme(),
    {},
    this->storage_ref()};
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename NewHash>
__host__ __device__ constexpr auto
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::
  rebind_hash_function(NewHash const& hash) const
{
  auto const probing_scheme = this->probing_scheme().rebind_hash_function(hash);
  return static_map_ref<Key,
                        T,
                        Scope,
                        KeyEqual,
                        cuda::std::decay_t<decltype(probing_scheme)>,
                        StorageRef,
                        Operators...>{cuco::empty_key<Key>{this->empty_key_sentinel()},
                                      cuco::empty_value<T>{this->empty_value_sentinel()},
                                      this->key_eq(),
                                      probing_scheme,
                                      {},
                                      this->storage_ref()};
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename CG, cuda::thread_scope NewScope>
__device__ constexpr auto
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::make_copy(
  CG tile,
  typename StorageRef::value_type* const memory_to_use,
  cuda_thread_scope<NewScope> scope) const noexcept
{
  this->impl_.make_copy(tile, memory_to_use);
  return static_map_ref<Key, T, NewScope, KeyEqual, ProbingScheme, StorageRef, Operators...>{
    cuco::empty_key<Key>{this->empty_key_sentinel()},
    cuco::empty_value<T>{this->empty_value_sentinel()},
    cuco::erased_key<Key>{this->erased_key_sentinel()},
    this->key_eq(),
    this->probing_scheme(),
    scope,
    storage_ref_type{this->extent(), memory_to_use}};
}

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
template <typename CG>
__device__ constexpr void
static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>::initialize(
  CG tile) noexcept
{
  this->impl_.initialize(tile);
}

namespace detail {

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::insert_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type  = typename base_type::value_type;
  using mapped_type = T;

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
    if (ref_.erased_key_sentinel() != ref_.empty_key_sentinel()) {
      return ref_.impl_.template insert<true>(group, value);
    } else {
      return ref_.impl_.template insert<false>(group, value);
    }
  }
};

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::insert_or_assign_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type  = typename base_type::value_type;
  using mapped_type = T;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, assigns `v`
   * to the mapped_type corresponding to the key `k`.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param value The element to insert
   */
  template <typename Value>
  __device__ void insert_or_assign(Value value) noexcept
  {
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");

    ref_type& ref_ = static_cast<ref_type&>(*this);

    auto const val            = ref_.impl_.heterogeneous_value(value);
    auto const key            = ref_.impl_.extract_key(val);
    auto const probing_scheme = ref_.impl_.probing_scheme();
    auto storage_ref          = ref_.impl_.storage_ref();
    auto probing_iter =
      probing_scheme.template make_iterator<bucket_size>(key, storage_ref.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref[*probing_iter];

      for (auto& slot_content : bucket_slots) {
        auto const eq_res =
          ref_.impl_.predicate_.template operator()<is_insert::YES>(key, slot_content.first);
        auto const intra_bucket_index = cuda::std::distance(bucket_slots.begin(), &slot_content);
        auto slot_ptr                 = ref_.impl_.get_slot_ptr(*probing_iter, intra_bucket_index);

        // If the key is already in the container, update the payload and return
        if (eq_res == detail::equal_result::EQUAL) {
          cuda::atomic_ref<mapped_type, Scope> payload_ref(slot_ptr->second);
          payload_ref.store(val.second, cuda::memory_order_relaxed);
          return;
        }
        if (eq_res == detail::equal_result::AVAILABLE) {
          if (attempt_insert_or_assign(slot_ptr, val)) { return; }
        }
      }
      ++probing_iter;
      if (*probing_iter == init_idx) { return; }
    }
  }

  /**
   * @brief Inserts an element.
   *
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, assigns `v`
   * to the mapped_type corresponding to the key `k`.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group insert
   * @param value The element to insert
   */
  template <typename Value, typename ParentCG>
  __device__ void insert_or_assign(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                                   Value value) noexcept
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);

    auto const val            = ref_.impl_.heterogeneous_value(value);
    auto const key            = ref_.impl_.extract_key(val);
    auto const probing_scheme = ref_.impl_.probing_scheme();
    auto storage_ref          = ref_.impl_.storage_ref();
    auto probing_iter =
      probing_scheme.template make_iterator<bucket_size>(group, key, storage_ref.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref[*probing_iter];

      auto const [state, intra_bucket_index] = [&]() {
        auto res = detail::equal_result::UNEQUAL;
        for (auto i = 0; i < bucket_size; ++i) {
          res =
            ref_.impl_.predicate_.template operator()<is_insert::YES>(key, bucket_slots[i].first);
          if (res != detail::equal_result::UNEQUAL) {
            return detail::bucket_probing_results{res, i};
          }
        }
        // returns dummy index `-1` for UNEQUAL
        return detail::bucket_probing_results{res, -1};
      }();

      auto slot_ptr = ref_.impl_.get_slot_ptr(*probing_iter, intra_bucket_index);

      auto const group_contains_equal = group.ballot(state == detail::equal_result::EQUAL);
      if (group_contains_equal) {
        auto const src_lane = __ffs(group_contains_equal) - 1;
        if (group.thread_rank() == src_lane) {
          cuda::atomic_ref<mapped_type, Scope> payload_ref(slot_ptr->second);
          payload_ref.store(val.second, cuda::memory_order_relaxed);
        }
        group.sync();
        return;
      }

      auto const group_contains_available = group.ballot(state == detail::equal_result::AVAILABLE);
      if (group_contains_available) {
        auto const src_lane = __ffs(group_contains_available) - 1;
        auto const status =
          (group.thread_rank() == src_lane) ? attempt_insert_or_assign(slot_ptr, val) : false;

        // Exit if inserted or assigned
        if (group.shfl(status, src_lane)) { return; }
      } else {
        ++probing_iter;
        if (*probing_iter == init_idx) { return; }
      }
    }
  }

 private:
  /**
   * @brief Attempts to insert an element into a slot or update the matching payload with the given
   * element
   *
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, assigns `v`
   * to the mapped_type corresponding to the key `k`.
   *
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   *
   * @param group The Cooperative Group used to perform group insert
   * @param value The element to insert
   *
   * @return Returns `true` if the given `value` is inserted or `value` has a match in the map.
   */
  template <typename Value>
  __device__ constexpr bool attempt_insert_or_assign(value_type* slot, Value value) noexcept
  {
    ref_type& ref_    = static_cast<ref_type&>(*this);
    auto expected_key = ref_.impl_.empty_slot_sentinel().first;

    cuda::atomic_ref<key_type, Scope> key_ref(slot->first);
    auto const success = key_ref.compare_exchange_strong(
      expected_key, static_cast<key_type>(value.first), cuda::memory_order_relaxed);

    // if key success or key was already present in the map
    if (success or (ref_.impl_.predicate().equal_to(value.first, expected_key) ==
                    detail::equal_result::EQUAL)) {
      // Update payload
      cuda::atomic_ref<mapped_type, Scope> payload_ref(slot->second);
      payload_ref.store(value.second, cuda::memory_order_relaxed);
      return true;
    }
    return false;
  }
};

// TODO use insert_or_apply internally
template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::insert_or_apply_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type = typename base_type::value_type;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, applies
   * `Op` binary function to the mapped_type corresponding to the key `k` and the value `v`.
   *
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.

   * @param value The element to insert
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   *
   * @return Returns `true` if the given `value` is inserted successfully.
   */

  template <typename Value, typename Op>
  __device__ bool insert_or_apply(Value value, Op op)
  {
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");

    static_assert(
      cuda::std::is_invocable_v<Op, cuda::atomic_ref<T, Scope>, T>,
      "insert_or_apply expects `Op` to be a callable as `Op(cuda::atomic_ref<T, Scope>, T)`");

    auto& ref_ = static_cast<ref_type&>(*this);

    // directly dispatch to implementation if no init element is given
    auto constexpr use_direct_apply = false;
    return ref_.insert_or_apply_impl<use_direct_apply>(value, op);
  }

  /**
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, applies
   * `Op` binary function to the mapped_type corresponding to the key `k` and the value `v`.
   *
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.

   * @param value The element to insert
   * @param init The init value of the op
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   *
   * @return Returns `true` if the given `value` is inserted successfully.
   */
  template <typename Value,
            typename Init,
            typename Op,
            typename = cuda::std::enable_if_t<std::is_convertible_v<Value, value_type>>>
  __device__ bool insert_or_apply(Value value, Init init, Op op)
  {
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");

    static_assert(
      cuda::std::is_invocable_v<Op, cuda::atomic_ref<T, Scope>, T>,
      "insert_or_apply expects `Op` to be a callable as `Op(cuda::atomic_ref<T, Scope>, T)`");

    auto& ref_ = static_cast<ref_type&>(*this);
    return ref_.dispatch_insert_or_apply(value, init, op);
  }

  /**
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, applies
   * `Op` binary function to the mapped_type corresponding to the key `k` and the value `v`.
   *
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group insert
   * @param value The element to insert
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   *
   * @return Returns `true` if the given `value` is inserted successfully.
   */

  template <typename Value, typename Op, typename ParentCG>
  __device__ bool insert_or_apply(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                                  Value value,
                                  Op op)
  {
    static_assert(
      cuda::std::is_invocable_v<Op, cuda::atomic_ref<T, Scope>, T>,
      "insert_or_apply expects `Op` to be a callable as `Op(cuda::atomic_ref<T, Scope>, T)`");
    auto& ref_ = static_cast<ref_type&>(*this);

    auto constexpr use_direct_apply = false;
    return ref_.insert_or_apply_impl<use_direct_apply>(group, value, op);
  }

  /**
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, applies
   * `Op` binary function to the mapped_type corresponding to the key `k` and the value `v`.
   *
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group insert
   * @param value The element to insert
   * @param init The init value of the op
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   *
   * @return Returns `true` if the given `value` is inserted successfully.
   */
  template <typename Value, typename Init, typename Op, typename ParentCG>
  __device__ bool insert_or_apply(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                                  Value value,
                                  Init init,
                                  Op op)
  {
    auto& ref_ = static_cast<ref_type&>(*this);
    static_assert(
      cuda::std::is_invocable_v<Op, cuda::atomic_ref<T, Scope>, T>,
      "insert_or_apply expects `Op` to be a callable as `Op(cuda::atomic_ref<T, Scope>, T)`");

    return ref_.dispatch_insert_or_apply(group, value, init, op);
  }

 private:
  /**
   * @brief dispatches `insert_or_apply_impl` based on condition `init == empty_value_sentinel`.
   *
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.
   *
   * @param value The element to insert
   * @param init The init value of the op
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   *
   * @return Returns `true` if the given `value` is inserted successfully.
   */
  template <typename Value, typename Init, typename Op>
  __device__ bool dispatch_insert_or_apply(Value value, Init init, Op op)
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);
    // if init equals sentinel value, then we can just `apply` op instead of write
    if (cuco::detail::bitwise_compare(init, ref_.empty_value_sentinel())) {
      return ref_.template insert_or_apply_impl<true>(value, op);
    } else {
      return ref_.template insert_or_apply_impl<false>(value, op);
    }
  }

  /**
   * @brief dispatches `insert_or_apply_impl` based on condition `init == empty_value_sentinel`.
   *
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group insert
   * @param value The element to insert
   * @param init The init value of the op
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   */
  template <typename Value, typename Init, typename Op, typename ParentCG>
  __device__ bool dispatch_insert_or_apply(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> group, Value value, Init init, Op op)
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);
    // if init equals sentinel value, then we can just `apply` op instead of write
    if (cuco::detail::bitwise_compare(init, ref_.empty_value_sentinel())) {
      return ref_.insert_or_apply_impl<true>(group, value, op);
    } else {
      return ref_.insert_or_apply_impl<false>(group, value, op);
    }
  }

  /**
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, applies
   * `Op` binary function to the mapped_type corresponding to the key `k` and the value `v`.
   *
   * @tparam UseDirectApply Boolean condition which enables direct apply `op` instead of
   * wait_for_payload with an atomic_store on an empty slot.
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.
   *
   * @param value The element to insert
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   *
   * @return Returns `true` if the given `value` is inserted successfully.
   */
  template <bool UseDirectApply, typename Value, typename Op>
  __device__ bool insert_or_apply_impl(Value value, Op op)
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);

    auto const val            = ref_.impl_.heterogeneous_value(value);
    auto const key            = ref_.impl_.extract_key(val);
    auto const probing_scheme = ref_.impl_.probing_scheme();
    auto storage_ref          = ref_.impl_.storage_ref();
    auto probing_iter =
      probing_scheme.template make_iterator<bucket_size>(key, storage_ref.extent());
    auto const init_idx    = *probing_iter;
    auto const empty_value = ref_.empty_value_sentinel();

    // wait for payload only when init != sentinel and insert strategy is not `packed_cas`
    auto constexpr wait_for_payload = (not UseDirectApply) and (sizeof(value_type) > 8);

    while (true) {
      auto const bucket_slots = storage_ref[*probing_iter];

      for (auto& slot_content : bucket_slots) {
        auto const eq_res =
          ref_.impl_.predicate_.template operator()<is_insert::YES>(key, slot_content.first);
        auto const intra_bucket_index = cuda::std::distance(bucket_slots.begin(), &slot_content);
        auto slot_ptr                 = ref_.impl_.get_slot_ptr(*probing_iter, intra_bucket_index);

        // If the key is already in the container, update the payload and return
        if (eq_res == detail::equal_result::EQUAL) {
          // wait for payload only when performing insert operation
          if constexpr (wait_for_payload) {
            ref_.impl_.wait_for_payload(slot_ptr->second, empty_value);
          }
          op(cuda::atomic_ref<T, Scope>{slot_ptr->second}, val.second);
          return false;
        }
        if (eq_res == detail::equal_result::AVAILABLE) {
          switch (ref_.template attempt_insert_or_apply<UseDirectApply>(
            slot_ptr, slot_content, val, op)) {
            case insert_result::SUCCESS: return true;
            case insert_result::DUPLICATE: {
              // wait for payload only when performing insert operation
              if constexpr (wait_for_payload) {
                ref_.impl_.wait_for_payload(slot_ptr->second, empty_value);
              }
              op(cuda::atomic_ref<T, Scope>{slot_ptr->second}, val.second);
              return false;
            }
            default: continue;
          }
        }
      }
      ++probing_iter;
      if (*probing_iter == init_idx) { return false; }
    }
  }

  /**
   * @brief Inserts a key-value pair `{k, v}` if it's not present in the map. Otherwise, applies
   * `Op` binary function to the mapped_type corresponding to the key `k` and the value `v`.
   *
   * @tparam UseDirectApply Boolean condition which enables direct apply `op` instead of
   * wait_for_payload with an atomic_store on an empty slot.
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group insert
   * @param value The element to insert
   * @param init The init value of the op
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   *
   * @return Returns `true` if the given `value` is inserted successfully.
   */
  template <bool UseDirectApply, typename Value, typename Op, typename ParentCG>
  __device__ bool insert_or_apply_impl(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> group, Value value, Op op)
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);

    auto const val            = ref_.impl_.heterogeneous_value(value);
    auto const key            = ref_.impl_.extract_key(val);
    auto const probing_scheme = ref_.impl_.probing_scheme();
    auto storage_ref          = ref_.impl_.storage_ref();
    auto probing_iter =
      probing_scheme.template make_iterator<bucket_size>(group, key, storage_ref.extent());
    auto const init_idx    = *probing_iter;
    auto const empty_value = ref_.empty_value_sentinel();

    // wait for payload only when init != sentinel and insert strategy is not `packed_cas`
    auto constexpr wait_for_payload = (not UseDirectApply) and (sizeof(value_type) > 8);

    while (true) {
      auto const bucket_slots = storage_ref[*probing_iter];

      auto const [state, intra_bucket_index] = [&]() {
        auto res = detail::equal_result::UNEQUAL;
        for (auto i = 0; i < bucket_size; ++i) {
          res =
            ref_.impl_.predicate_.template operator()<is_insert::YES>(key, bucket_slots[i].first);
          if (res != detail::equal_result::UNEQUAL) {
            return detail::bucket_probing_results{res, i};
          }
        }
        // returns dummy index `-1` for UNEQUAL
        return detail::bucket_probing_results{res, -1};
      }();

      auto* slot_ptr = ref_.impl_.get_slot_ptr(*probing_iter, intra_bucket_index);

      auto const group_contains_equal = group.ballot(state == detail::equal_result::EQUAL);
      if (group_contains_equal) {
        auto const src_lane = __ffs(group_contains_equal) - 1;
        if (group.thread_rank() == src_lane) {
          if constexpr (wait_for_payload) {
            ref_.impl_.wait_for_payload(slot_ptr->second, empty_value);
          }
          op(cuda::atomic_ref<T, Scope>{slot_ptr->second}, val.second);
        }
        return false;
      }

      auto const group_contains_available = group.ballot(state == detail::equal_result::AVAILABLE);
      if (group_contains_available) {
        auto const src_lane = __ffs(group_contains_available) - 1;
        auto const status   = [&, target_idx = intra_bucket_index]() {
          if (group.thread_rank() != src_lane) { return insert_result::CONTINUE; }
          return ref_.attempt_insert_or_apply<UseDirectApply>(
            slot_ptr, bucket_slots[target_idx], val, op);
        }();

        switch (group.shfl(status, src_lane)) {
          case insert_result::SUCCESS: return true;
          case insert_result::DUPLICATE: {
            if (group.thread_rank() == src_lane) {
              if constexpr (wait_for_payload) {
                ref_.impl_.wait_for_payload(slot_ptr->second, empty_value);
              }
              op(cuda::atomic_ref<T, Scope>{slot_ptr->second}, val.second);
            }
            return false;
          }
          default: continue;
        }
      } else {
        ++probing_iter;
        if (*probing_iter == init_idx) { return false; }
      }
    }
  }

  /**
   * @brief Attempts to insert a key-value pair `{k, v}` if it's not present in the map. Otherwise,
   * applies `Op` binary function to the mapped_type corresponding to the key `k` and the value `v`.
   *
   * @tparam UseDirectApply Boolean condition which enables direct apply `op` instead of
   * wait_for_payload with an atomic_store on an empty slot when insert strategy is not
   * `packed_cas`.
   * @tparam Value Input type which is implicitly convertible to 'value_type'
   * @tparam Op Callable type which is used as apply operation and can be
   *   called with arguments as Op(cuda::atomic_ref<T, Scope>, T). Op strictly must
   *   have this signature to atomically apply the operation.
   *
   * @param address Pointer to the slot in memory
   * @param expected Element to compare against
   * @param desired Element to insert
   * @param op The callable object to perform binary operation between existing value at the slot
   *  and the element to insert.
   */
  template <bool UseDirectApply, typename Value, typename Op>
  [[nodiscard]] __device__ insert_result
  attempt_insert_or_apply(value_type* address, value_type expected, Value desired, Op op) noexcept
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);

    if constexpr (sizeof(value_type) <= 8) {
      return ref_.impl_.packed_cas(address, expected, desired);  // no need to wait for payload
    } else {
      using mapped_type = T;

      cuda::atomic_ref<key_type, Scope> key_ref(address->first);
      auto expected_key  = expected.first;
      auto const success = key_ref.compare_exchange_strong(
        expected_key, static_cast<key_type>(desired.first), cuda::memory_order_relaxed);

      // if key success
      if (success) {
        cuda::atomic_ref<mapped_type, Scope> payload_ref(address->second);
        // if init values == sentinel then directly apply the `op`
        if constexpr (UseDirectApply) {
          op(payload_ref, desired.second);
        } else {
          payload_ref.store(desired.second, cuda::memory_order_relaxed);
        }
        return insert_result::SUCCESS;
      }

      // Our key was already present in the slot, so our key is a duplicate
      // Shouldn't use `predicate` operator directly since it includes a redundant bitwise compare
      if (ref_.impl_.predicate_.equal_to(desired.first, expected_key) ==
          detail::equal_result::EQUAL) {
        return insert_result::DUPLICATE;
      }

      return insert_result::CONTINUE;
    }
  }
};

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::insert_and_find_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type     = typename base_type::value_type;
  using mapped_type    = T;
  using iterator       = typename base_type::iterator;
  using const_iterator = typename base_type::const_iterator;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Inserts the given element into the map.
   *
   * @note This API returns a pair consisting of an iterator to the inserted element (or to the
   * element that prevented the insertion) and a `bool` denoting whether the insertion took place or
   * not.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param value The element to insert
   *
   * @return a pair consisting of an iterator to the element and a bool indicating whether the
   * insertion is successful or not.
   */
  template <typename Value>
  __device__ cuda::std::pair<iterator, bool> insert_and_find(Value value) noexcept
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);
    return ref_.impl_.insert_and_find(value);
  }

  /**
   * @brief Inserts the given element into the map.
   *
   * @note This API returns a pair consisting of an iterator to the inserted element (or to the
   * element that prevented the insertion) and a `bool` denoting whether the insertion took place or
   * not.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group insert_and_find
   * @param value The element to insert
   *
   * @return a pair consisting of an iterator to the element and a bool indicating whether the
   * insertion is successful or not.
   */
  template <typename Value, typename ParentCG>
  __device__ cuda::std::pair<iterator, bool> insert_and_find(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> group, Value value) noexcept
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);
    return ref_.impl_.insert_and_find(group, value);
  }
};

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::erase_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type = typename base_type::value_type;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Erases an element.
   *
   * @tparam ProbeKey Input key type which is convertible to 'key_type'
   *
   * @param key The element to erase
   *
   * @return True if the given element is successfully erased
   */
  template <typename ProbeKey>
  __device__ bool erase(ProbeKey key) noexcept
  {
    ref_type& ref_ = static_cast<ref_type&>(*this);
    return ref_.impl_.erase(key);
  }

  /**
   * @brief Erases an element.
   *
   * @tparam ProbeKey Input key type which is convertible to 'key_type'
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group insert
   * @param key The element to erase
   *
   * @return True if the given element is successfully erased
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ bool erase(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                        ProbeKey key) noexcept
  {
    auto& ref_ = static_cast<ref_type&>(*this);
    return ref_.impl_.erase(group, key);
  }
};

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::contains_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type = typename base_type::value_type;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Indicates whether the probe key `key` was inserted into the container.
   *
   * @note If the probe key `key` was inserted into the container, returns
   * true. Otherwise, returns false.
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
    // CRTP: cast `this` to the actual ref type
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.contains(key);
  }

  /**
   * @brief Indicates whether the probe key `key` was inserted into the container.
   *
   * @note If the probe key `key` was inserted into the container, returns
   * true. Otherwise, returns false.
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
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::find_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type     = typename base_type::value_type;
  using iterator       = typename base_type::iterator;
  using const_iterator = typename base_type::const_iterator;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Finds an element in the map with key equivalent to the probe key.
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
  [[nodiscard]] __device__ iterator find(ProbeKey key) const noexcept
  {
    // CRTP: cast `this` to the actual ref type
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.find(key);
  }

  /**
   * @brief Finds an element in the map with key equivalent to the probe key.
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
  [[nodiscard]] __device__ iterator
  find(cooperative_groups::thread_block_tile<cg_size, ParentCG> group, ProbeKey key) const noexcept
  {
    auto const& ref_ = static_cast<ref_type const&>(*this);
    return ref_.impl_.find(group, key);
  }
};

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::for_each_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type     = typename base_type::value_type;
  using iterator       = typename base_type::iterator;
  using const_iterator = typename base_type::const_iterator;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief For a given key, applies the function object `callback_op` to its match found in the
   * container.
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @tparam ProbeKey Probe key type
   * @tparam CallbackOp Type of unary callback function object
   *
   * @param key The key to search for
   * @param callback_op Function to apply to the copy of the matched key-value pair
   */
  template <class ProbeKey, class CallbackOp>
  __device__ void for_each(ProbeKey key, CallbackOp&& callback_op) const noexcept
  {
    // CRTP: cast `this` to the actual ref type
    auto const& ref_ = static_cast<ref_type const&>(*this);
    ref_.impl_.for_each(key, cuda::std::forward<CallbackOp>(callback_op));
  }

  /**
   * @brief For a given key, applies the function object `callback_op` to its match found in the
   * container.
   *
   * @note This function uses cooperative group semantics, meaning that any thread may call the
   * callback if it finds a matching key-value pair.
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @note Synchronizing `group` within `callback_op` is undefined behavior.
   *
   * @tparam ProbeKey Probe key type
   * @tparam CallbackOp Type of unary callback function object
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform this operation
   * @param key The key to search for
   * @param callback_op Function to apply to the copy of the matched key-value pair
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
};

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::count_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
  using value_type = typename base_type::value_type;
  using size_type  = typename base_type::size_type;

  static constexpr auto cg_size     = base_type::cg_size;
  static constexpr auto bucket_size = base_type::bucket_size;

 public:
  /**
   * @brief Counts the occurrence of a given key contained in map
   *
   * @tparam ProbeKey Input type
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
   * @brief Counts the occurrence of a given key contained in map
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

template <typename Key,
          typename T,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class operator_impl<
  op::retrieve_tag,
  static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>> {
  using base_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef>;
  using ref_type = static_map_ref<Key, T, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>;
  using key_type = typename base_type::key_type;
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
   * If key `k = *(first + i)` exists in the container, copies `k` to `
   * output_probe` and associated slot content to `output_match`, respectively. The output order is
   * unspecified.
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
  template <int BlockSize,
            class InputProbeIt,
            class OutputProbeIt,
            class OutputMatchIt,
            class AtomicCounter>
  __device__ void retrieve(cooperative_groups::thread_block block,
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

}  // namespace detail
}  // namespace cuco