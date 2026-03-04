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

#include <cuco/detail/open_addressing/open_addressing_ref_impl.cuh>
#include <cuco/hash_functions.cuh>
#include <cuco/operator.hpp>
#include <cuco/probing_scheme.cuh>
#include <cuco/storage.cuh>
#include <cuco/types.cuh>
#include <cuco/utility/cuda_thread_scope.cuh>

#include <cuda/std/atomic>

#include <memory>

namespace cuco {

/**
 * @brief Device non-owning "ref" of `static_multiset` that can be used in device code to perform
 * arbitrary operations defined in `include/cuco/operator.hpp`
 *
 * @note Concurrent modify and lookup will be supported if both kinds of operators are specified
 * during the ref construction.
 * @note cuCollections data structures always place the slot keys on the left-hand
 * side when invoking the key comparison predicate.
 * @note Ref types are trivially-copyable and are intended to be passed by value.
 * @note `ProbingScheme::cg_size` indicates how many threads are used to handle one independent
 * device operation. `cg_size == 1` uses the scalar (or non-CG) code paths.
 *
 * @throw If the size of the given key type is larger than 8 bytes
 * @throw If the given key type doesn't have unique object representations, i.e.,
 * `cuco::bitwise_comparable_v<Key> == false`
 * @throw If the probing scheme type is not inherited from `cuco::detail::probing_scheme_base`
 *
 * @tparam Key Type used for keys. Requires `cuco::is_bitwise_comparable_v<Key>` returning true
 * @tparam Scope The scope in which operations will be performed by individual threads.
 * @tparam KeyEqual Binary callable type used to compare two keys for equality
 * @tparam ProbingScheme Probing scheme (see `include/cuco/probing_scheme.cuh` for options)
 * @tparam StorageRef Storage ref type
 * @tparam Operators Device operator options defined in `include/cuco/operator.hpp`
 */
template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          typename... Operators>
class static_multiset_ref
  : public detail::operator_impl<
      Operators,
      static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, Operators...>>... {
  /// Flag indicating whether duplicate keys are allowed or not
  static constexpr auto allows_duplicates = true;

  /// Implementation type
  using impl_type = detail::
    open_addressing_ref_impl<Key, Scope, KeyEqual, ProbingScheme, StorageRef, allows_duplicates>;

 public:
  using key_type            = Key;                                     ///< Key Type
  using probing_scheme_type = ProbingScheme;                           ///< Type of probing scheme
  using hasher              = typename probing_scheme_type::hasher;    ///< Hash function type
  using storage_ref_type    = StorageRef;                              ///< Type of storage ref
  using bucket_type         = typename storage_ref_type::bucket_type;  ///< Bucket type
  using value_type          = typename storage_ref_type::value_type;   ///< Storage element type
  using extent_type         = typename storage_ref_type::extent_type;  ///< Extent type
  using size_type           = typename storage_ref_type::size_type;    ///< Probing scheme size type
  using key_equal           = KeyEqual;  ///< Type of key equality binary callable
  using iterator            = typename storage_ref_type::iterator;   ///< Slot iterator type
  using const_iterator = typename storage_ref_type::const_iterator;  ///< Const slot iterator type

  static constexpr auto cg_size = probing_scheme_type::cg_size;  ///< Cooperative group size
  static constexpr auto bucket_size =
    storage_ref_type::bucket_size;  ///< Number of elements handled per bucket
  static constexpr auto thread_scope = impl_type::thread_scope;  ///< CUDA thread scope

  /**
   * @brief Constructs static_multiset_ref
   *
   * @param empty_key_sentinel Sentinel indicating empty key
   * @param predicate Key equality binary callable
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage_ref Non-owning ref of slot storage
   */
  __host__ __device__ explicit constexpr static_multiset_ref(
    cuco::empty_key<Key> empty_key_sentinel,
    KeyEqual const& predicate,
    ProbingScheme const& probing_scheme,
    cuda_thread_scope<Scope> scope,
    StorageRef storage_ref) noexcept;

  /**
   * @brief Constructs static_multiset_ref
   *
   * @param empty_key_sentinel Sentinel indicating empty key
   * @param erased_key_sentinel Sentinel indicating erased key
   * @param predicate Key equality binary callable
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage_ref Non-owning ref of slot storage
   */
  __host__ __device__ explicit constexpr static_multiset_ref(
    cuco::empty_key<Key> empty_key_sentinel,
    cuco::erased_key<Key> erased_key_sentinel,
    KeyEqual const& predicate,
    ProbingScheme const& probing_scheme,
    cuda_thread_scope<Scope> scope,
    StorageRef storage_ref) noexcept;

  /**
   * @brief Operator-agnostic move constructor.
   *
   * @tparam OtherOperators Operator set of the `other` object
   *
   * @param other Object to construct `*this` from
   */
  template <typename... OtherOperators>
  __host__ __device__ explicit constexpr static_multiset_ref(
    static_multiset_ref<Key, Scope, KeyEqual, ProbingScheme, StorageRef, OtherOperators...>&&
      other) noexcept;

  /**
   * @brief Gets the maximum number of elements the container can hold.
   *
   * @return The maximum number of elements the container can hold
   */
  [[nodiscard]] __host__ __device__ constexpr auto capacity() const noexcept;

  /**
   * @brief Gets the extent of the current storage.
   *
   * @return The bucket extent.
   */
  [[nodiscard]] __host__ __device__ constexpr extent_type extent() const noexcept;

  /**
   * @brief Gets the bucket extent of the current storage.
   *
   * @return The bucket extent.
   */
  [[nodiscard]] __host__ __device__ constexpr extent_type bucket_extent() const noexcept;

  /**
   * @brief Gets the sentinel value used to represent an empty key slot.
   *
   * @return The sentinel value used to represent an empty key slot
   */
  [[nodiscard]] __host__ __device__ constexpr key_type empty_key_sentinel() const noexcept;

  /**
   * @brief Gets the sentinel value used to represent an erased key slot.
   *
   * @return The sentinel value used to represent an erased key slot
   */
  [[nodiscard]] __host__ __device__ constexpr key_type erased_key_sentinel() const noexcept;

  /**
   * @brief Gets the key comparator.
   *
   * @return The comparator used to compare keys
   */
  [[nodiscard]] __host__ __device__ constexpr key_equal key_eq() const noexcept;

  /**
   * @brief Gets the function(s) used to hash keys
   *
   * @return The function(s) used to hash keys
   */
  [[nodiscard]] __host__ __device__ constexpr hasher hash_function() const noexcept;

  /**
   * @brief Returns a const_iterator to one past the last slot.
   *
   * @return A const_iterator to one past the last slot
   */
  [[nodiscard]] __host__ __device__ constexpr const_iterator end() const noexcept;

  /**
   * @brief Returns an iterator to one past the last slot.
   *
   * @return An iterator to one past the last slot
   */
  [[nodiscard]] __host__ __device__ constexpr iterator end() noexcept;

  /**
   * @brief Gets the non-owning storage ref.
   *
   * @return The non-owning storage ref of the container
   */
  [[nodiscard]] __host__ __device__ constexpr auto storage_ref() const noexcept;

  /**
   * @brief Gets the probing scheme.
   *
   * @return The probing scheme used for the container
   */
  [[nodiscard]] __host__ __device__ constexpr auto probing_scheme() const noexcept;

  /**
   * @brief Creates a copy of the current non-owning reference using the given operators
   *
   * @tparam NewOperators List of `cuco::op::*_tag` types
   *
   * @param ops List of operators, e.g., `cuco::op::insert`
   *
   * @return Copy of the current device ref
   */
  template <typename... NewOperators>
  [[nodiscard]] __host__ __device__ constexpr auto rebind_operators(
    NewOperators... ops) const noexcept;

  /**
   * @brief Makes a copy of the current device reference with the given key comparator
   *
   * @tparam NewKeyEqual The new key equal type
   *
   * @param key_equal New key comparator
   *
   * @return Copy of the current device ref
   */
  template <typename NewKeyEqual>
  [[nodiscard]] __host__ __device__ constexpr auto rebind_key_eq(
    NewKeyEqual const& key_equal) const noexcept;

  /**
   * @brief Makes a copy of the current device reference with the given hasher
   *
   * @tparam NewHash The new hasher type
   *
   * @param hash New hasher
   *
   * @return Copy of the current device ref
   */
  template <typename NewHash>
  [[nodiscard]] __host__ __device__ constexpr auto rebind_hash_function(NewHash const& hash) const;

  /**
   * @brief Makes a copy of the current device reference using non-owned memory
   *
   * This function is intended to be used to create shared memory copies of small static sets,
   * although global memory can be used as well.
   *
   * @note This function synchronizes the group `tile`.
   * @note By-default the thread scope of the copy will be the same as the scope of the parent ref.
   *
   * @tparam CG The type of the cooperative thread group
   * @tparam NewScope The thread scope of the newly created device ref
   *
   * @param tile The cooperative thread group used to copy the data structure
   * @param memory_to_use Array large enough to support `capacity` elements. Object does not take
   * the ownership of the memory
   * @param scope The thread scope of the newly created device ref
   *
   * @return Copy of the current device ref
   */
  template <typename CG, cuda::thread_scope NewScope = thread_scope>
  [[nodiscard]] __device__ constexpr auto make_copy(
    CG tile,
    bucket_type* const memory_to_use,
    cuda_thread_scope<NewScope> scope = {}) const noexcept;

  /**
   * @brief Initializes the set storage using the threads in the group `tile`.
   *
   * @note This function synchronizes the group `tile`.
   *
   * @tparam CG The type of the cooperative thread group
   *
   * @param tile The cooperative thread group used to initialize the set
   */
  template <typename CG>
  __device__ constexpr void initialize(CG tile) noexcept;

 private:
  impl_type impl_;

  // Mixins need to be friends with this class in order to access private members
  template <typename Op, typename Ref>
  friend class detail::operator_impl;

  // Refs with other operator sets need to be friends too
  template <typename Key_,
            cuda::thread_scope Scope_,
            typename KeyEqual_,
            typename ProbingScheme_,
            typename StorageRef_,
            typename... Operators_>
  friend class static_multiset_ref;
};

}  // namespace cuco

#include <cuco/detail/static_multiset/static_multiset_ref.inl>
