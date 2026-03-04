/*
 * Copyright (c) 2021-2024, NVIDIA CORPORATION.
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

#include <cuco/detail/utils.cuh>
#include <cuco/pair.cuh>

#include <cuda/std/atomic>

#include <cooperative_groups.h>

#include <utility>

namespace cuco::legacy::detail {

/**
 * @brief Base class of public probe sequence. This class should not be used directly.
 *
 * @tparam CGSize Size of CUDA Cooperative Groups
 */
template <uint32_t CGSize>
class probe_sequence_base {
 protected:
  /**
   * @brief Returns the size of the CUDA cooperative thread group.
   */
  static constexpr std::size_t cg_size = CGSize;

  /**
   * @brief Returns the number of elements loaded with each vector load.
   *
   * @return The number of elements loaded with each vector load
   */
  static __host__ __device__ constexpr uint32_t vector_width() noexcept { return 2u; }
};

/**
 * @brief Base class of probe sequence implementation.
 *
 * Hash map operations are generally memory-bandwidth bound. A vector-load loads two consecutive
 * slots instead of one to fully utilize the 16B memory load supported by SASS/hardware thus
 * improve memory throughput. This method (flagged by `uses_vector_load` logic) is implicitly
 * applied to all hash map operations (e.g. `insert`, `count`, and `retrieve`, etc.) when pairs
 * are packable (see `cuco::detail::is_packable` logic).
 *
 * @tparam Key Type used for keys
 * @tparam Value Type of the mapped values
 * @tparam Scope The scope in which multimap operations will be performed by
 * individual threads
 * @tparam VectorWidth Length of vector load
 * @tparam CGSize Size of CUDA Cooperative Groups
 */
template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          uint32_t VectorWidth,
          uint32_t CGSize>
class probe_sequence_impl_base {
 protected:
  using value_type         = cuco::pair<Key, Value>;            ///< Type of key/value pairs
  using key_type           = Key;                               ///< Key type
  using mapped_type        = Value;                             ///< Type of mapped values
  using atomic_key_type    = cuda::atomic<key_type, Scope>;     ///< Type of atomic keys
  using atomic_mapped_type = cuda::atomic<mapped_type, Scope>;  ///< Type of atomic mapped values
  /// Pair type of atomic key and atomic mapped value
  using pair_atomic_type = cuco::pair<atomic_key_type, atomic_mapped_type>;
  /// Type of the forward iterator to `pair_atomic_type`
  using iterator = pair_atomic_type*;
  /// Type of the forward iterator to `const pair_atomic_type`
  using const_iterator = pair_atomic_type const*;

  /**
   * @brief Returns the number of elements loaded with each vector-load.
   */
  static constexpr uint32_t vector_width = VectorWidth;

  /**
   * @brief Returns the size of the CUDA cooperative thread group.
   */
  static constexpr std::size_t cg_size = CGSize;

  /**
   * @brief Indicates if vector-load is used.
   *
   * Users have no explicit control on whether vector-load is used.
   *
   * @return Boolean indicating if vector-load is used.
   */
  __host__ __device__ static constexpr bool uses_vector_load() noexcept
  {
    return cuco::detail::is_packable<value_type>();
  }

  /**
   * @brief Constructs a probe sequence based on the given hash map features.
   *
   * @param slots Pointer to beginning of the hash map slots
   * @param capacity Capacity of the hash map
   */
  __host__ __device__ explicit probe_sequence_impl_base(iterator slots, std::size_t capacity)
    : slots_{slots}, capacity_{capacity}
  {
  }

 public:
  /**
   * @brief Returns the capacity of the hash map.
   *
   * @return The capacity of the hash map
   */
  __host__ __device__ __forceinline__ std::size_t get_capacity() const noexcept
  {
    return capacity_;
  }

  /**
   * @brief Returns slots array.
   *
   * @return Slots array
   */
  __device__ __forceinline__ iterator get_slots() noexcept { return slots_; }

  /**
   * @brief Returns slots array.
   *
   * @return Slots array
   */
  __device__ __forceinline__ const_iterator get_slots() const noexcept { return slots_; }

 protected:
  iterator slots_;              ///< Pointer to beginning of the hash map slots
  const std::size_t capacity_;  ///< Total number of slots
};  // class probe_sequence_impl_base

/**
 * @brief Cooperative Groups based Linear probing scheme.
 *
 * @tparam Key Type used for keys
 * @tparam Value Type of the mapped values
 * @tparam Scope The scope in which multimap operations will be performed by
 * individual threads
 * @tparam VectorWidth Length of vector load
 * @tparam CGSize Size of CUDA Cooperative Groups
 * @tparam Hash Unary callable type
 */
template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          uint32_t VectorWidth,
          int32_t CGSize,
          typename Hash>
class linear_probing_impl
  : public probe_sequence_impl_base<Key, Value, Scope, VectorWidth, CGSize> {
 public:
  using probe_sequence_impl_base_type =
    probe_sequence_impl_base<Key, Value, Scope, VectorWidth, CGSize>;
  using value_type         = typename probe_sequence_impl_base_type::value_type;
  using key_type           = typename probe_sequence_impl_base_type::key_type;
  using mapped_type        = typename probe_sequence_impl_base_type::mapped_type;
  using atomic_key_type    = typename probe_sequence_impl_base_type::atomic_key_type;
  using atomic_mapped_type = typename probe_sequence_impl_base_type::atomic_mapped_type;
  using pair_atomic_type   = typename probe_sequence_impl_base_type::pair_atomic_type;
  using iterator           = typename probe_sequence_impl_base_type::iterator;
  using const_iterator     = typename probe_sequence_impl_base_type::const_iterator;

  using probe_sequence_impl_base_type::capacity_;
  using probe_sequence_impl_base_type::cg_size;
  using probe_sequence_impl_base_type::slots_;
  using probe_sequence_impl_base_type::uses_vector_load;
  using probe_sequence_impl_base_type::vector_width;

  /**
   * @brief Constructs a linear probing scheme based on the given hash map features.
   *
   * @param slots Pointer to beginning of the hash map slots
   * @param capacity Capacity of the hash map
   * @param hash Unary function to hash each key
   */
  __host__ __device__ explicit linear_probing_impl(iterator slots, std::size_t capacity)
    : probe_sequence_impl_base_type{slots, capacity}, hash_{Hash{}}
  {
  }

  /**
   * @brief Returns the initial slot for a given key `k`.
   *
   * If vector-load is enabled, the return slot is always even to avoid illegal memory access.
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g the Cooperative Group for which the initial slot is needed
   * @param k The key to get the slot for
   * @return Pointer to the initial slot for `k`
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ __forceinline__ iterator initial_slot(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> g, ProbeKey const& k) noexcept
  {
    return const_cast<iterator>(cuda::std::as_const(*this).initial_slot(g, k));
  }

  /**
   * @brief Returns the initial slot for a given key `k`.
   *
   * If vector-load is enabled, the return slot is always even to avoid illegal memory access.
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g the Cooperative Group for which the initial slot is needed
   * @param k The key to get the slot for
   * @return Pointer to the initial slot for `k`
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ __forceinline__ const_iterator initial_slot(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> g, ProbeKey const& k) const noexcept
  {
    auto const hash_value = [&]() {
      auto const tmp = hash_(k);
      if constexpr (uses_vector_load()) {
        // initial hash value is always even
        return tmp + tmp % 2;
      }
      if constexpr (not uses_vector_load()) { return tmp; }
    }();

    auto const offset = [&]() {
      if constexpr (uses_vector_load()) { return g.thread_rank() * vector_width; }
      if constexpr (not uses_vector_load()) { return g.thread_rank(); }
    }();

    // Each CG accesses to a bucket of (`cg_size` * `vector_width`)
    // slots if vector-load is used or `cg_size` slots otherwise
    return &slots_[(hash_value + offset) % capacity_];
  }

  /**
   * @brief Given a slot `s`, returns the next slot.
   *
   * If `s` is the last slot, wraps back around to the first slot.
   *
   * @param s The slot to advance
   * @return The next slot after `s`
   */
  __device__ __forceinline__ iterator next_slot(iterator s) noexcept
  {
    return const_cast<iterator>(cuda::std::as_const(*this).next_slot(s));
  }

  /**
   * @brief Given a slot `s`, returns the next slot.
   *
   * If `s` is the last slot, wraps back around to the first slot.
   *
   * @param s The slot to advance
   * @return The next slot after `s`
   */
  __device__ __forceinline__ const_iterator next_slot(const_iterator s) const noexcept
  {
    std::size_t index = s - slots_;
    std::size_t offset;
    if constexpr (uses_vector_load()) {
      offset = cg_size * vector_width;
    } else {
      offset = cg_size;
    }
    return &slots_[(index + offset) % capacity_];
  }

 private:
  Hash hash_;  ///< The unary callable used to hash the key
};  // class linear_probing

/**
 * @brief Cooperative Groups based double hashing scheme.
 *
 * Default probe sequence for `cuco::static_multimap`. Double hashing shows superior
 * performance when dealing with high multiplicty and/or high occupancy use cases. Performance
 * hints:
 * - `CGSize` = 1 or 2 when hash map is small (10'000'000 or less), 4 or 8 otherwise.
 *
 * `Hash1` and `Hash2` should be callable object type.
 *
 * @tparam Key Type used for keys
 * @tparam Value Type of the mapped values
 * @tparam Scope The scope in which multimap operations will be performed by
 * individual threads
 * @tparam VectorWidth Length of vector load
 * @tparam CGSize Size of CUDA Cooperative Groups
 * @tparam Hash1 Unary callable type
 * @tparam Hash2 Unary callable type
 */
template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          uint32_t VectorWidth,
          uint32_t CGSize,
          typename Hash1,
          typename Hash2>
class double_hashing_impl
  : public probe_sequence_impl_base<Key, Value, Scope, VectorWidth, CGSize> {
 public:
  using probe_sequence_impl_base_type =
    probe_sequence_impl_base<Key, Value, Scope, VectorWidth, CGSize>;
  using value_type         = typename probe_sequence_impl_base_type::value_type;
  using key_type           = typename probe_sequence_impl_base_type::key_type;
  using mapped_type        = typename probe_sequence_impl_base_type::mapped_type;
  using atomic_key_type    = typename probe_sequence_impl_base_type::atomic_key_type;
  using atomic_mapped_type = typename probe_sequence_impl_base_type::atomic_mapped_type;
  using pair_atomic_type   = typename probe_sequence_impl_base_type::pair_atomic_type;
  using iterator           = typename probe_sequence_impl_base_type::iterator;
  using const_iterator     = typename probe_sequence_impl_base_type::const_iterator;

  using probe_sequence_impl_base_type::capacity_;
  using probe_sequence_impl_base_type::cg_size;
  using probe_sequence_impl_base_type::slots_;
  using probe_sequence_impl_base_type::uses_vector_load;
  using probe_sequence_impl_base_type::vector_width;

  /**
   * @brief Constructs a double hashing scheme based on the given hash map features.
   *
   * `hash2` takes a different seed to reduce the chance of secondary clustering.
   *
   * @param slots Pointer to beginning of the hash map slots
   * @param capacity Capacity of the hash map
   * @param hash1 First hasher to hash each key
   * @param hash2 Second hasher to determine step size
   */
  __host__ __device__ explicit double_hashing_impl(iterator slots, std::size_t capacity)
    : probe_sequence_impl_base_type{slots, capacity},
      hash1_{Hash1{}},
      hash2_{Hash2{1}},
      step_size_{}
  {
  }

  /**
   * @brief Returns the initial slot for a given key `k`.
   *
   * If vector-load is enabled, the return slot is always a multiple of (`cg_size` * `vector_width`)
   * to avoid illegal memory access.
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g the Cooperative Group for which the initial slot is needed
   * @param k The key to get the slot for
   * @return Pointer to the initial slot for `k`
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ __forceinline__ iterator initial_slot(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> g, ProbeKey const& k) noexcept
  {
    return const_cast<iterator>(cuda::std::as_const(*this).initial_slot(g, k));
  }

  /**
   * @brief Returns the initial slot for a given key `k`.
   *
   * If vector-load is enabled, the return slot is always a multiple of (`cg_size` * `vector_width`)
   * to avoid illegal memory access.
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g the Cooperative Group for which the initial slot is needed
   * @param k The key to get the slot for
   * @return Pointer to the initial slot for `k`
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ __forceinline__ const_iterator initial_slot(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> g, ProbeKey const& k) const noexcept
  {
    std::size_t index;
    auto const hash_value = hash1_(k);
    if constexpr (uses_vector_load()) {
      // step size in range [1, prime - 1] * cg_size * vector_width
      step_size_ =
        (hash2_(k) % (capacity_ / (cg_size * vector_width) - 1) + 1) * cg_size * vector_width;
      index = hash_value % (capacity_ / (cg_size * vector_width)) * cg_size * vector_width +
              g.thread_rank() * vector_width;
    } else {
      // step size in range [1, prime - 1] * cg_size
      step_size_ = (hash2_(k) % (capacity_ / cg_size - 1) + 1) * cg_size;
      index      = (hash_value + g.thread_rank()) % capacity_;
    }
    return slots_ + index;
  }

  /**
   * @brief Given a slot `s`, returns the next slot.
   *
   * If `s` is the last slot, wraps back around to the first slot.
   *
   * @param s The slot to advance
   * @return The next slot after `s`
   */
  __device__ __forceinline__ iterator next_slot(iterator s) noexcept
  {
    return const_cast<iterator>(cuda::std::as_const(*this).next_slot(s));
  }

  /**
   * @brief Given a slot `s`, returns the next slot.
   *
   * If `s` is the last slot, wraps back around to the first slot.
   *
   * @param s The slot to advance
   * @return The next slot after `s`
   */
  __device__ __forceinline__ const_iterator next_slot(const_iterator s) const noexcept
  {
    std::size_t index = s - slots_;
    return &slots_[(index + step_size_) % capacity_];
  }

 private:
  Hash1 hash1_;                    ///< The first unary callable used to hash the key
  Hash2 hash2_;                    ///< The second unary callable used to determine step size
  mutable std::size_t step_size_;  ///< The step stride when searching for the next slot
};  // class double_hashing

/**
 * @brief Probe sequence used internally by hash map.
 *
 * @tparam ProbeImpl Type of probe sequence implementation
 * @tparam Key Type used for keys
 * @tparam Value Type of the mapped values
 * @tparam Scope The scope in which multimap operations will be performed by
 * individual threads
 */
template <typename ProbeImpl, typename Key, typename Value, cuda::thread_scope Scope>
class probe_sequence : public ProbeImpl::template impl<Key, Value, Scope> {
 public:
  using impl_type =
    typename ProbeImpl::template impl<Key, Value, Scope>;  ///< Type of implementation details

  /**
   * @brief Constructs a probe sequence based on the given hash map features.
   *
   * @param slots Pointer to beginning of the hash map slots
   * @param capacity Capacity of the hash map
   */
  __host__ __device__ explicit probe_sequence(typename impl_type::iterator slots,
                                              std::size_t capacity)
    : impl_type{slots, capacity}
  {
  }
};  // class probe_sequence

}  // namespace cuco::legacy::detail
