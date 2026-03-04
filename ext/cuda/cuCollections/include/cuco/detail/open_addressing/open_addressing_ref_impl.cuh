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

#include <cuco/detail/equal_wrapper.cuh>
#include <cuco/detail/probing_scheme/probing_scheme_base.cuh>
#include <cuco/detail/utility/cuda.cuh>
#include <cuco/extent.cuh>
#include <cuco/pair.cuh>
#include <cuco/probing_scheme.cuh>

#include <cuda/atomic>
#include <cuda/std/cstdint>
#include <cuda/std/functional>
#include <cuda/std/iterator>
#include <cuda/std/type_traits>
#include <thrust/execution_policy.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/logical.h>
#include <thrust/reduce.h>
#if defined(CUCO_HAS_CUDA_BARRIER)
#include <cuda/barrier>
#endif

#include <cooperative_groups.h>

namespace cuco {
namespace detail {

/// Three-way insert result enum
enum class insert_result : cuda::std::int32_t { CONTINUE = 0, SUCCESS = 1, DUPLICATE = 2 };

/**
 * @brief Helper struct to store intermediate bucket probing results.
 */
struct bucket_probing_results {
  detail::equal_result state_;             ///< Equal result
  cuda::std::int32_t intra_bucket_index_;  ///< Intra-bucket index

  /**
   * @brief Constructs bucket_probing_results.
   *
   * @param state The three way equality result
   * @param index Intra-bucket index
   */
  __device__ explicit constexpr bucket_probing_results(detail::equal_result state,
                                                       cuda::std::int32_t index) noexcept
    : state_{state}, intra_bucket_index_{index}
  {
  }
};

/**
 * @brief Common device non-owning "ref" implementation class.
 *
 * @note This class should NOT be used directly.
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
 * @tparam AllowsDuplicates Flag indicating whether duplicate keys are allowed or not
 */
template <typename Key,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename StorageRef,
          bool AllowsDuplicates>
class open_addressing_ref_impl {
  static_assert(sizeof(Key) <= 8, "Container does not support key types larger than 8 bytes.");

  static_assert(
    cuco::is_bitwise_comparable_v<Key>,
    "Key type must have unique object representations or have been explicitly declared as safe for "
    "bitwise comparison via specialization of cuco::is_bitwise_comparable_v<Key>.");

  static_assert(cuda::std::is_base_of_v<cuco::detail::probing_scheme_base<ProbingScheme::cg_size>,
                                        ProbingScheme>,
                "ProbingScheme must inherit from cuco::detail::probing_scheme_base");

  /// Determines if the container is a key/value or key-only store
  static constexpr auto has_payload =
    not cuda::std::is_same_v<Key, typename StorageRef::value_type>;

  /// Flag indicating whether duplicate keys are allowed or not
  static constexpr auto allows_duplicates = AllowsDuplicates;

  // TODO: how to re-enable this check?
  // static_assert(is_bucket_extent_v<typename StorageRef::extent_type>,
  // "Extent is not a valid cuco::bucket_extent");

 public:
  using key_type            = Key;                                     ///< Key type
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
    storage_ref_type::bucket_size;             ///< Number of elements handled per bucket
  static constexpr auto thread_scope = Scope;  ///< CUDA thread scope

  /**
   * @brief Constructs open_addressing_ref_impl.
   *
   * @param empty_slot_sentinel Sentinel indicating an empty slot
   * @param predicate Key equality binary callable
   * @param probing_scheme Probing scheme
   * @param storage_ref Non-owning ref of slot storage
   */
  __host__ __device__ explicit constexpr open_addressing_ref_impl(
    value_type empty_slot_sentinel,
    key_equal const& predicate,
    probing_scheme_type const& probing_scheme,
    storage_ref_type storage_ref) noexcept
    : empty_slot_sentinel_{empty_slot_sentinel},
      predicate_{
        this->extract_key(empty_slot_sentinel), this->extract_key(empty_slot_sentinel), predicate},
      probing_scheme_{probing_scheme},
      storage_ref_{storage_ref}
  {
  }

  /**
   * @brief Constructs open_addressing_ref_impl.
   *
   * @param empty_slot_sentinel Sentinel indicating an empty slot
   * @param erased_key_sentinel Sentinel indicating an erased key
   * @param predicate Key equality binary callable
   * @param probing_scheme Probing scheme
   * @param storage_ref Non-owning ref of slot storage
   */
  __host__ __device__ explicit constexpr open_addressing_ref_impl(
    value_type empty_slot_sentinel,
    key_type erased_key_sentinel,
    key_equal const& predicate,
    probing_scheme_type const& probing_scheme,
    storage_ref_type storage_ref) noexcept
    : empty_slot_sentinel_{empty_slot_sentinel},
      predicate_{this->extract_key(empty_slot_sentinel), erased_key_sentinel, predicate},
      probing_scheme_{probing_scheme},
      storage_ref_{storage_ref}
  {
  }

  /**
   * @brief Gets the sentinel value used to represent an empty key slot.
   *
   * @return The sentinel value used to represent an empty key slot
   */
  [[nodiscard]] __host__ __device__ constexpr key_type empty_key_sentinel() const noexcept
  {
    return this->predicate_.empty_sentinel_;
  }

  /**
   * @brief Gets the sentinel value used to represent an empty payload slot.
   *
   * @return The sentinel value used to represent an empty payload slot
   */
  template <bool Dummy = true, typename Enable = cuda::std::enable_if_t<has_payload and Dummy>>
  [[nodiscard]] __host__ __device__ constexpr auto empty_value_sentinel() const noexcept
  {
    return this->extract_payload(this->empty_slot_sentinel());
  }

  /**
   * @brief Gets the sentinel value used to represent an erased key slot.
   *
   * @return The sentinel value used to represent an erased key slot
   */
  [[nodiscard]] __host__ __device__ constexpr key_type erased_key_sentinel() const noexcept
  {
    return this->predicate_.erased_sentinel_;
  }

  /**
   * @brief Gets the sentinel used to represent an empty slot.
   *
   * @return The sentinel value used to represent an empty slot
   */
  [[nodiscard]] __host__ __device__ constexpr value_type empty_slot_sentinel() const noexcept
  {
    return empty_slot_sentinel_;
  }

  /**
   * @brief Returns the function that compares keys for equality.
   *
   * @return The key equality predicate
   */
  [[nodiscard]] __host__ __device__ constexpr detail::equal_wrapper<key_type, key_equal> predicate()
    const noexcept
  {
    return this->predicate_;
  }

  /**
   * @brief Gets the key comparator.
   *
   * @return The comparator used to compare keys
   */
  [[nodiscard]] __host__ __device__ constexpr key_equal key_eq() const noexcept
  {
    return this->predicate().equal_;
  }

  /**
   * @brief Gets the probing scheme.
   *
   * @return The probing scheme used for the container
   */
  [[nodiscard]] __host__ __device__ constexpr probing_scheme_type probing_scheme() const noexcept
  {
    return probing_scheme_;
  }

  /**
   * @brief Gets the function(s) used to hash keys
   *
   * @return The function(s) used to hash keys
   */
  [[nodiscard]] __host__ __device__ constexpr hasher hash_function() const noexcept
  {
    return this->probing_scheme().hash_function();
  }

  /**
   * @brief Gets the non-owning storage ref.
   *
   * @return The non-owning storage ref of the container
   */
  [[nodiscard]] __host__ __device__ constexpr storage_ref_type storage_ref() const noexcept
  {
    return storage_ref_;
  }

  /**
   * @brief Gets the maximum number of elements the container can hold.
   *
   * @return The maximum number of elements the container can hold
   */
  [[nodiscard]] __host__ __device__ constexpr auto capacity() const noexcept
  {
    return storage_ref_.capacity();
  }

  /**
   * @brief Gets the bucket extent of the current storage.
   *
   * @return The bucket extent.
   */
  [[nodiscard]] __host__ __device__ constexpr extent_type extent() const noexcept
  {
    return storage_ref_.extent();
  }

  /**
   * @brief Returns an iterator to one past the last slot.
   *
   * @return An iterator to one past the last slot
   */
  [[nodiscard]] __host__ __device__ constexpr iterator end() const noexcept
  {
    return storage_ref_.end();
  }

  /**
   * @brief Returns an iterator to one past the last slot.
   *
   * @return An iterator to one past the last slot
   */
  [[nodiscard]] __host__ __device__ constexpr iterator end() noexcept { return storage_ref_.end(); }

  /**
   * @brief Makes a copy of the current device reference using non-owned memory.
   *
   * This function is intended to be used to create shared memory copies of small static data
   * structures, although global memory can be used as well.
   *
   * @tparam CG The type of the cooperative thread group
   *
   * @param g The cooperative thread group used to copy the data structure
   * @param memory_to_use Array large enough to support `capacity` elements. Object does not take
   * the ownership of the memory
   */
  template <typename CG>
  __device__ void make_copy(CG g, value_type* const memory_to_use) const noexcept
  {
    auto const num_slots = this->capacity();
#if defined(CUCO_HAS_CUDA_BARRIER)
#pragma nv_diagnostic push
// Disables `barrier` initialization warning.
#pragma nv_diag_suppress static_var_with_dynamic_init
    __shared__ cuda::barrier<cuda::thread_scope::thread_scope_block> barrier;
#pragma nv_diagnostic pop
    if (g.thread_rank() == 0) { init(&barrier, g.size()); }
    g.sync();

    cuda::memcpy_async(
      g, memory_to_use, this->storage_ref().data(), sizeof(value_type) * num_slots, barrier);

    barrier.arrive_and_wait();
#else
    value_type const* const slots_ptr = this->storage_ref().data();
    for (size_type i = g.thread_rank(); i < num_slots; i += g.size()) {
      memory_to_use[i] = slots_ptr[i];
    }
    g.sync();
#endif
  }

  /**
   * @brief Initializes the container storage.
   *
   * @note This function synchronizes the group `tile`.
   *
   * @tparam CG The type of the cooperative thread group
   *
   * @param tile The cooperative thread group used to initialize the container
   */
  template <typename CG>
  __device__ constexpr void initialize(CG tile) noexcept
  {
    auto tid          = tile.thread_rank();
    auto const extent = static_cast<size_type>(this->extent());

    auto* const slots_ptr = this->storage_ref().data();
    while (tid < extent) {
      slots_ptr[tid] = this->empty_slot_sentinel();
      tid += tile.size();
    }

    tile.sync();
  }

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
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");

    auto const val = this->heterogeneous_value(value);
    auto const key = this->extract_key(val);

    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref_[*probing_iter];

      for (auto& slot_content : bucket_slots) {
        auto const eq_res = this->predicate_.template operator()<is_insert::YES>(
          key, this->extract_key(slot_content));

        if constexpr (not allows_duplicates) {
          // If the key is already in the container, return false
          if (eq_res == detail::equal_result::EQUAL) { return false; }
        }
        if (eq_res == detail::equal_result::AVAILABLE) {
          auto const intra_bucket_index = cuda::std::distance(bucket_slots.begin(), &slot_content);
          switch (attempt_insert(
            this->get_slot_ptr(*probing_iter, intra_bucket_index), slot_content, val)) {
            case insert_result::DUPLICATE: {
              if constexpr (allows_duplicates) {
                [[fallthrough]];
              } else {
                return false;
              }
            }
            case insert_result::CONTINUE: continue;
            case insert_result::SUCCESS: return true;
          }
        }
      }
      ++probing_iter;
      if (*probing_iter == init_idx) { return false; }
    }
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
  template <bool SupportsErase, typename Value, typename ParentCG>
  __device__ bool insert(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                         Value value) noexcept
  {
    auto const val = this->heterogeneous_value(value);
    auto const key = this->extract_key(val);
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(group, key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref_[*probing_iter];

      auto const [state, intra_bucket_index] = [&]() {
        for (auto i = 0; i < bucket_size; ++i) {
          switch (this->predicate_.template operator()<is_insert::YES>(
            key, this->extract_key(bucket_slots[i]))) {
            case detail::equal_result::AVAILABLE:
              return bucket_probing_results{detail::equal_result::AVAILABLE, i};
            case detail::equal_result::EQUAL: {
              if constexpr (allows_duplicates) {
                continue;
              } else {
                return bucket_probing_results{detail::equal_result::EQUAL, i};
              }
            }
            default: continue;
          }
        }
        // returns dummy index `-1` for UNEQUAL
        return bucket_probing_results{detail::equal_result::UNEQUAL, -1};
      }();

      if constexpr (not allows_duplicates) {
        // If the key is already in the container, return false
        if (group.any(state == detail::equal_result::EQUAL)) { return false; }
      }

      auto const group_contains_available = group.ballot(state == detail::equal_result::AVAILABLE);
      if (group_contains_available) {
        auto const src_lane = __ffs(group_contains_available) - 1;
        auto status         = insert_result::CONTINUE;
        if (group.thread_rank() == src_lane) {
          if constexpr (SupportsErase) {
            status = attempt_insert(this->get_slot_ptr(*probing_iter, intra_bucket_index),
                                    bucket_slots[intra_bucket_index],
                                    val);
          } else {
            status = attempt_insert(this->get_slot_ptr(*probing_iter, intra_bucket_index),
                                    this->empty_slot_sentinel(),
                                    val);
          }
        }

        switch (group.shfl(status, src_lane)) {
          case insert_result::SUCCESS: return true;
          case insert_result::DUPLICATE: {
            if constexpr (allows_duplicates) {
              [[fallthrough]];
            } else {
              return false;
            }
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
   * @brief Inserts the given element into the container.
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
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");
#if __CUDA_ARCH__ < 700
    // Spinning to ensure that the write to the value part took place requires
    // independent thread scheduling introduced with the Volta architecture.
    static_assert(
      cuco::detail::is_packable<value_type>(),
      "insert_and_find is not supported for pair types larger than 8 bytes on pre-Volta GPUs.");
#endif

    auto const val = this->heterogeneous_value(value);
    auto const key = this->extract_key(val);
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref_[*probing_iter];

      for (auto i = 0; i < bucket_size; ++i) {
        auto const eq_res = this->predicate_.template operator()<is_insert::YES>(
          key, this->extract_key(bucket_slots[i]));
        auto* slot_ptr = this->get_slot_ptr(*probing_iter, i);

        // If the key is already in the container, return false
        if (eq_res == detail::equal_result::EQUAL) {
          if constexpr (has_payload) {
            // wait to ensure that the write to the value part also took place
            this->wait_for_payload(slot_ptr->second, this->empty_value_sentinel());
          }
          return {iterator{slot_ptr}, false};
        }
        if (eq_res == detail::equal_result::AVAILABLE) {
          switch (this->attempt_insert_stable(slot_ptr, bucket_slots[i], val)) {
            case insert_result::SUCCESS: {
              if constexpr (has_payload) {
                // wait to ensure that the write to the value part also took place
                this->wait_for_payload(slot_ptr->second, this->empty_value_sentinel());
              }
              return {iterator{slot_ptr}, true};
            }
            case insert_result::DUPLICATE: {
              if constexpr (has_payload) {
                // wait to ensure that the write to the value part also took place
                this->wait_for_payload(slot_ptr->second, this->empty_value_sentinel());
              }
              return {iterator{slot_ptr}, false};
            }
            default: continue;
          }
        }
      }
      ++probing_iter;
      if (*probing_iter == init_idx) { return {this->end(), false}; }
    };
  }

  /**
   * @brief Inserts the given element into the container.
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
#if __CUDA_ARCH__ < 700
    // Spinning to ensure that the write to the value part took place requires
    // independent thread scheduling introduced with the Volta architecture.
    static_assert(
      cuco::detail::is_packable<value_type>(),
      "insert_and_find is not supported for pair types larger than 8 bytes on pre-Volta GPUs.");
#endif

    auto const val = this->heterogeneous_value(value);
    auto const key = this->extract_key(val);
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(group, key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref_[*probing_iter];

      auto const [state, intra_bucket_index] = [&]() {
        auto res = detail::equal_result::UNEQUAL;
        for (auto i = 0; i < bucket_size; ++i) {
          res = this->predicate_.template operator()<is_insert::YES>(
            key, this->extract_key(bucket_slots[i]));
          if (res != detail::equal_result::UNEQUAL) { return bucket_probing_results{res, i}; }
        }
        // returns dummy index `-1` for UNEQUAL
        return bucket_probing_results{res, -1};
      }();

      auto* slot_ptr = this->get_slot_ptr(*probing_iter, intra_bucket_index);

      // If the key is already in the container, return false
      auto const group_finds_equal = group.ballot(state == detail::equal_result::EQUAL);
      if (group_finds_equal) {
        auto const src_lane = __ffs(group_finds_equal) - 1;
        auto const res      = group.shfl(reinterpret_cast<intptr_t>(slot_ptr), src_lane);
        if (group.thread_rank() == src_lane) {
          if constexpr (has_payload) {
            // wait to ensure that the write to the value part also took place
            this->wait_for_payload(slot_ptr->second, this->empty_value_sentinel());
          }
        }
        group.sync();
        return {iterator{reinterpret_cast<value_type*>(res)}, false};
      }

      auto const group_contains_available = group.ballot(state == detail::equal_result::AVAILABLE);
      if (group_contains_available) {
        auto const src_lane = __ffs(group_contains_available) - 1;
        auto const res      = group.shfl(reinterpret_cast<intptr_t>(slot_ptr), src_lane);
        auto const status   = [&, target_idx = intra_bucket_index]() {
          if (group.thread_rank() != src_lane) { return insert_result::CONTINUE; }
          return this->attempt_insert_stable(slot_ptr, bucket_slots[target_idx], val);
        }();

        switch (group.shfl(status, src_lane)) {
          case insert_result::SUCCESS: {
            if (group.thread_rank() == src_lane) {
              if constexpr (has_payload) {
                // wait to ensure that the write to the value part also took place
                this->wait_for_payload(slot_ptr->second, this->empty_value_sentinel());
              }
            }
            group.sync();
            return {iterator{reinterpret_cast<value_type*>(res)}, true};
          }
          case insert_result::DUPLICATE: {
            if (group.thread_rank() == src_lane) {
              if constexpr (has_payload) {
                // wait to ensure that the write to the value part also took place
                this->wait_for_payload(slot_ptr->second, this->empty_value_sentinel());
              }
            }
            group.sync();
            return {iterator{reinterpret_cast<value_type*>(res)}, false};
          }
          default: continue;
        }
      } else {
        ++probing_iter;
        if (*probing_iter == init_idx) { return {this->end(), false}; }
      }
    }
  }

  /**
   * @brief Erases an element.
   *
   * @tparam ProbeKey Input type which is convertible to 'key_type'
   *
   * @param key The element to erase
   *
   * @return True if the given element is successfully erased
   */
  template <typename ProbeKey>
  __device__ bool erase(ProbeKey key) noexcept
  {
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");

    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref_[*probing_iter];

      for (auto& slot_content : bucket_slots) {
        auto const eq_res =
          this->predicate_.template operator()<is_insert::NO>(key, this->extract_key(slot_content));

        // Key doesn't exist, return false
        if (eq_res == detail::equal_result::EMPTY) { return false; }
        // Key exists, return true if successfully deleted
        if (eq_res == detail::equal_result::EQUAL) {
          auto const intra_bucket_index = cuda::std::distance(bucket_slots.begin(), &slot_content);
          switch (attempt_insert_stable(this->get_slot_ptr(*probing_iter, intra_bucket_index),
                                        slot_content,
                                        this->erased_slot_sentinel())) {
            case insert_result::SUCCESS: return true;
            case insert_result::DUPLICATE: return false;
            default: continue;
          }
        }
      }
      ++probing_iter;
      if (*probing_iter == init_idx) { return false; }
    }
  }

  /**
   * @brief Erases an element.
   *
   * @tparam ProbeKey Input type which is convertible to 'key_type'
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform group erase
   * @param key The element to erase
   *
   * @return True if the given element is successfully erased
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ bool erase(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                        ProbeKey key) noexcept
  {
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(group, key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref_[*probing_iter];

      auto const [state, intra_bucket_index] = [&]() {
        auto res = detail::equal_result::UNEQUAL;
        for (auto i = 0; i < bucket_size; ++i) {
          res = this->predicate_.template operator()<is_insert::NO>(
            key, this->extract_key(bucket_slots[i]));
          if (res != detail::equal_result::UNEQUAL) { return bucket_probing_results{res, i}; }
        }
        // returns dummy index `-1` for UNEQUAL
        return bucket_probing_results{res, -1};
      }();

      auto const group_contains_equal = group.ballot(state == detail::equal_result::EQUAL);
      if (group_contains_equal) {
        auto const src_lane = __ffs(group_contains_equal) - 1;
        auto const status =
          (group.thread_rank() == src_lane)
            ? attempt_insert_stable(this->get_slot_ptr(*probing_iter, intra_bucket_index),
                                    bucket_slots[intra_bucket_index],
                                    this->erased_slot_sentinel())
            : insert_result::CONTINUE;

        switch (group.shfl(status, src_lane)) {
          case insert_result::SUCCESS: return true;
          case insert_result::DUPLICATE: return false;
          default: continue;
        }
      }

      // Key doesn't exist, return false
      if (group.any(state == detail::equal_result::EMPTY)) { return false; }

      ++probing_iter;
      if (*probing_iter == init_idx) { return false; }
    }
  }

  /**
   * @brief Indicates whether the probe key `key` was inserted into the container.
   *
   * @note If the probe key `key` was inserted into the container, returns true. Otherwise, returns
   * false.
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
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      // TODO atomic_ref::load if insert operator is present
      auto const bucket_slots = storage_ref_[*probing_iter];

      for (auto i = 0; i < bucket_size; ++i) {
        switch (this->predicate_.template operator()<is_insert::NO>(
          key, this->extract_key(bucket_slots[i]))) {
          case detail::equal_result::UNEQUAL: continue;
          case detail::equal_result::EMPTY: return false;
          case detail::equal_result::EQUAL: return true;
        }
      }
      ++probing_iter;
      if (*probing_iter == init_idx) { return false; }
    }
  }

  /**
   * @brief Indicates whether the probe key `key` was inserted into the container.
   *
   * @note If the probe key `key` was inserted into the container, returns true. Otherwise, returns
   * false.
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
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(group, key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref_[*probing_iter];

      auto const state = [&]() {
        auto res = detail::equal_result::UNEQUAL;
        for (auto i = 0; i < bucket_size; ++i) {
          res = this->predicate_.template operator()<is_insert::NO>(
            key, this->extract_key(bucket_slots[i]));
          if (res != detail::equal_result::UNEQUAL) { return res; }
        }
        return res;
      }();

      if (group.any(state == detail::equal_result::EQUAL)) { return true; }
      if (group.any(state == detail::equal_result::EMPTY)) { return false; }

      ++probing_iter;
      if (*probing_iter == init_idx) { return false; }
    }
  }

  /**
   * @brief Finds an element in the container with key equivalent to the probe key.
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
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      // TODO atomic_ref::load if insert operator is present
      auto const bucket_slots = storage_ref_[*probing_iter];

      for (auto i = 0; i < bucket_size; ++i) {
        switch (this->predicate_.template operator()<is_insert::NO>(
          key, this->extract_key(bucket_slots[i]))) {
          case detail::equal_result::EMPTY: {
            return this->end();
          }
          case detail::equal_result::EQUAL: {
            return iterator{this->get_slot_ptr(*probing_iter, i)};
          }
          default: continue;
        }
      }
      ++probing_iter;
      if (*probing_iter == init_idx) { return this->end(); }
    }
  }

  /**
   * @brief Finds an element in the container with key equivalent to the probe key.
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
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(group, key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      auto const bucket_slots = storage_ref_[*probing_iter];

      auto const [state, intra_bucket_index] = [&]() {
        auto res = detail::equal_result::UNEQUAL;
        for (auto i = 0; i < bucket_size; ++i) {
          res = this->predicate_.template operator()<is_insert::NO>(
            key, this->extract_key(bucket_slots[i]));
          if (res != detail::equal_result::UNEQUAL) { return bucket_probing_results{res, i}; }
        }
        // returns dummy index `-1` for UNEQUAL
        return bucket_probing_results{res, -1};
      }();

      // Find a match for the probe key, thus return an iterator to the entry
      auto const group_finds_match = group.ballot(state == detail::equal_result::EQUAL);
      if (group_finds_match) {
        auto const src_lane = __ffs(group_finds_match) - 1;
        auto const res      = group.shfl(
          reinterpret_cast<intptr_t>(this->get_slot_ptr(*probing_iter, intra_bucket_index)),
          src_lane);
        return iterator{reinterpret_cast<value_type*>(res)};
      }

      // Find an empty slot, meaning that the probe key isn't present in the container
      if (group.any(state == detail::equal_result::EMPTY)) { return this->end(); }

      ++probing_iter;
      if (*probing_iter == init_idx) { return this->end(); }
    }
  }

  /**
   * @brief Counts the occurrence of a given key contained in the container
   *
   * @tparam ProbeKey Probe key type
   *
   * @param key The key to count for
   *
   * @return Number of occurrences found by the current thread
   */
  template <typename ProbeKey>
  [[nodiscard]] __device__ size_type count(ProbeKey key) const noexcept
  {
    if constexpr (not allows_duplicates) {
      return static_cast<size_type>(this->contains(key));
    } else {
      auto probing_iter =
        probing_scheme_.template make_iterator<bucket_size>(key, storage_ref_.extent());
      auto const init_idx = *probing_iter;
      size_type count     = 0;

      while (true) {
        auto const bucket_slots                = storage_ref_[*probing_iter];
        cuda::std::int32_t equals[bucket_size] = {0};
        bool empty_found                       = false;

#pragma unroll bucket_size
        for (cuda::std::int32_t i = 0; i < bucket_size; ++i) {
          auto const result =
            predicate_.template operator()<is_insert::NO>(key, this->extract_key(bucket_slots[i]));
          equals[i] = (result == detail::equal_result::EQUAL);
          if (result == detail::equal_result::EMPTY) {
            empty_found = true;
            break;
          }
        }

        count += thrust::reduce(thrust::seq, equals, equals + bucket_size);

        if (empty_found) { return count; }

        ++probing_iter;
        if (*probing_iter == init_idx) { return count; }
      }
    }
  }

  /**
   * @brief Counts the occurrence of a given key contained in the container
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
  [[nodiscard]] __device__ size_type
  count(cooperative_groups::thread_block_tile<cg_size, ParentCG> group, ProbeKey key) const noexcept
  {
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(group, key, storage_ref_.extent());
    auto const init_idx = *probing_iter;
    size_type count     = 0;

    while (true) {
      auto const bucket_slots                = storage_ref_[*probing_iter];
      cuda::std::int32_t equals[bucket_size] = {0};
      bool empty_found                       = false;

#pragma unroll bucket_size
      for (cuda::std::int32_t i = 0; i < bucket_size; ++i) {
        auto const result =
          predicate_.template operator()<is_insert::NO>(key, this->extract_key(bucket_slots[i]));
        equals[i] = (result == detail::equal_result::EQUAL);
        if (result == detail::equal_result::EMPTY) {
          empty_found = true;
          break;
        }
      }

      count += thrust::reduce(thrust::seq, equals, equals + bucket_size);

      if (group.any(empty_found)) { return count; }

      ++probing_iter;
      if (*probing_iter == init_idx) { return count; }
    }
  }

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
   * @tparam InputProbeIt Device accessible input iterator
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the `InputProbeIt`'s `value_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam AtomicCounter Integral atomic counter type that follows the same semantics as
   * `cuda::(std::)atomic(_ref)`
   *
   * @param block Thread block this operation is executed in
   * @param input_probe_begin Beginning of the input sequence of keys
   * @param input_probe_end End of the input sequence of keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param atomic_counter Atomic object of integral type that is used to count the
   * number of output elements
   */
  template <int BlockSize,
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
    auto constexpr is_outer = false;
    auto const n = cuco::detail::distance(input_probe_begin, input_probe_end);  // TODO include
    auto const always_true_stencil = thrust::constant_iterator<bool>(true);
    auto const identity_predicate  = cuda::std::identity{};
    this->retrieve_impl<is_outer, BlockSize>(block,
                                             input_probe_begin,
                                             n,
                                             always_true_stencil,
                                             identity_predicate,
                                             output_probe,
                                             output_match,
                                             atomic_counter);
  }

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
   * If a key `k` has no matches in the container, then `{key, empty_slot_sentinel}` will be added
   * to the output sequence.
   *
   * @tparam BlockSize Size of the thread block this operation is executed in
   * @tparam InputProbeIt Device accessible input iterator
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the `InputProbeIt`'s `value_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam AtomicCounter Integral atomic counter type that follows the same semantics as
   * `cuda::(std::)atomic(_ref)`
   *
   * @param block Thread block this operation is executed in
   * @param input_probe_begin Beginning of the input sequence of keys
   * @param input_probe_end End of the input sequence of keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param atomic_counter Atomic object of integral type that is used to count the
   * number of output elements
   */
  template <int BlockSize,
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
    auto constexpr is_outer = true;
    auto const n = cuco::detail::distance(input_probe_begin, input_probe_end);  // TODO include
    auto const always_true_stencil = thrust::constant_iterator<bool>(true);
    auto const identity_predicate  = cuda::std::identity{};
    this->retrieve_impl<is_outer, BlockSize>(block,
                                             input_probe_begin,
                                             n,
                                             always_true_stencil,
                                             identity_predicate,
                                             output_probe,
                                             output_match,
                                             atomic_counter);
  }

  /**
   * @brief Retrieves all the slots corresponding to all keys in the range `[input_probe_begin,
   * input_probe_end)` if `pred` of the corresponding stencil returns true.
   *
   * If key `k = *(first + i)` exists in the container and `pred( *(stencil + i) )` returns true,
   * copies `k` to `output_probe` and associated slot contents to `output_match`,
   * respectively. The output order is unspecified.
   *
   * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
   * Use `count()` to determine the size of the output range.
   *
   * @tparam BlockSize Size of the thread block this operation is executed in
   * @tparam InputProbeIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool`
   * and argument type is convertible from `std::iterator_traits<StencilIt>::value_type`
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the `InputProbeIt`'s `value_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam AtomicCounter Integral atomic counter type that follows the same semantics as
   * `cuda::(std::)atomic(_ref)`
   *
   * @param block Thread block this operation is executed in
   * @param input_probe_begin Beginning of the input sequence of keys
   * @param input_probe_end End of the input sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil + n)`
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param atomic_counter Atomic object of integral type that is used to count the
   * number of output elements
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
    auto constexpr is_outer = false;
    auto const n            = cuco::detail::distance(input_probe_begin, input_probe_end);
    this->retrieve_impl<is_outer, BlockSize>(
      block, input_probe_begin, n, stencil, pred, output_probe, output_match, atomic_counter);
  }

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
   * If `IsOuter == true` and a key `k` has no matches in the container, then `{key,
   * empty_slot_sentinel}` will be added to the output sequence.
   *
   * @tparam IsOuter Flag indicating if an inner or outer retrieve operation should be performed
   * @tparam BlockSize Size of the thread block this operation is executed in
   * @tparam InputProbeIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool`
   * and argument type is convertible from `std::iterator_traits<StencilIt>::value_type`
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the `InputProbeIt`'s `value_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam AtomicCounter Integral atomic type that follows the same semantics as
   * `cuda::(std::)atomic(_ref)`
   *
   * @param block Thread block this operation is executed in
   * @param input_probe Beginning of the input sequence of keys
   * @param n Number of input keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil + n)`
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param atomic_counter Atomic object of integral type that is used to count the
   * number of output elements
   */
  template <bool IsOuter,
            int BlockSize,
            class InputProbeIt,
            class StencilIt,
            class Predicate,
            class OutputProbeIt,
            class OutputMatchIt,
            class AtomicCounter>
  __device__ void retrieve_impl(cooperative_groups::thread_block const& block,
                                InputProbeIt input_probe,
                                cuco::detail::index_type n,
                                StencilIt stencil,
                                Predicate pred,
                                OutputProbeIt output_probe,
                                OutputMatchIt output_match,
                                AtomicCounter& atomic_counter) const
  {
    namespace cg = cooperative_groups;

    if (n == 0) { return; }

    using probe_type = typename cuda::std::iterator_traits<InputProbeIt>::value_type;

    // tuning parameter
    auto constexpr buffer_multiplier = 1;
    static_assert(buffer_multiplier > 0);

    auto constexpr probing_tile_size  = cg_size;
    auto constexpr flushing_tile_size = cuco::detail::warp_size();
    static_assert(flushing_tile_size >= probing_tile_size);

    auto constexpr num_flushing_tiles   = BlockSize / flushing_tile_size;
    auto constexpr max_matches_per_step = flushing_tile_size * bucket_size;
    auto constexpr buffer_size = buffer_multiplier * max_matches_per_step + flushing_tile_size;

    auto const flushing_tile = cg::tiled_partition<flushing_tile_size, cg::thread_block>(block);
    auto const probing_tile  = cg::tiled_partition<probing_tile_size, cg::thread_block>(block);

    auto const flushing_tile_id = flushing_tile.meta_group_rank();
    auto const stride           = probing_tile.meta_group_size();
    auto idx                    = probing_tile.meta_group_rank();

    __shared__ cuco::pair<probe_type, value_type> buffers[num_flushing_tiles][buffer_size];
    __shared__ cuda::std::int32_t counters[num_flushing_tiles];

    if (flushing_tile.thread_rank() == 0) { counters[flushing_tile_id] = 0; }
    flushing_tile.sync();

    auto flush_buffers = [&](auto tile) {
      size_type offset = 0;
      auto const count = counters[flushing_tile_id];
      auto const rank  = tile.thread_rank();
      if (rank == 0) { offset = atomic_counter.fetch_add(count, cuda::memory_order_relaxed); }
      offset = tile.shfl(offset, 0);

      // flush_buffers
      for (auto i = rank; i < count; i += tile.size()) {
        *(output_probe + offset + i) = buffers[flushing_tile_id][i].first;
        *(output_match + offset + i) = buffers[flushing_tile_id][i].second;
      }
    };

    while (flushing_tile.any(idx < n)) {
      bool active_flag = idx < n and pred(*(stencil + idx));
      auto const active_flushing_tile =
        cg::binary_partition<flushing_tile_size>(flushing_tile, active_flag);

      if (active_flag) {
        // perform probing
        // make sure the flushing_tile is converged at this point to get a coalesced load
        auto const probe_key = *(input_probe + idx);

        auto probing_iter = probing_scheme_.template make_iterator<bucket_size>(
          probing_tile, probe_key, storage_ref_.extent());
        auto const init_idx = *probing_iter;

        bool running                      = true;
        [[maybe_unused]] bool found_match = false;

        bool equals[bucket_size];
        cuda::std::uint32_t exists[bucket_size];

        while (active_flushing_tile.any(running)) {
          if (running) {
            // TODO atomic_ref::load if insert operator is present
            auto const bucket_slots = this->storage_ref_[*probing_iter];

#pragma unroll bucket_size
            for (cuda::std::int32_t i = 0; i < bucket_size; ++i) {
              equals[i] = false;
              if (running) {
                // inspect slot content
                switch (this->predicate_.template operator()<is_insert::NO>(
                  probe_key, this->extract_key(bucket_slots[i]))) {
                  case detail::equal_result::EMPTY: {
                    running = false;
                    break;
                  }
                  case detail::equal_result::EQUAL: {
                    if constexpr (!AllowsDuplicates) { running = false; }
                    equals[i] = true;
                    break;
                  }
                  default: {
                    break;
                  }
                }
              }
            }

            probing_tile.sync();
            running = probing_tile.all(running);
#pragma unroll bucket_size
            for (cuda::std::int32_t i = 0; i < bucket_size; ++i) {
              exists[i] = probing_tile.ballot(equals[i]);
            }

            // Fill the buffer if any matching keys are found
            auto const lane_id = probing_tile.thread_rank();
            if (thrust::any_of(thrust::seq, exists, exists + bucket_size, cuda::std::identity{})) {
              if constexpr (IsOuter) { found_match = true; }

              cuda::std::int32_t num_matches[bucket_size];

              for (cuda::std::int32_t i = 0; i < bucket_size; ++i) {
                num_matches[i] = __popc(exists[i]);
              }

              cuda::std::int32_t output_idx;
              if (lane_id == 0) {
                auto const total_matches =
                  thrust::reduce(thrust::seq, num_matches, num_matches + bucket_size);
                auto ref = cuda::atomic_ref<cuda::std::int32_t, cuda::thread_scope_block>{
                  counters[flushing_tile_id]};
                output_idx = ref.fetch_add(total_matches, cuda::memory_order_relaxed);
              }
              output_idx = probing_tile.shfl(output_idx, 0);

              cuda::std::int32_t matches_offset = 0;
#pragma unroll bucket_size
              for (cuda::std::int32_t i = 0; i < bucket_size; ++i) {
                if (equals[i]) {
                  auto const lane_offset = detail::count_least_significant_bits(exists[i], lane_id);
                  buffers[flushing_tile_id][output_idx + matches_offset + lane_offset] = {
                    probe_key, bucket_slots[i]};
                }
                matches_offset += num_matches[i];
              }
            }
            // Special handling for outer cases where no match is found
            if constexpr (IsOuter) {
              if (!running) {
                if (!found_match and lane_id == 0) {
                  auto ref = cuda::atomic_ref<cuda::std::int32_t, cuda::thread_scope_block>{
                    counters[flushing_tile_id]};
                  auto const output_idx = ref.fetch_add(1, cuda::memory_order_relaxed);
                  buffers[flushing_tile_id][output_idx] = {probe_key, this->empty_slot_sentinel()};
                }
              }
            }
          }  // if running

          active_flushing_tile.sync();
          // if the buffer has not enough empty slots for the next iteration
          if (counters[flushing_tile_id] > (buffer_size - max_matches_per_step)) {
            flush_buffers(active_flushing_tile);
            active_flushing_tile.sync();

            // reset buffer counter
            if (active_flushing_tile.thread_rank() == 0) { counters[flushing_tile_id] = 0; }
            active_flushing_tile.sync();
          }

          // onto the next probing bucket
          ++probing_iter;
          if (*probing_iter == init_idx) { running = false; }
        }  // while running
      }  // if active_flag

      // onto the next key
      idx += stride;
    }

    flushing_tile.sync();
    // entire flusing_tile has finished; flush remaining elements
    if (counters[flushing_tile_id] > 0) { flush_buffers(flushing_tile); }
  }

  /**
   * @brief For a given key, applies the function object `callback_op` to the copy of all
   * corresponding matches found in the container.
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @tparam ProbeKey Probe key type
   * @tparam CallbackOp Type of unary callback function object
   *
   * @param key The key to search for
   * @param callback_op Function to apply to every matched slot
   */
  template <class ProbeKey, class CallbackOp>
  __device__ void for_each(ProbeKey key, CallbackOp&& callback_op) const noexcept
  {
    static_assert(cg_size == 1, "Non-CG operation is incompatible with the current probing scheme");
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(key, storage_ref_.extent());
    auto const init_idx = *probing_iter;

    while (true) {
      // TODO atomic_ref::load if insert operator is present
      auto const bucket_slots = this->storage_ref_[*probing_iter];

      for (cuda::std::int32_t i = 0; i < bucket_size; ++i) {
        switch (this->predicate_.template operator()<is_insert::NO>(
          key, this->extract_key(bucket_slots[i]))) {
          case detail::equal_result::EMPTY: {
            return;
          }
          case detail::equal_result::EQUAL: {
            callback_op(bucket_slots[i]);
            continue;
          }
          default: continue;
        }
      }
      ++probing_iter;
      if (*probing_iter == init_idx) { return; }
    }
  }

  /**
   * @brief For a given key, applies the function object `callback_op` to the copy of all
   * corresponding matches found in the container.
   *
   * @note This function uses cooperative group semantics, meaning that any thread may call the
   * callback if it finds a matching element. If multiple elements are found within the same group,
   * each thread with a match will call the callback with its associated element.
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
   * @param callback_op Function to apply to every matched slot
   */
  template <class ProbeKey, class CallbackOp, typename ParentCG>
  __device__ void for_each(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                           ProbeKey key,
                           CallbackOp&& callback_op) const noexcept
  {
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(group, key, storage_ref_.extent());
    auto const init_idx = *probing_iter;
    bool empty          = false;

    while (true) {
      // TODO atomic_ref::load if insert operator is present
      auto const bucket_slots = this->storage_ref_[*probing_iter];

      for (cuda::std::int32_t i = 0; i < bucket_size and !empty; ++i) {
        switch (this->predicate_.template operator()<is_insert::NO>(
          key, this->extract_key(bucket_slots[i]))) {
          case detail::equal_result::EMPTY: {
            empty = true;
            continue;
          }
          case detail::equal_result::EQUAL: {
            callback_op(bucket_slots[i]);
            continue;
          }
          default: {
            continue;
          }
        }
      }
      if (group.any(empty)) { return; }

      ++probing_iter;
      if (*probing_iter == init_idx) { return; }
    }
  }

  /**
   * @brief Applies the function object `callback_op` to the copy of every slot in the container
   * with key equivalent to the probe key and can additionally perform work that requires
   * synchronizing the Cooperative Group performing this operation.
   *
   * @note This function uses cooperative group semantics, meaning that any thread may call the
   * callback if it finds a matching element. If multiple elements are found within the same group,
   * each thread with a match will call the callback with its associated element.
   *
   * @note Synchronizing `group` within `callback_op` is undefined behavior.
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @note The `sync_op` function can be used to perform work that requires synchronizing threads in
   * `group` inbetween probing steps, where the number of probing steps performed between
   * synchronization points is capped by `bucket_size * cg_size`. The functor will be called right
   * after the current probing bucket has been traversed.
   *
   * @tparam ProbeKey Probe key type
   * @tparam CallbackOp Type of unary callback function object
   * @tparam SyncOp Type of function object which accepts the current `group` object
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param group The Cooperative Group used to perform this operation
   * @param key The key to search for
   * @param callback_op Function to apply to every matched slot
   * @param sync_op Function that is allowed to synchronize `group` inbetween probing buckets
   */
  template <class ProbeKey, class CallbackOp, class SyncOp, typename ParentCG>
  __device__ void for_each(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
                           ProbeKey key,
                           CallbackOp&& callback_op,
                           SyncOp&& sync_op) const noexcept
  {
    auto probing_iter =
      probing_scheme_.template make_iterator<bucket_size>(group, key, storage_ref_.extent());
    auto const init_idx = *probing_iter;
    bool empty          = false;

    while (true) {
      // TODO atomic_ref::load if insert operator is present
      auto const bucket_slots = this->storage_ref_[*probing_iter];

      for (cuda::std::int32_t i = 0; i < bucket_size and !empty; ++i) {
        switch (this->predicate_.template operator()<is_insert::NO>(
          key, this->extract_key(bucket_slots[i]))) {
          case detail::equal_result::EMPTY: {
            empty = true;
            continue;
          }
          case detail::equal_result::EQUAL: {
            callback_op(bucket_slots[i]);
            continue;
          }
          default: {
            continue;
          }
        }
      }
      sync_op(group);
      if (group.any(empty)) { return; }

      ++probing_iter;
      if (*probing_iter == init_idx) { return; }
    }
  }

  /**
   * @brief Gets a pointer to the slot at the given probing index and intra-bucket index.
   *
   * @param probing_idx The current probing index
   * @param intra_bucket_idx The index within the bucket (0 for flat storage)
   * @return Pointer to the slot
   */
  __device__ value_type* get_slot_ptr(size_type probing_idx,
                                      cuda::std::int32_t intra_bucket_idx) const noexcept
  {
    return storage_ref_.data() + probing_idx + intra_bucket_idx;
  }

  /**
   * @brief Extracts the key from a given value type.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param value The input value
   *
   * @return The key
   */
  template <typename Value>
  [[nodiscard]] __host__ __device__ constexpr auto extract_key(Value value) const noexcept
  {
    if constexpr (has_payload) {
      return thrust::raw_reference_cast(value).first;
    } else {
      return thrust::raw_reference_cast(value);
    }
  }

  /**
   * @brief Extracts the payload from a given value type.
   *
   * @note This function is only available if `this->has_payload == true`
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param value The input value
   *
   * @return The payload
   */
  template <typename Value, typename Enable = cuda::std::enable_if_t<has_payload and sizeof(Value)>>
  [[nodiscard]] __host__ __device__ constexpr auto extract_payload(Value value) const noexcept
  {
    return thrust::raw_reference_cast(value).second;
  }

  /**
   * @brief Converts the given type to the container's native `value_type`.
   *
   * @tparam T Input type which is convertible to 'value_type'
   *
   * @param value The input value
   *
   * @return The converted object
   */
  template <typename T>
  [[nodiscard]] __device__ constexpr value_type native_value(T value) const noexcept
  {
    if constexpr (has_payload) {
      return {static_cast<key_type>(this->extract_key(value)), this->extract_payload(value)};
    } else {
      return static_cast<value_type>(value);
    }
  }

  /**
   * @brief Converts the given type to the container's native `value_type` while maintaining the
   * heterogeneous key type.
   *
   * @tparam T Input type which is convertible to 'value_type'
   *
   * @param value The input value
   *
   * @return The converted object
   */
  template <typename T>
  [[nodiscard]] __device__ constexpr auto heterogeneous_value(T value) const noexcept
  {
    if constexpr (has_payload and not cuda::std::is_same_v<T, value_type>) {
      using mapped_type = decltype(this->empty_value_sentinel());
      if constexpr (cuco::detail::is_cuda_std_pair_like<T>::value) {
        return cuco::pair{cuda::std::get<0>(value),
                          static_cast<mapped_type>(cuda::std::get<1>(value))};
      } else {
        // hail mary (convert using .first/.second members)
        return cuco::pair{thrust::raw_reference_cast(value.first),
                          static_cast<mapped_type>(value.second)};
      }
    } else {
      return thrust::raw_reference_cast(value);
    }
  }

  /**
   * @brief Gets the sentinel used to represent an erased slot.
   *
   * @return The sentinel value used to represent an erased slot
   */
  [[nodiscard]] __device__ constexpr value_type erased_slot_sentinel() const noexcept
  {
    if constexpr (has_payload) {
      return cuco::pair{this->erased_key_sentinel(), this->empty_value_sentinel()};
    } else {
      return this->erased_key_sentinel();
    }
  }

  /**
   * @brief Inserts the specified element with one single CAS operation.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param address Pointer to the slot in memory
   * @param expected Element to compare against
   * @param desired Element to insert
   *
   * @return Result of this operation, i.e., success/continue/duplicate
   */
  template <typename Value>
  [[nodiscard]] __device__ constexpr insert_result packed_cas(value_type* address,
                                                              value_type expected,
                                                              Value desired) noexcept
  {
    using packed_type =
      cuda::std::conditional_t<sizeof(value_type) == 4, cuda::std::uint32_t, cuda::std::uint64_t>;

    auto* slot_ptr     = reinterpret_cast<packed_type*>(address);
    auto* expected_ptr = reinterpret_cast<packed_type*>(&expected);
    auto* desired_ptr  = reinterpret_cast<packed_type*>(&desired);

    auto slot_ref = cuda::atomic_ref<packed_type, Scope>{*slot_ptr};

    auto const success =
      slot_ref.compare_exchange_strong(*expected_ptr, *desired_ptr, cuda::memory_order_relaxed);

    if (success) {
      return insert_result::SUCCESS;
    } else {
      return this->predicate_.equal_to(this->extract_key(desired), this->extract_key(expected)) ==
                 detail::equal_result::EQUAL
               ? insert_result::DUPLICATE
               : insert_result::CONTINUE;
    }
  }

  /**
   * @brief Inserts the specified element with two back-to-back CAS operations.
   *
   * @note This CAS can be used exclusively for `cuco::op::insert` operations.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param address Pointer to the slot in memory
   * @param expected Element to compare against
   * @param desired Element to insert
   *
   * @return Result of this operation, i.e., success/continue/duplicate
   */
  template <typename Value>
  [[nodiscard]] __device__ constexpr insert_result back_to_back_cas(value_type* address,
                                                                    value_type expected,
                                                                    Value desired) noexcept
  {
    using mapped_type = cuda::std::decay_t<decltype(this->empty_value_sentinel())>;

    auto expected_key     = expected.first;
    auto expected_payload = this->empty_value_sentinel();

    cuda::atomic_ref<key_type, Scope> key_ref(address->first);
    cuda::atomic_ref<mapped_type, Scope> payload_ref(address->second);

    auto const key_cas_success = key_ref.compare_exchange_strong(
      expected_key, static_cast<key_type>(desired.first), cuda::memory_order_relaxed);
    auto payload_cas_success = payload_ref.compare_exchange_strong(
      expected_payload, desired.second, cuda::memory_order_relaxed);

    // if key success
    if (key_cas_success) {
      while (not payload_cas_success) {
        payload_cas_success =
          payload_ref.compare_exchange_strong(expected_payload = this->empty_value_sentinel(),
                                              desired.second,
                                              cuda::memory_order_relaxed);
      }
      return insert_result::SUCCESS;
    } else if (payload_cas_success) {
      // This is insert-specific, cannot for `erase` operations
      payload_ref.store(this->empty_value_sentinel(), cuda::memory_order_relaxed);
    }

    // Our key was already present in the slot, so our key is a duplicate
    // Shouldn't use `predicate` operator directly since it includes a redundant bitwise compare
    if (this->predicate_.equal_to(desired.first, expected_key) == detail::equal_result::EQUAL) {
      return insert_result::DUPLICATE;
    }

    return insert_result::CONTINUE;
  }

  /**
   * @brief Inserts the specified element with CAS-dependent write operations.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param address Pointer to the slot in memory
   * @param expected Element to compare against
   * @param desired Element to insert
   *
   * @return Result of this operation, i.e., success/continue/duplicate
   */
  template <typename Value>
  [[nodiscard]] __device__ constexpr insert_result cas_dependent_write(value_type* address,
                                                                       value_type expected,
                                                                       Value desired) noexcept
  {
    using mapped_type = cuda::std::decay_t<decltype(this->empty_value_sentinel())>;

    cuda::atomic_ref<key_type, Scope> key_ref(address->first);
    auto expected_key  = expected.first;
    auto const success = key_ref.compare_exchange_strong(
      expected_key, static_cast<key_type>(desired.first), cuda::memory_order_relaxed);

    // if key success
    if (success) {
      cuda::atomic_ref<mapped_type, Scope> payload_ref(address->second);
      payload_ref.store(desired.second, cuda::memory_order_relaxed);
      return insert_result::SUCCESS;
    }

    // Our key was already present in the slot, so our key is a duplicate
    // Shouldn't use `predicate` operator directly since it includes a redundant bitwise compare
    if (this->predicate_.equal_to(desired.first, expected_key) == detail::equal_result::EQUAL) {
      return insert_result::DUPLICATE;
    }

    return insert_result::CONTINUE;
  }

  /**
   * @brief Attempts to insert an element into a slot.
   *
   * @note Dispatches the correct implementation depending on the container
   * type and presence of other operator mixins.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param address Pointer to the slot in memory
   * @param expected Element to compare against
   * @param desired Element to insert
   *
   * @return Result of this operation, i.e., success/continue/duplicate
   */
  template <typename Value>
  [[nodiscard]] __device__ insert_result attempt_insert(value_type* address,
                                                        value_type expected,
                                                        Value desired) noexcept
  {
    if constexpr (sizeof(value_type) <= 8) {
      return packed_cas(address, expected, desired);
    } else {
#if (__CUDA_ARCH__ < 700)
      return cas_dependent_write(address, expected, desired);
#else
      return back_to_back_cas(address, expected, desired);
#endif
    }
  }

  /**
   * @brief Attempts to insert an element into a slot.
   *
   * @note Dispatches the correct implementation depending on the container
   * type and presence of other operator mixins.
   *
   * @note `stable` indicates that the payload will only be updated once from the sentinel value to
   * the desired value, meaning there can be no ABA situations.
   *
   * @tparam Value Input type which is convertible to 'value_type'
   *
   * @param address Pointer to the slot in memory
   * @param expected Element to compare against
   * @param desired Element to insert
   *
   * @return Result of this operation, i.e., success/continue/duplicate
   */
  template <typename Value>
  [[nodiscard]] __device__ insert_result attempt_insert_stable(value_type* address,
                                                               value_type expected,
                                                               Value desired) noexcept
  {
    if constexpr (sizeof(value_type) <= 8) {
      return packed_cas(address, expected, desired);
    } else {
      return cas_dependent_write(address, expected, desired);
    }
  }

  /**
   * @brief Waits until the slot payload has been updated
   *
   * @note The function will return once the slot payload is no longer equal to the sentinel
   * value.
   *
   * @tparam T Map slot type
   *
   * @param slot The target slot to check payload with
   * @param sentinel The slot sentinel value
   */
  template <typename T>
  __device__ void wait_for_payload(T& slot, T sentinel) const noexcept
  {
    auto ref = cuda::atomic_ref<T, Scope>{slot};
    T current;
    // TODO exponential backoff strategy
    do {
      current = ref.load(cuda::std::memory_order_relaxed);
    } while (cuco::detail::bitwise_compare(current, sentinel));
  }

  // TODO: Clean up the sentinel handling since it's duplicated in ref and equal wrapper
  value_type empty_slot_sentinel_;  ///< Sentinel value indicating an empty slot
  detail::equal_wrapper<key_type, key_equal> predicate_;  ///< Key equality binary callable
  probing_scheme_type probing_scheme_;                    ///< Probing scheme
  storage_ref_type storage_ref_;                          ///< Slot storage ref
};

}  // namespace detail
}  // namespace cuco
