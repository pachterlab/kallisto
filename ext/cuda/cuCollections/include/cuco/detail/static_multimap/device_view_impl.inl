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

#include <cuco/detail/bitwise_compare.cuh>
#include <cuco/detail/static_multimap/kernels.cuh>
#include <cuco/detail/utils.cuh>

#include <thrust/type_traits/is_contiguous_iterator.h>

#include <cooperative_groups.h>
#include <cooperative_groups/memcpy_async.h>

namespace cuco {
template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
class static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view_impl_base {
 protected:
  // Import member type definitions from `static_multimap`
  using value_type          = value_type;
  using key_type            = Key;
  using mapped_type         = Value;
  using iterator            = pair_atomic_type*;
  using const_iterator      = pair_atomic_type const*;
  using probe_sequence_type = probe_sequence_type;

  /**
   * @brief Indicates if vector-load is used.
   *
   * Users have no explicit control on whether vector-load is used.
   *
   * @return Boolean indicating if vector-load is used.
   */
  static constexpr bool uses_vector_load() noexcept
  {
    return probe_sequence_type::uses_vector_load();
  }

  /**
   * @brief Returns the number of pairs loaded with each vector-load
   */
  static constexpr uint32_t vector_width() noexcept { return probe_sequence_type::vector_width(); }

  __host__ __device__ device_view_impl_base(pair_atomic_type* slots,
                                            std::size_t capacity,
                                            Key empty_key_sentinel,
                                            Value empty_value_sentinel) noexcept
    : probe_sequence_{slots, capacity},
      empty_key_sentinel_{empty_key_sentinel},
      empty_value_sentinel_{empty_value_sentinel}
  {
  }

  /**
   * @brief Returns the initial slot for a given key `k`
   *
   * To be used for Cooperative Group based probing.
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g the Cooperative Group for which the initial slot is needed
   * @param k The key to get the slot for
   * @return Pointer to the initial slot for `k`
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ __forceinline__ iterator
  initial_slot(cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
               ProbeKey const& k) noexcept
  {
    return probe_sequence_.initial_slot(g, k);
  }

  /**
   * @brief Returns the initial slot for a given key `k`
   *
   * To be used for Cooperative Group based probing.
   *
   * @tparam ProbeKey Probe key type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g the Cooperative Group for which the initial slot is needed
   * @param k The key to get the slot for
   * @return Pointer to the initial slot for `k`
   */
  template <typename ProbeKey, typename ParentCG>
  __device__ __forceinline__ const_iterator
  initial_slot(cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
               ProbeKey const& k) const noexcept
  {
    return probe_sequence_.initial_slot(g, k);
  }

  /**
   * @brief Given a slot `s`, returns the next slot.
   *
   * If `s` is the last slot, wraps back around to the first slot. To
   * be used for Cooperative Group based probing.
   *
   * @param s The slot to advance
   * @return The next slot after `s`
   */
  __device__ __forceinline__ iterator next_slot(iterator s) noexcept
  {
    return probe_sequence_.next_slot(s);
  }

  /**
   * @brief Given a slot `s`, returns the next slot.
   *
   * If `s` is the last slot, wraps back around to the first slot. To
   * be used for Cooperative Group based probing.
   *
   * @param s The slot to advance
   * @return The next slot after `s`
   */
  __device__ __forceinline__ const_iterator next_slot(const_iterator s) const noexcept
  {
    return probe_sequence_.next_slot(s);
  }

  /**
   * @brief Load two key/value pairs from the given slot to the target pair array.
   *
   * @param arr The pair array to be loaded
   * @param current_slot The given slot to load from
   */
  __device__ __forceinline__ void load_pair_array(value_type* arr,
                                                  const_iterator current_slot) const noexcept
  {
    if constexpr (sizeof(value_type) == 4) {
      auto const tmp = *reinterpret_cast<ushort4 const*>(current_slot);
      memcpy(&arr[0], &tmp, 2 * sizeof(value_type));
    } else {
      auto const tmp = *reinterpret_cast<uint4 const*>(current_slot);
      memcpy(&arr[0], &tmp, 2 * sizeof(value_type));
    }
  }

 public:
  /**
   * @brief Gets the sentinel value used to represent an empty key slot.
   *
   * @return The sentinel value used to represent an empty key slot
   */
  __host__ __device__ __forceinline__ Key get_empty_key_sentinel() const noexcept
  {
    return empty_key_sentinel_;
  }

  /**
   * @brief Gets the sentinel value used to represent an empty value slot.
   *
   * @return The sentinel value used to represent an empty value slot
   */
  __host__ __device__ __forceinline__ Value get_empty_value_sentinel() const noexcept
  {
    return empty_value_sentinel_;
  }

  /**
   * @brief Gets slots array.
   *
   * @return Slots array
   */
  __device__ __forceinline__ pair_atomic_type* get_slots() noexcept
  {
    return probe_sequence_.get_slots();
  }

  /**
   * @brief Gets slots array.
   *
   * @return Slots array
   */
  __device__ __forceinline__ pair_atomic_type const* get_slots() const noexcept
  {
    return probe_sequence_.get_slots();
  }

  /**
   * @brief Gets the maximum number of elements the hash map can hold.
   *
   * @return The maximum number of elements the hash map can hold
   */
  __host__ __device__ __forceinline__ std::size_t get_capacity() const noexcept
  {
    return probe_sequence_.get_capacity();
  }

 private:
  probe_sequence_type probe_sequence_;  ///< Probe sequence used to probe the hash map
  Key empty_key_sentinel_{};            ///< Key value that represents an empty slot
  Value empty_value_sentinel_{};        ///< Initial Value of empty slot
};  // class device_view_impl_base

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
class static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_mutable_view_impl
  : public device_view_impl_base {
 public:
  using value_type     = typename device_view_impl_base::value_type;
  using key_type       = typename device_view_impl_base::key_type;
  using mapped_type    = typename device_view_impl_base::mapped_type;
  using iterator       = typename device_view_impl_base::iterator;
  using const_iterator = typename device_view_impl_base::const_iterator;

 private:
  /**
   * @brief Enumeration of the possible results of attempting to insert into a hash slot.
   */
  enum class insert_result {
    CONTINUE,  ///< Insert did not succeed, continue trying to insert
    SUCCESS,   ///< New pair inserted successfully
    DUPLICATE  ///< Insert did not succeed, key is already present
  };

  /**
   * @brief Inserts the specified key/value pair with one single CAS operation.
   *
   * @param current_slot The slot to insert
   * @param insert_pair The pair to insert
   * @param key_equal The binary callable used to compare two keys for
   * equality
   * @return An insert result from the `insert_resullt` enumeration.
   */
  __device__ __forceinline__ insert_result packed_cas(iterator current_slot,
                                                      value_type const& insert_pair) noexcept
  {
    auto expected_key   = this->get_empty_key_sentinel();
    auto expected_value = this->get_empty_value_sentinel();

    cuco::detail::pair_converter<value_type> expected_pair{
      cuco::make_pair(expected_key, expected_value)};
    cuco::detail::pair_converter<value_type> new_pair{insert_pair};

    auto slot = reinterpret_cast<
      cuda::atomic<typename cuco::detail::pair_converter<value_type>::packed_type, Scope>*>(
      current_slot);

    bool success = slot->compare_exchange_strong(
      expected_pair.packed, new_pair.packed, cuda::std::memory_order_relaxed);
    if (success) { return insert_result::SUCCESS; }

    return insert_result::CONTINUE;
  }

  /**
   * @brief Inserts the specified key/value pair with two back-to-back CAS operations.
   *
   * @param current_slot The slot to insert
   * @param insert_pair The pair to insert
   * @return An insert result from the `insert_resullt` enumeration.
   */
  __device__ __forceinline__ insert_result back_to_back_cas(iterator current_slot,
                                                            value_type const& insert_pair) noexcept
  {
    using cuda::std::memory_order_relaxed;

    auto expected_key   = this->get_empty_key_sentinel();
    auto expected_value = this->get_empty_value_sentinel();

    // Back-to-back CAS for 8B/8B key/value pairs
    auto& slot_key   = current_slot->first;
    auto& slot_value = current_slot->second;

    bool key_success =
      slot_key.compare_exchange_strong(expected_key, insert_pair.first, memory_order_relaxed);
    bool value_success =
      slot_value.compare_exchange_strong(expected_value, insert_pair.second, memory_order_relaxed);

    if (key_success) {
      while (not value_success) {
        value_success =
          slot_value.compare_exchange_strong(expected_value = this->get_empty_value_sentinel(),
                                             insert_pair.second,
                                             memory_order_relaxed);
      }
      return insert_result::SUCCESS;
    } else if (value_success) {
      slot_value.store(this->get_empty_value_sentinel(), memory_order_relaxed);
    }

    return insert_result::CONTINUE;
  }

  /**
   * @brief Inserts the specified key/value pair with a CAS of the key and a dependent write
   * of the value.
   *
   * @param current_slot The slot to insert
   * @param insert_pair The pair to insert
   * @return An insert result from the `insert_resullt` enumeration.
   */
  __device__ __forceinline__ insert_result
  cas_dependent_write(iterator current_slot, value_type const& insert_pair) noexcept
  {
    using cuda::std::memory_order_relaxed;
    auto expected_key = this->get_empty_key_sentinel();

    auto& slot_key = current_slot->first;

    auto const key_success =
      slot_key.compare_exchange_strong(expected_key, insert_pair.first, memory_order_relaxed);

    if (key_success) {
      auto& slot_value = current_slot->second;
      slot_value.store(insert_pair.second, memory_order_relaxed);
      return insert_result::SUCCESS;
    }

    return insert_result::CONTINUE;
  }

 public:
  __host__ __device__ device_mutable_view_impl(pair_atomic_type* slots,
                                               std::size_t capacity,
                                               Key empty_key_sentinel,
                                               Value empty_value_sentinel) noexcept
    : device_view_impl_base{slots, capacity, empty_key_sentinel, empty_value_sentinel}
  {
  }

  /**
   * @brief Inserts the specified key/value pair into the map using vector loads.
   *
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam CG Cooperative Group type
   *
   * @param g The Cooperative Group that performs the insert
   * @param insert_pair The pair to insert
   * @return void.
   */
  template <bool uses_vector_load, typename CG>
  __device__ __forceinline__ cuda::std::enable_if_t<uses_vector_load, void> insert(
    CG g, value_type const& insert_pair) noexcept
  {
    auto current_slot = initial_slot(g, insert_pair.first);
    while (true) {
      value_type arr[2];
      load_pair_array(&arr[0], current_slot);

      // The user provide `key_equal` can never be used to compare against `empty_key_sentinel` as
      // the sentinel is not a valid key value. Therefore, first check for the sentinel
      auto const first_slot_is_empty =
        (detail::bitwise_compare(arr[0].first, this->get_empty_key_sentinel()));
      auto const second_slot_is_empty =
        (detail::bitwise_compare(arr[1].first, this->get_empty_key_sentinel()));
      auto const bucket_contains_empty = g.ballot(first_slot_is_empty or second_slot_is_empty);

      if (bucket_contains_empty) {
        // the first lane in the group with an empty slot will attempt the insert
        insert_result status{insert_result::CONTINUE};
        uint32_t src_lane = __ffs(bucket_contains_empty) - 1;
        if (g.thread_rank() == src_lane) {
          auto insert_location = first_slot_is_empty ? current_slot : current_slot + 1;
          // One single CAS operation since vector loads are dedicated to packable pairs
          status = packed_cas(insert_location, insert_pair);
        }

        // successful insert
        if (g.any(status == insert_result::SUCCESS)) { return; }
        // if we've gotten this far, a different key took our spot
        // before we could insert. We need to retry the insert on the
        // same bucket
      }
      // if there are no empty slots in the current bucket,
      // we move onto the next bucket
      else {
        current_slot = next_slot(current_slot);
      }
    }  // while true
  }

  /**
   * @brief Inserts the specified key/value pair into the map using scalar loads.
   *
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam CG Cooperative Group type
   *
   * @param g The Cooperative Group that performs the insert
   * @param insert_pair The pair to insert
   * @return void.
   */
  template <bool uses_vector_load, typename CG>
  __device__ __forceinline__ cuda::std::enable_if_t<not uses_vector_load, void> insert(
    CG g, value_type const& insert_pair) noexcept
  {
    auto current_slot = this->initial_slot(g, insert_pair.first);

    while (true) {
      value_type slot_contents = *reinterpret_cast<value_type const*>(current_slot);
      auto const& existing_key = slot_contents.first;

      // The user provide `key_equal` can never be used to compare against `empty_key_sentinel` as
      // the sentinel is not a valid key value. Therefore, first check for the sentinel
      auto const slot_is_empty =
        detail::bitwise_compare(existing_key, this->get_empty_key_sentinel());
      auto const bucket_contains_empty = g.ballot(slot_is_empty);

      if (bucket_contains_empty) {
        // the first lane in the group with an empty slot will attempt the insert
        insert_result status{insert_result::CONTINUE};
        uint32_t src_lane = __ffs(bucket_contains_empty) - 1;

        if (g.thread_rank() == src_lane) {
#if (__CUDA_ARCH__ < 700)
          status = cas_dependent_write(current_slot, insert_pair);
#else
          status = back_to_back_cas(current_slot, insert_pair);
#endif
        }

        // successful insert
        if (g.any(status == insert_result::SUCCESS)) { return; }
        // if we've gotten this far, a different key took our spot
        // before we could insert. We need to retry the insert on the
        // same bucket
      }
      // if there are no empty slots in the current bucket,
      // we move onto the next bucket
      else {
        current_slot = this->next_slot(current_slot);
      }
    }  // while true
  }
};  // class device_mutable_view_impl

template <typename Key,
          typename Value,
          cuda::thread_scope Scope,
          typename Allocator,
          class ProbeSequence>
class static_multimap<Key, Value, Scope, Allocator, ProbeSequence>::device_view_impl
  : public device_view_impl_base {
 public:
  using value_type     = typename device_view_impl_base::value_type;
  using key_type       = typename device_view_impl_base::key_type;
  using mapped_type    = typename device_view_impl_base::mapped_type;
  using iterator       = typename device_view_impl_base::iterator;
  using const_iterator = typename device_view_impl_base::const_iterator;

  __host__ __device__ device_view_impl(pair_atomic_type* slots,
                                       std::size_t capacity,
                                       Key empty_key_sentinel,
                                       Value empty_value_sentinel) noexcept
    : device_view_impl_base{slots, capacity, empty_key_sentinel, empty_value_sentinel}
  {
  }

  /**
   * @brief Flushes per-CG buffer into the output sequence.
   *
   * A given CUDA Cooperative Group, `g`, loads `num_outputs` key-value pairs from `output_buffer`
   * and writes them into global memory in a coalesced fashion. CG-wide `memcpy_sync` is used if
   * `thrust::is_contiguous_iterator_v<OutputIt>` returns true. All threads of `g` must be active
   * due to implicit CG-wide synchronization during flushing.
   *
   * @tparam CG Cooperative Group type
   * @tparam atomicT Type of atomic storage
   * @tparam OutputIt Device accessible output iterator whose `value_type` is
   * constructible from the map's `value_type`
   * @param g The Cooperative Group used to flush output buffer
   * @param num_outputs Number of valid output in the buffer
   * @param output_buffer Buffer of the key/value pair sequence
   * @param num_matches Size of the output sequence
   * @param output_begin Beginning of the output sequence of key/value pairs
   */
  template <typename CG, typename atomicT, typename OutputIt>
  __device__ __forceinline__ void flush_output_buffer(CG g,
                                                      uint32_t const num_outputs,
                                                      value_type* output_buffer,
                                                      atomicT* num_matches,
                                                      OutputIt output_begin) noexcept
  {
    std::size_t offset;
    const auto lane_id = g.thread_rank();
    if (0 == lane_id) {
      offset = num_matches->fetch_add(num_outputs, cuda::std::memory_order_relaxed);
    }
    offset = g.shfl(offset, 0);

    if constexpr (thrust::is_contiguous_iterator_v<OutputIt>) {
#if defined(CUCO_HAS_CUDA_BARRIER)
      cooperative_groups::memcpy_async(
        g,
        &thrust::raw_reference_cast(*(output_begin + offset)),
        output_buffer,
        cuda::aligned_size_t<alignof(value_type)>(sizeof(value_type) * num_outputs));
#else
      cooperative_groups::memcpy_async(g,
                                       &thrust::raw_reference_cast(*(output_begin + offset)),
                                       output_buffer,
                                       sizeof(value_type) * num_outputs);
#endif  // end CUCO_HAS_CUDA_BARRIER
    } else {
      for (auto index = lane_id; index < num_outputs; index += g.size()) {
        cuda::std::get<0>(*(output_begin + offset + index)) = output_buffer[index].first;
        cuda::std::get<1>(*(output_begin + offset + index)) = output_buffer[index].second;
      }
    }
  }

  /**
   * @brief Flushes per-CG buffer into the output sequences.
   *
   * A given CUDA Cooperative Group, `g`, loads `num_outputs` elements from `probe_output_buffer`
   * and `num_outputs` elements from `contained_output_buffer`, then writes them into global
   * memory started from `probe_output_begin` and `contained_output_begin` respectively. All
   * threads of `g` must be active due to implicit CG-wide synchronization during flushing.
   *
   * @tparam CG Cooperative Group type
   * @tparam atomicT Type of atomic storage
   * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
   * `InputIt`s `value_type`.
   * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
   * the map's `value_type`.
   * @param g The Cooperative Group used to flush output buffer
   * @param num_outputs Number of valid output in the buffer
   * @param probe_output_buffer Buffer of the matched probe pair sequence
   * @param contained_output_buffer Buffer of the matched contained pair sequence
   * @param num_matches Size of the output sequence
   * @param probe_output_begin Beginning of the output sequence of the matched probe pairs
   * @param contained_output_begin Beginning of the output sequence of the matched contained
   * pairs
   */
  template <typename CG, typename atomicT, typename OutputIt1, typename OutputIt2>
  __device__ __forceinline__ void flush_output_buffer(CG g,
                                                      uint32_t const num_outputs,
                                                      value_type* probe_output_buffer,
                                                      value_type* contained_output_buffer,
                                                      atomicT* num_matches,
                                                      OutputIt1 probe_output_begin,
                                                      OutputIt2 contained_output_begin) noexcept
  {
    std::size_t offset;
    const auto lane_id = g.thread_rank();
    if (0 == lane_id) {
      offset = num_matches->fetch_add(num_outputs, cuda::std::memory_order_relaxed);
    }
    offset = g.shfl(offset, 0);

    for (auto index = lane_id; index < num_outputs; index += g.size()) {
      auto& probe_pair                                          = probe_output_buffer[index];
      auto& contained_pair                                      = contained_output_buffer[index];
      cuda::std::get<0>(*(probe_output_begin + offset + index)) = probe_pair.first;
      cuda::std::get<1>(*(probe_output_begin + offset + index)) = probe_pair.second;
      cuda::std::get<0>(*(contained_output_begin + offset + index)) = contained_pair.first;
      cuda::std::get<1>(*(contained_output_begin + offset + index)) = contained_pair.second;
    }
  }

  /**
   * @brief Indicates whether the probe `element` exists in the map using vector loads.
   *
   * If `element` was inserted into the map, `contains` returns true. Otherwise, it returns false.
   * Uses the CUDA Cooperative Groups API to leverage multiple threads to perform a single
   * `contains` operation. This provides a significant boost in throughput compared to the non
   * Cooperative Group based `contains` at moderate to high load factors.
   *
   * @tparam is_pair_contains `true` if it's a `pair_contains` implementation
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam ProbeT Probe data type
   * @tparam Equal Binary callable type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g The Cooperative Group used to perform the contains operation
   * @param element The probe element to search for
   * @param equal The binary function to compare input element and slot content for equality
   * @return A boolean indicating whether the key/value pair represented by `element` was inserted
   */
  template <bool is_pair_contains,
            bool uses_vector_load,
            typename ProbeT,
            typename Equal,
            typename ParentCG>
  __device__ __forceinline__ cuda::std::enable_if_t<uses_vector_load, bool> contains(
    cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
    ProbeT const& element,
    Equal equal) const noexcept
  {
    auto current_slot = [&]() {
      if constexpr (is_pair_contains) { return initial_slot(g, element.first); }
      if constexpr (not is_pair_contains) { return initial_slot(g, element); }
    }();

    while (true) {
      value_type arr[2];
      load_pair_array(&arr[0], current_slot);

      auto const first_slot_is_empty =
        detail::bitwise_compare(arr[0].first, this->get_empty_key_sentinel());
      auto const second_slot_is_empty =
        detail::bitwise_compare(arr[1].first, this->get_empty_key_sentinel());
      auto const first_equals = [&]() {
        if constexpr (is_pair_contains) {
          return not first_slot_is_empty and equal(arr[0], element);
        }
        if constexpr (not is_pair_contains) {
          return not first_slot_is_empty and equal(arr[0].first, element);
        }
      }();
      auto const second_equals = [&]() {
        if constexpr (is_pair_contains) {
          return not second_slot_is_empty and equal(arr[1], element);
        }
        if constexpr (not is_pair_contains) {
          return not second_slot_is_empty and equal(arr[1].first, element);
        }
      }();

      // the key we were searching for was found by one of the threads, so we return true
      if (g.any(first_equals or second_equals)) { return true; }

      // we found an empty slot, meaning that the key we're searching for isn't present
      if (g.any(first_slot_is_empty or second_slot_is_empty)) { return false; }

      // otherwise, all slots in the current bucket are full with other keys, so we move onto the
      // next bucket
      current_slot = next_slot(current_slot);
    }
  }

  /**
   * @brief Indicates whether `element` exists in the map using scalar loads.
   *
   * If `element` was inserted into the map, `contains` returns true. Otherwise, it returns false.
   * Uses the CUDA Cooperative Groups API to leverage multiple threads to perform a single
   * `contains` operation. This provides a significant boost in throughput compared to the non
   * Cooperative Group `contains` at moderate to high load factors.
   *
   * @tparam is_pair_contains `true` if it's a `pair_contains` implementation
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam ProbeT Probe data type
   * @tparam Equal Binary callable type
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g The Cooperative Group used to perform the contains operation
   * @param element The probe element to search for
   * @param equal The binary function to compare input element and slot content for equality
   * @return A boolean indicating whether the key/value pair represented by `element` was inserted
   */
  template <bool is_pair_contains,
            bool uses_vector_load,
            typename ProbeT,
            typename Equal,
            typename ParentCG>
  __device__ __forceinline__ cuda::std::enable_if_t<not uses_vector_load, bool> contains(
    cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
    ProbeT const& element,
    Equal equal) const noexcept
  {
    auto current_slot = [&]() {
      if constexpr (is_pair_contains) { return this->initial_slot(g, element.first); }
      if constexpr (not is_pair_contains) { return this->initial_slot(g, element); }
    }();

    while (true) {
      value_type slot_contents = *reinterpret_cast<value_type const*>(current_slot);
      auto const& existing_key = slot_contents.first;

      // The user provide `key_equal` can never be used to compare against `empty_key_sentinel` as
      // the sentinel is not a valid key value. Therefore, first check for the sentinel
      auto const slot_is_empty =
        detail::bitwise_compare(existing_key, this->get_empty_key_sentinel());

      auto const equals = [&]() {
        if constexpr (is_pair_contains) {
          return not slot_is_empty and equal(slot_contents, element);
        }
        if constexpr (not is_pair_contains) {
          return not slot_is_empty and equal(existing_key, element);
        }
      }();

      // the key we were searching for was found by one of the threads, so we return true
      if (g.any(equals)) { return true; }

      // we found an empty slot, meaning that the key we're searching for isn't present
      if (g.any(slot_is_empty)) { return false; }

      // otherwise, all slots in the current bucket are full with other keys, so we move onto the
      // next bucket
      current_slot = this->next_slot(current_slot);
    }
  }

  /**
   * @brief Counts the occurrence of a given key contained in multimap using vector loads.
   *
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam CG Cooperative Group type
   * @tparam KeyEqual Binary callable type
   * @param g The Cooperative Group used to perform the count operation
   * @param k The key to search for
   * @param key_equal The binary callable used to compare two keys
   * for equality
   * @return Number of matches found by the current thread
   */
  template <bool uses_vector_load, bool is_outer, typename CG, typename KeyEqual>
  __device__ __forceinline__ cuda::std::enable_if_t<uses_vector_load, std::size_t> count(
    CG g, Key const& k, KeyEqual key_equal) noexcept
  {
    std::size_t count = 0;
    auto current_slot = initial_slot(g, k);

    [[maybe_unused]] bool found_match = false;

    while (true) {
      value_type arr[2];
      load_pair_array(&arr[0], current_slot);

      auto const first_slot_is_empty =
        detail::bitwise_compare(arr[0].first, this->get_empty_key_sentinel());
      auto const second_slot_is_empty =
        detail::bitwise_compare(arr[1].first, this->get_empty_key_sentinel());
      auto const first_equals  = (not first_slot_is_empty and key_equal(arr[0].first, k));
      auto const second_equals = (not second_slot_is_empty and key_equal(arr[1].first, k));

      if constexpr (is_outer) {
        if (g.any(first_equals or second_equals)) { found_match = true; }
      }

      count += (first_equals + second_equals);

      if (g.any(first_slot_is_empty or second_slot_is_empty)) {
        if constexpr (is_outer) {
          if ((not found_match) && (g.thread_rank() == 0)) { count++; }
        }
        return count;
      }

      current_slot = next_slot(current_slot);
    }
  }

  /**
   * @brief Counts the occurrence of a given key contained in multimap using scalar loads.
   *
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam CG Cooperative Group type
   * @tparam KeyEqual Binary callable type
   * @param g The Cooperative Group used to perform the count operation
   * @param k The key to search for
   * @param key_equal The binary callable used to compare two keys
   * for equality
   * @return Number of matches found by the current thread
   */
  template <bool uses_vector_load, bool is_outer, typename CG, typename KeyEqual>
  __device__ __forceinline__ cuda::std::enable_if_t<not uses_vector_load, std::size_t> count(
    CG g, Key const& k, KeyEqual key_equal) noexcept
  {
    std::size_t count = 0;
    auto current_slot = initial_slot(g, k);

    [[maybe_unused]] bool found_match = false;

    while (true) {
      value_type slot_contents = *reinterpret_cast<value_type const*>(current_slot);
      auto const& current_key  = slot_contents.first;

      auto const slot_is_empty =
        detail::bitwise_compare(current_key, this->get_empty_key_sentinel());
      auto const equals = not slot_is_empty and key_equal(current_key, k);

      if constexpr (is_outer) {
        if (g.any(equals)) { found_match = true; }
      }

      count += equals;

      if (g.any(slot_is_empty)) {
        if constexpr (is_outer) {
          if ((not found_match) && (g.thread_rank() == 0)) { count++; }
        }
        return count;
      }

      current_slot = next_slot(current_slot);
    }
  }

  /**
   * @brief Counts the occurrence of a given key/value pair contained in multimap using vector
   * loads.
   *
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam CG Cooperative Group type
   * @tparam PairEqual Binary callable type
   * @param g The Cooperative Group used to perform the pair_count operation
   * @param pair The pair to search for
   * @param pair_equal The binary callable used to compare two pairs
   * for equality
   * @return Number of matches found by the current thread
   */
  template <bool uses_vector_load, bool is_outer, typename CG, typename PairEqual>
  __device__ __forceinline__ cuda::std::enable_if_t<uses_vector_load, std::size_t> pair_count(
    CG g, value_type const& pair, PairEqual pair_equal) noexcept
  {
    std::size_t count = 0;
    auto key          = pair.first;
    auto current_slot = initial_slot(g, key);

    [[maybe_unused]] bool found_match = false;

    while (true) {
      value_type arr[2];
      load_pair_array(&arr[0], current_slot);

      auto const first_slot_is_empty =
        detail::bitwise_compare(arr[0].first, this->get_empty_key_sentinel());
      auto const second_slot_is_empty =
        detail::bitwise_compare(arr[1].first, this->get_empty_key_sentinel());

      auto const first_slot_equals  = (not first_slot_is_empty and pair_equal(arr[0], pair));
      auto const second_slot_equals = (not second_slot_is_empty and pair_equal(arr[1], pair));

      if constexpr (is_outer) {
        if (g.any(first_slot_equals or second_slot_equals)) { found_match = true; }
      }

      count += (first_slot_equals + second_slot_equals);

      if (g.any(first_slot_is_empty or second_slot_is_empty)) {
        if constexpr (is_outer) {
          if ((not found_match) && (g.thread_rank() == 0)) { count++; }
        }
        return count;
      }

      current_slot = next_slot(current_slot);
    }
  }

  /**
   * @brief Counts the occurrence of a given key/value pair contained in multimap using scalar
   * loads.
   *
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam CG Cooperative Group type
   * @tparam PairEqual Binary callable type
   * @param g The Cooperative Group used to perform the pair_count operation
   * @param pair The pair to search for
   * @param pair_equal The binary callable used to compare two pairs
   * for equality
   * @return Number of matches found by the current thread
   */
  template <bool uses_vector_load, bool is_outer, typename CG, typename PairEqual>
  __device__ __forceinline__ cuda::std::enable_if_t<not uses_vector_load, std::size_t> pair_count(
    CG g, value_type const& pair, PairEqual pair_equal) noexcept
  {
    std::size_t count = 0;
    auto key          = pair.first;
    auto current_slot = initial_slot(g, key);

    [[maybe_unused]] bool found_match = false;

    while (true) {
      auto slot_contents = *reinterpret_cast<value_type const*>(current_slot);

      auto const slot_is_empty =
        detail::bitwise_compare(slot_contents.first, this->get_empty_key_sentinel());

      auto const equals = not slot_is_empty and pair_equal(slot_contents, pair);

      if constexpr (is_outer) {
        if (g.any(equals)) { found_match = true; }
      }

      count += equals;

      if (g.any(slot_is_empty)) {
        if constexpr (is_outer) {
          if ((not found_match) && (g.thread_rank() == 0)) { count++; }
        }
        return count;
      }

      current_slot = next_slot(current_slot);
    }
  }

  /**
   * @brief Retrieves all the matches of a given key contained in multimap using vector
   * loads with per-flushing-CG shared memory buffer.
   *
   * For key `k` existing in the map, copies `k` and all associated values to unspecified
   * locations in `[output_begin, output_end)`. If `k` does not have any matches, copies `k` and
   * `empty_value_sentinel()` into the output only if `is_outer` is true.
   *
   * @tparam buffer_size Size of the output buffer
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam FlushingCG Type of Cooperative Group used to flush output buffer
   * @tparam ProbingCG Type of Cooperative Group used to retrieve
   * @tparam atomicT Type of atomic storage
   * @tparam OutputIt Device accessible output iterator whose `value_type` is
   * constructible from the map's `value_type`
   * @tparam KeyEqual Binary callable type
   * @param flushing_cg The Cooperative Group used to flush output buffer
   * @param probing_cg The Cooperative Group used to retrieve
   * @param k The key to search for
   * @param flushing_cg_counter Pointer to the flushing cg counter
   * @param output_buffer Shared memory buffer of the key/value pair sequence
   * @param num_matches Size of the output sequence
   * @param output_begin Beginning of the output sequence of key/value pairs
   * @param key_equal The binary callable used to compare two keys
   * for equality
   */
  template <uint32_t buffer_size,
            bool is_outer,
            typename FlushingCG,
            typename ProbingCG,
            typename atomicT,
            typename OutputIt,
            typename KeyEqual>
  __device__ __forceinline__ void retrieve(FlushingCG flushing_cg,
                                           ProbingCG probing_cg,
                                           Key const& k,
                                           uint32_t* flushing_cg_counter,
                                           value_type* output_buffer,
                                           atomicT* num_matches,
                                           OutputIt output_begin,
                                           KeyEqual key_equal) noexcept
  {
    const uint32_t cg_lane_id = probing_cg.thread_rank();

    auto current_slot = initial_slot(probing_cg, k);

    bool running                      = true;
    [[maybe_unused]] bool found_match = false;

    while (flushing_cg.any(running)) {
      if (running) {
        value_type arr[2];
        load_pair_array(&arr[0], current_slot);

        auto const first_slot_is_empty =
          detail::bitwise_compare(arr[0].first, this->get_empty_key_sentinel());
        auto const second_slot_is_empty =
          detail::bitwise_compare(arr[1].first, this->get_empty_key_sentinel());
        auto const first_equals  = (not first_slot_is_empty and key_equal(arr[0].first, k));
        auto const second_equals = (not second_slot_is_empty and key_equal(arr[1].first, k));
        auto const first_exists  = probing_cg.ballot(first_equals);
        auto const second_exists = probing_cg.ballot(second_equals);

        if (first_exists or second_exists) {
          if constexpr (is_outer) { found_match = true; }

          auto const num_first_matches  = __popc(first_exists);
          auto const num_second_matches = __popc(second_exists);

          uint32_t output_idx;
          if (0 == cg_lane_id) {
            output_idx = atomicAdd(flushing_cg_counter, (num_first_matches + num_second_matches));
          }
          output_idx = probing_cg.shfl(output_idx, 0);

          if (first_equals) {
            auto const lane_offset = detail::count_least_significant_bits(first_exists, cg_lane_id);
            output_buffer[output_idx + lane_offset] = cuco::make_pair(k, arr[0].second);
          }
          if (second_equals) {
            auto const lane_offset =
              detail::count_least_significant_bits(second_exists, cg_lane_id);
            output_buffer[output_idx + num_first_matches + lane_offset] =
              cuco::make_pair(k, arr[1].second);
          }
        }
        if (probing_cg.any(first_slot_is_empty or second_slot_is_empty)) {
          running = false;
          if constexpr (is_outer) {
            if ((not found_match) && (cg_lane_id == 0)) {
              auto const output_idx     = atomicAdd(flushing_cg_counter, 1);
              output_buffer[output_idx] = cuco::make_pair(k, this->get_empty_value_sentinel());
            }
          }
        }
      }  // if running

      flushing_cg.sync();
      if (*flushing_cg_counter + flushing_cg.size() * vector_width() > buffer_size) {
        flush_output_buffer(
          flushing_cg, *flushing_cg_counter, output_buffer, num_matches, output_begin);
        // Everyone in the group reads the counter when flushing, so
        // sync before writing.
        flushing_cg.sync();
        // First lane reset warp-level counter
        if (flushing_cg.thread_rank() == 0) { *flushing_cg_counter = 0; }
        flushing_cg.sync();
      }

      current_slot = next_slot(current_slot);
    }  // while running
  }

  /**
   * @brief Retrieves all the matches of a given key contained in multimap using scalar
   * loads with per-CG shared memory buffer.
   *
   * For key `k` existing in the map, copies `k` and all associated values to unspecified
   * locations in `[output_begin, output_end)`. If `k` does not have any matches, copies `k` and
   * `empty_value_sentinel()` into the output only if `is_outer` is true.
   *
   * @tparam buffer_size Size of the output buffer
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam CG Cooperative Group type
   * @tparam atomicT Type of atomic storage
   * @tparam OutputIt Device accessible output iterator whose `value_type` is
   * constructible from the map's `value_type`
   * @tparam KeyEqual Binary callable type
   * @param g The Cooperative Group used to retrieve
   * @param k The key to search for
   * @param cg_counter Pointer to the CG counter
   * @param output_buffer Shared memory buffer of the key/value pair sequence
   * @param num_matches Size of the output sequence
   * @param output_begin Beginning of the output sequence of key/value pairs
   * @param key_equal The binary callable used to compare two keys
   * for equality
   */
  template <uint32_t buffer_size,
            bool is_outer,
            typename CG,
            typename atomicT,
            typename OutputIt,
            typename KeyEqual>
  __device__ __forceinline__ void retrieve(CG g,
                                           Key const& k,
                                           uint32_t* cg_counter,
                                           value_type* output_buffer,
                                           atomicT* num_matches,
                                           OutputIt output_begin,
                                           KeyEqual key_equal) noexcept
  {
    const uint32_t lane_id = g.thread_rank();

    auto current_slot = initial_slot(g, k);

    bool running                      = true;
    [[maybe_unused]] bool found_match = false;

    while (running) {
      // TODO: Replace reinterpret_cast with atomic ref when possible. The current implementation
      // is unsafe!
      static_assert(sizeof(Key) == sizeof(cuda::atomic<Key>));
      static_assert(sizeof(Value) == sizeof(cuda::atomic<Value>));
      value_type slot_contents = *reinterpret_cast<value_type const*>(current_slot);

      auto const slot_is_empty =
        detail::bitwise_compare(slot_contents.first, this->get_empty_key_sentinel());
      auto const equals = (not slot_is_empty and key_equal(slot_contents.first, k));
      auto const exists = g.ballot(equals);

      uint32_t output_idx = *cg_counter;

      if (exists) {
        if constexpr (is_outer) { found_match = true; }
        auto const num_matches = __popc(exists);
        if (equals) {
          // Each match computes its lane-level offset
          auto const lane_offset = detail::count_least_significant_bits(exists, lane_id);
          output_buffer[output_idx + lane_offset] = cuco::make_pair(k, slot_contents.second);
        }
        if (0 == lane_id) { (*cg_counter) += num_matches; }
      }
      if (g.any(slot_is_empty)) {
        running = false;
        if constexpr (is_outer) {
          if ((not found_match) && (lane_id == 0)) {
            output_idx                = (*cg_counter)++;
            output_buffer[output_idx] = cuco::make_pair(k, this->get_empty_value_sentinel());
          }
        }
      }

      g.sync();

      // Flush if the next iteration won't fit into buffer
      if ((*cg_counter + g.size()) > buffer_size) {
        flush_output_buffer(g, *cg_counter, output_buffer, num_matches, output_begin);
        // Everyone in the group reads the counter when flushing, so
        // sync before writing.
        g.sync();
        // First lane reset CG-level counter
        if (lane_id == 0) { *cg_counter = 0; }
        g.sync();
      }
      current_slot = next_slot(current_slot);
    }  // while running
  }

  /**
   * @brief Retrieves all the matches of a given pair using vector loads.
   *
   * For pair `p` with `n` matching pairs, if `pair_equal(p, slot)` returns true, stores
   * `probe_key_begin[j] = p.first`, `probe_val_begin[j] = p.second`, `contained_key_begin[j] =
   * slot.first`, and `contained_val_begin[j] = slot.second` for an unspecified value of `j` where
   * `0 <= j < n`. If `p` does not have any matches, stores `probe_key_begin[0] = p.first`,
   * `probe_val_begin[0] = p.second`, `contained_key_begin[0] = empty_key_sentinel`, and
   * `contained_val_begin[0] = empty_value_sentinel` only if `is_outer` is true.
   *
   * Concurrent reads or writes to any of the output ranges results in undefined behavior.
   *
   * Behavior is undefined if the extent of any of the output ranges is less than `n`.
   *
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam ProbingCG Type of Cooperative Group used to retrieve
   * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
   * `pair`'s `Key` type.
   * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
   * `pair`'s `Value` type.
   * @tparam OutputIt3 Device accessible output iterator whose `value_type` is constructible from
   * the map's `key_type`.
   * @tparam OutputIt4 Device accessible output iterator whose `value_type` is constructible from
   * the map's `mapped_type`.
   * @tparam PairEqual Binary callable type
   * @param probing_cg The Cooperative Group used to retrieve
   * @param pair The pair to search for
   * @param probe_key_begin Beginning of the output sequence of the matched probe keys
   * @param probe_val_begin Beginning of the output sequence of the matched probe values
   * @param contained_key_begin Beginning of the output sequence of the matched contained keys
   * @param contained_val_begin Beginning of the output sequence of the matched contained values
   * @param pair_equal The binary callable used to compare two pairs for equality
   */
  template <bool is_outer,
            bool uses_vector_load,
            typename ProbingCG,
            typename OutputIt1,
            typename OutputIt2,
            typename OutputIt3,
            typename OutputIt4,
            typename PairEqual>
  __device__ __forceinline__ cuda::std::enable_if_t<uses_vector_load, void> pair_retrieve(
    ProbingCG probing_cg,
    value_type const& pair,
    OutputIt1 probe_key_begin,
    OutputIt2 probe_val_begin,
    OutputIt3 contained_key_begin,
    OutputIt4 contained_val_begin,
    PairEqual pair_equal) noexcept
  {
    auto const lane_id                = probing_cg.thread_rank();
    auto current_slot                 = initial_slot(probing_cg, pair.first);
    [[maybe_unused]] auto found_match = false;

    auto num_matches = 0;

    while (true) {
      value_type arr[2];
      load_pair_array(&arr[0], current_slot);

      auto const first_slot_is_empty =
        detail::bitwise_compare(arr[0].first, this->get_empty_key_sentinel());
      auto const second_slot_is_empty =
        detail::bitwise_compare(arr[1].first, this->get_empty_key_sentinel());
      auto const first_equals  = (not first_slot_is_empty and pair_equal(arr[0], pair));
      auto const second_equals = (not second_slot_is_empty and pair_equal(arr[1], pair));
      auto const first_exists  = probing_cg.ballot(first_equals);
      auto const second_exists = probing_cg.ballot(second_equals);

      if (first_exists or second_exists) {
        if constexpr (is_outer) { found_match = true; }

        auto const num_first_matches = __popc(first_exists);

        if (first_equals) {
          auto lane_offset      = detail::count_least_significant_bits(first_exists, lane_id);
          auto const output_idx = num_matches + lane_offset;

          *(probe_key_begin + output_idx)     = pair.first;
          *(probe_val_begin + output_idx)     = pair.second;
          *(contained_key_begin + output_idx) = arr[0].first;
          *(contained_val_begin + output_idx) = arr[0].second;
        }
        if (second_equals) {
          auto const lane_offset = detail::count_least_significant_bits(second_exists, lane_id);
          auto const output_idx  = num_matches + num_first_matches + lane_offset;

          *(probe_key_begin + output_idx)     = pair.first;
          *(probe_val_begin + output_idx)     = pair.second;
          *(contained_key_begin + output_idx) = arr[1].first;
          *(contained_val_begin + output_idx) = arr[1].second;
        }
        num_matches += (num_first_matches + __popc(second_exists));
      }
      if (probing_cg.any(first_slot_is_empty or second_slot_is_empty)) {
        if constexpr (is_outer) {
          if ((not found_match) and lane_id == 0) {
            *(probe_key_begin)     = pair.first;
            *(probe_val_begin)     = pair.second;
            *(contained_key_begin) = this->get_empty_key_sentinel();
            *(contained_val_begin) = this->get_empty_value_sentinel();
          }
        }
        return;  // exit if any slot in the current bucket is empty
      }

      current_slot = next_slot(current_slot);
    }  // while
  }

  /**
   * @brief Retrieves all the matches of a given pair using scalar loads.
   *
   * For pair `p` with `n` matching pairs, if `pair_equal(p, slot)` returns true, stores
   * `probe_key_begin[j] = p.first`, `probe_val_begin[j] = p.second`, `contained_key_begin[j] =
   * slot.first`, and `contained_val_begin[j] = slot.second` for an unspecified value of `j` where
   * `0 <= j < n`. If `p` does not have any matches, stores `probe_key_begin[0] = p.first`,
   * `probe_val_begin[0] = p.second`, `contained_key_begin[0] = empty_key_sentinel`, and
   * `contained_val_begin[0] = empty_value_sentinel` only if `is_outer` is true.
   *
   * Concurrent reads or writes to any of the output ranges results in undefined behavior.
   *
   * Behavior is undefined if the extent of any of the output ranges is less than `n`.
   *
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam uses_vector_load Boolean flag indicating whether vector loads are used
   * @tparam ProbingCG Type of Cooperative Group used to retrieve
   * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
   * `pair`'s `Key` type.
   * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
   * `pair`'s `Value` type.
   * @tparam OutputIt3 Device accessible output iterator whose `value_type` is constructible from
   * the map's `key_type`.
   * @tparam OutputIt4 Device accessible output iterator whose `value_type` is constructible from
   * the map's `mapped_type`.
   * @tparam PairEqual Binary callable type
   * @param probing_cg The Cooperative Group used to retrieve
   * @param pair The pair to search for
   * @param probe_key_begin Beginning of the output sequence of the matched probe keys
   * @param probe_val_begin Beginning of the output sequence of the matched probe values
   * @param contained_key_begin Beginning of the output sequence of the matched contained keys
   * @param contained_val_begin Beginning of the output sequence of the matched contained values
   * @param pair_equal The binary callable used to compare two pairs for equality
   */
  template <bool is_outer,
            bool uses_vector_load,
            typename ProbingCG,
            typename OutputIt1,
            typename OutputIt2,
            typename OutputIt3,
            typename OutputIt4,
            typename PairEqual>
  __device__ __forceinline__ cuda::std::enable_if_t<not uses_vector_load, void> pair_retrieve(
    ProbingCG probing_cg,
    value_type const& pair,
    OutputIt1 probe_key_begin,
    OutputIt2 probe_val_begin,
    OutputIt3 contained_key_begin,
    OutputIt4 contained_val_begin,
    PairEqual pair_equal) noexcept
  {
    auto const lane_id                = probing_cg.thread_rank();
    auto current_slot                 = initial_slot(probing_cg, pair.first);
    [[maybe_unused]] auto found_match = false;

    auto num_matches = 0;

    while (true) {
      // TODO: Replace reinterpret_cast with atomic ref when possible. The current implementation
      // is unsafe!
      static_assert(sizeof(Key) == sizeof(cuda::atomic<Key>));
      static_assert(sizeof(Value) == sizeof(cuda::atomic<Value>));
      value_type slot_contents = *reinterpret_cast<value_type const*>(current_slot);

      auto const slot_is_empty =
        detail::bitwise_compare(slot_contents.first, this->get_empty_key_sentinel());
      auto const equals = (not slot_is_empty and pair_equal(slot_contents, pair));
      auto const exists = probing_cg.ballot(equals);

      if (exists) {
        if constexpr (is_outer) { found_match = true; }

        if (equals) {
          auto const lane_offset = detail::count_least_significant_bits(exists, lane_id);
          auto const output_idx  = num_matches + lane_offset;

          *(probe_key_begin + output_idx)     = pair.first;
          *(probe_val_begin + output_idx)     = pair.second;
          *(contained_key_begin + output_idx) = slot_contents.first;
          *(contained_val_begin + output_idx) = slot_contents.second;
        }
        num_matches += __popc(exists);
      }
      if (probing_cg.any(slot_is_empty)) {
        if constexpr (is_outer) {
          if ((not found_match) and lane_id == 0) {
            *(probe_key_begin)     = pair.first;
            *(probe_val_begin)     = pair.second;
            *(contained_key_begin) = this->get_empty_key_sentinel();
            *(contained_val_begin) = this->get_empty_value_sentinel();
          }
        }
        return;  // exit if any slot in the current bucket is empty
      }

      current_slot = next_slot(current_slot);
    }  // while
  }

  /**
   * @brief Retrieves all the matches of a given pair contained in multimap using vector
   * loads with per-flushing-CG shared memory buffer.
   *
   * For pair `p`, if pair_equal(p, slot[j]) returns true, copies `p` to unspecified locations
   * in `[probe_output_begin, probe_output_end)` and copies slot[j] to unspecified locations in
   * `[contained_output_begin, contained_output_end)`. If `p` does not have any matches, copies
   * `p` and a pair of `empty_key_sentinel` and `empty_value_sentinel` into the output only if
   * `is_outer` is true.
   *
   * @tparam buffer_size Size of the output buffer
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam FlushingCG Type of Cooperative Group used to flush output buffer
   * @tparam ProbingCG Type of Cooperative Group used to retrieve
   * @tparam atomicT Type of atomic storage
   * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
   * `InputIt`s `value_type`.
   * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
   * the map's `value_type`.
   * @tparam PairEqual Binary callable type
   * @param flushing_cg The Cooperative Group used to flush output buffer
   * @param probing_cg The Cooperative Group used to retrieve
   * @param pair The pair to search for
   * @param flushing_cg_counter Pointer to the flushing CG counter
   * @param probe_output_buffer Buffer of the matched probe pair sequence
   * @param contained_output_buffer Buffer of the matched contained pair sequence
   * @param num_matches Size of the output sequence
   * @param probe_output_begin Beginning of the output sequence of the matched probe pairs
   * @param contained_output_begin Beginning of the output sequence of the matched contained
   * pairs
   * @param pair_equal The binary callable used to compare two pairs for equality
   */
  template <uint32_t buffer_size,
            bool is_outer,
            typename FlushingCG,
            typename ProbingCG,
            typename atomicT,
            typename OutputIt1,
            typename OutputIt2,
            typename PairEqual>
  __device__ __forceinline__ void pair_retrieve(FlushingCG flushing_cg,
                                                ProbingCG probing_cg,
                                                value_type const& pair,
                                                uint32_t* flushing_cg_counter,
                                                value_type* probe_output_buffer,
                                                value_type* contained_output_buffer,
                                                atomicT* num_matches,
                                                OutputIt1 probe_output_begin,
                                                OutputIt2 contained_output_begin,
                                                PairEqual pair_equal) noexcept
  {
    const uint32_t cg_lane_id = probing_cg.thread_rank();

    auto key          = pair.first;
    auto current_slot = initial_slot(probing_cg, key);

    bool running                      = true;
    [[maybe_unused]] bool found_match = false;

    while (flushing_cg.any(running)) {
      if (running) {
        value_type arr[2];
        load_pair_array(&arr[0], current_slot);

        auto const first_slot_is_empty =
          detail::bitwise_compare(arr[0].first, this->get_empty_key_sentinel());
        auto const second_slot_is_empty =
          detail::bitwise_compare(arr[1].first, this->get_empty_key_sentinel());
        auto const first_equals  = (not first_slot_is_empty and pair_equal(arr[0], pair));
        auto const second_equals = (not second_slot_is_empty and pair_equal(arr[1], pair));
        auto const first_exists  = probing_cg.ballot(first_equals);
        auto const second_exists = probing_cg.ballot(second_equals);

        if (first_exists or second_exists) {
          if constexpr (is_outer) { found_match = true; }

          auto const num_first_matches  = __popc(first_exists);
          auto const num_second_matches = __popc(second_exists);

          uint32_t output_idx;
          if (0 == cg_lane_id) {
            output_idx = atomicAdd(flushing_cg_counter, (num_first_matches + num_second_matches));
          }
          output_idx = probing_cg.shfl(output_idx, 0);

          if (first_equals) {
            auto const lane_offset = detail::count_least_significant_bits(first_exists, cg_lane_id);
            probe_output_buffer[output_idx + lane_offset]     = pair;
            contained_output_buffer[output_idx + lane_offset] = arr[0];
          }
          if (second_equals) {
            auto const lane_offset =
              detail::count_least_significant_bits(second_exists, cg_lane_id);
            probe_output_buffer[output_idx + num_first_matches + lane_offset]     = pair;
            contained_output_buffer[output_idx + num_first_matches + lane_offset] = arr[1];
          }
        }
        if (probing_cg.any(first_slot_is_empty or second_slot_is_empty)) {
          running = false;
          if constexpr (is_outer) {
            if ((not found_match) && (cg_lane_id == 0)) {
              auto const output_idx           = atomicAdd(flushing_cg_counter, 1);
              probe_output_buffer[output_idx] = pair;
              contained_output_buffer[output_idx] =
                cuco::make_pair(this->get_empty_key_sentinel(), this->get_empty_value_sentinel());
            }
          }
        }
      }  // if running

      flushing_cg.sync();
      if (*flushing_cg_counter + flushing_cg.size() * vector_width() > buffer_size) {
        flush_output_buffer(flushing_cg,
                            *flushing_cg_counter,
                            probe_output_buffer,
                            contained_output_buffer,
                            num_matches,
                            probe_output_begin,
                            contained_output_begin);
        // Everyone in the group reads the counter when flushing, so
        // sync before writing.
        flushing_cg.sync();
        // First lane reset warp-level counter
        if (flushing_cg.thread_rank() == 0) { *flushing_cg_counter = 0; }
        flushing_cg.sync();
      }

      current_slot = next_slot(current_slot);
    }  // while running
  }

  /**
   * @brief Retrieves all the matches of a given pair contained in multimap using scalar
   * loads with per-CG shared memory buffer.
   *
   * For pair `p`, if pair_equal(p, slot[j]) returns true, copies `p` to unspecified locations
   * in `[probe_output_begin, probe_output_end)` and copies slot[j] to unspecified locations in
   * `[contained_output_begin, contained_output_end)`. If `p` does not have any matches, copies
   * `p` and a pair of `empty_key_sentinel` and `empty_value_sentinel` into the output only if
   * `is_outer` is true.
   *
   * @tparam buffer_size Size of the output buffer
   * @tparam is_outer Boolean flag indicating whether outer join is peformed
   * @tparam CG Cooperative Group type
   * @tparam atomicT Type of atomic storage
   * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
   * `InputIt`s `value_type`.
   * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
   * the map's `value_type`.
   * @tparam PairEqual Binary callable type
   * @param g The Cooperative Group used to retrieve
   * @param pair The pair to search for
   * @param cg_counter Pointer to the CG counter
   * @param probe_output_buffer Buffer of the matched probe pair sequence
   * @param contained_output_buffer Buffer of the matched contained pair sequence
   * @param num_matches Size of the output sequence
   * @param probe_output_begin Beginning of the output sequence of the matched probe pairs
   * @param contained_output_begin Beginning of the output sequence of the matched contained
   * pairs
   * @param pair_equal The binary callable used to compare two pairs for equality
   */
  template <uint32_t buffer_size,
            bool is_outer,
            typename CG,
            typename atomicT,
            typename OutputIt1,
            typename OutputIt2,
            typename PairEqual>
  __device__ __forceinline__ void pair_retrieve(CG g,
                                                value_type const& pair,
                                                uint32_t* cg_counter,
                                                value_type* probe_output_buffer,
                                                value_type* contained_output_buffer,
                                                atomicT* num_matches,
                                                OutputIt1 probe_output_begin,
                                                OutputIt2 contained_output_begin,
                                                PairEqual pair_equal) noexcept
  {
    const uint32_t lane_id = g.thread_rank();

    auto key          = pair.first;
    auto current_slot = initial_slot(g, key);

    bool running                      = true;
    [[maybe_unused]] bool found_match = false;

    while (running) {
      // TODO: Replace reinterpret_cast with atomic ref when possible. The current implementation
      // is unsafe!
      static_assert(sizeof(Key) == sizeof(cuda::atomic<Key>));
      static_assert(sizeof(Value) == sizeof(cuda::atomic<Value>));
      value_type slot_contents = *reinterpret_cast<value_type const*>(current_slot);

      auto const slot_is_empty =
        detail::bitwise_compare(slot_contents.first, this->get_empty_key_sentinel());
      auto const equals = (not slot_is_empty and pair_equal(slot_contents, pair));
      auto const exists = g.ballot(equals);

      uint32_t output_idx = *cg_counter;

      if (exists) {
        if constexpr (is_outer) { found_match = true; }
        auto const num_matches = __popc(exists);
        if (equals) {
          // Each match computes its lane-level offset
          auto const lane_offset = detail::count_least_significant_bits(exists, lane_id);
          probe_output_buffer[output_idx + lane_offset]     = pair;
          contained_output_buffer[output_idx + lane_offset] = slot_contents;
        }
        if (0 == lane_id) { (*cg_counter) += num_matches; }
      }
      if (g.any(slot_is_empty)) {
        running = false;
        if constexpr (is_outer) {
          if ((not found_match) && (lane_id == 0)) {
            output_idx                      = (*cg_counter)++;
            probe_output_buffer[output_idx] = pair;
            contained_output_buffer[output_idx] =
              cuco::make_pair(this->get_empty_key_sentinel(), this->get_empty_value_sentinel());
          }
        }
      }

      g.sync();

      // Flush if the next iteration won't fit into buffer
      if ((*cg_counter + g.size()) > buffer_size) {
        flush_output_buffer(g,
                            *cg_counter,
                            probe_output_buffer,
                            contained_output_buffer,
                            num_matches,
                            probe_output_begin,
                            contained_output_begin);
        // Everyone in the group reads the counter when flushing, so
        // sync before writing.
        g.sync();
        // First lane reset CG-level counter
        if (lane_id == 0) { *cg_counter = 0; }
        g.sync();
      }
      current_slot = next_slot(current_slot);
    }  // while running
  }
};  // class device_view_impl

}  // namespace cuco
