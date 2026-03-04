/*
 * SPDX-FileCopyrightText: Copyright (c) 2023-2025 NVIDIA CORPORATION & AFFILIATES. All rights
 * reserved. SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#pragma once

#include <cuda/std/array>
#include <cuda/std/cstddef>
#include <cuda/stream_ref>
#include <thrust/device_malloc_allocator.h>
#include <thrust/device_vector.h>

#include <climits>
#include <cstddef>

namespace cuco {
namespace experimental {
namespace detail {

/**
 * @brief Struct to store ranks of bits at 256-bit intervals (or blocks)
 *
 * This struct encodes a list of four rank values using base + offset format
 * e.g. [1000, 1005, 1006, 1009] is stored as base = 1000, offsets = [5, 6, 9]
 * base uses 40 bits, split between one uint32_t and one uint8_t
 * each offset uses 8 bits
 */
struct rank {
  uint32_t base_hi_;                      ///< Upper 32 bits of base
  uint8_t base_lo_;                       ///< Lower 8 bits of base
  cuda::std::array<uint8_t, 3> offsets_;  ///< Offsets for 64-bit sub-intervals, relative to base

  /**
   * @brief Gets base rank of current 256-bit interval
   *
   * @return The base rank
   */
  __host__ __device__ constexpr uint64_t base() const noexcept
  {
    return (static_cast<uint64_t>(base_hi_) << CHAR_BIT) | base_lo_;
  }

  /**
   * @brief Sets base rank of current 256-bit interval
   *
   * @param base Base rank
   */
  __host__ __device__ constexpr void set_base(uint64_t base) noexcept
  {
    base_hi_ = static_cast<uint32_t>(base >> CHAR_BIT);
    base_lo_ = static_cast<uint8_t>(base);
  }
};

/**
 * @brief Bitset class with rank and select index structures
 *
 * In addition to standard bitset set/test operations, this class provides
 * rank and select operation API. It maintains index structures to make both these
 * new operations close to constant time.
 *
 * Current limitations:
 * - Stream controls are partially supported due to the use of `thrust::device_vector` as storage
 * - Device ref doesn't support modifiers like `set`, `reset`, etc.
 *
 * @tparam Allocator Type of allocator used for device storage
 */
// TODO: have to use device_malloc_allocator for now otherwise the container cannot grow
template <class Allocator = thrust::device_malloc_allocator<cuda::std::byte>>
class dynamic_bitset {
 public:
  using size_type = std::size_t;  ///< size type to specify bit index
  using word_type = uint64_t;     ///< word type
  /// Type of the allocator to (de)allocate words
  using allocator_type =
    typename std::allocator_traits<Allocator>::template rebind_alloc<word_type>;

  /// Number of bits per block. Note this is a tradeoff between space efficiency and perf.
  static constexpr size_type words_per_block = 4;
  /// Number of bits in a word
  static constexpr size_type bits_per_word = sizeof(word_type) * CHAR_BIT;
  /// Number of bits in a block
  static constexpr size_type bits_per_block = words_per_block * bits_per_word;

  /**
   * @brief Constructs an empty bitset
   *
   * @param allocator Allocator used for allocating device storage
   */
  constexpr dynamic_bitset(Allocator const& allocator = Allocator{});

  /**
   * @brief Appends the given element `value` to the end of the bitset
   *
   * This API may involve data reallocation if the current storage is exhausted.
   *
   * @param value Boolean value of the new bit to be added
   */
  constexpr void push_back(bool value) noexcept;

  /**
   * @brief Sets the target bit indexed by `index` to a specified `value`.
   *
   * @param index Position of bit to be modified
   * @param value New value of the target bit
   */
  constexpr void set(size_type index, bool value) noexcept;

  /**
   * @brief Sets the last bit to a specified value
   *
   * @param value New value of the last bit
   */
  constexpr void set_last(bool value) noexcept;

  /**
   * @brief For any element `keys_begin[i]` in the range `[keys_begin, keys_end)`, stores the
   * boolean value at position `keys_begin[i]` to `output_begin[i]`.
   *
   * @tparam KeyIt Device-accessible iterator whose `value_type` can be converted to bitset's
   * `size_type`
   * @tparam OutputIt Device-accessible iterator whose `value_type` can be constructed from boolean
   * type
   *
   * @param keys_begin Begin iterator to keys list whose values are queried
   * @param keys_end End iterator to keys list
   * @param outputs_begin Begin iterator to outputs of test operation
   * @param stream Stream to execute test kernel
   */
  template <typename KeyIt, typename OutputIt>
  constexpr void test(KeyIt keys_begin,
                      KeyIt keys_end,
                      OutputIt outputs_begin,
                      cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief For any element `keys_begin[i]` in the range `[keys_begin, keys_end)`, stores total
   * count of `1` bits preceeding (but not including) position `keys_begin[i]` to `output_begin[i]`.
   *
   * @tparam KeyIt Device-accessible iterator whose `value_type` can be converted to bitset's
   * `size_type`
   * @tparam OutputIt Device-accessible iterator whose `value_type` can be constructed from bitset's
   * `size_type`
   *
   * @param keys_begin Begin iterator to keys list whose ranks are queried
   * @param keys_end End iterator to keys list
   * @param outputs_begin Begin iterator to outputs ranks list
   * @param stream Stream to execute ranks kernel
   */
  template <typename KeyIt, typename OutputIt>
  constexpr void rank(KeyIt keys_begin,
                      KeyIt keys_end,
                      OutputIt outputs_begin,
                      cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief For any element `keys_begin[i]` in the range `[keys_begin, keys_end)`, stores the
   * position of `keys_begin[i]`th `1` bit to `output_begin[i]`.
   *
   * @tparam KeyIt Device-accessible iterator whose `value_type` can be converted to bitset's
   * `size_type`
   * @tparam OutputIt Device-accessible iterator whose `value_type` can be constructed from bitset's
   * `size_type`
   *
   * @param keys_begin Begin iterator to keys list whose select values are queried
   * @param keys_end End iterator to keys list
   * @param outputs_begin Begin iterator to outputs selects list
   * @param stream Stream to execute selects kernel
   */
  template <typename KeyIt, typename OutputIt>
  constexpr void select(KeyIt keys_begin,
                        KeyIt keys_end,
                        OutputIt outputs_begin,
                        cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) noexcept;

  using rank_type = cuco::experimental::detail::rank;  ///< Rank type

  /**
   *@brief Struct to hold all storage refs needed by reference
   */
  // TODO: this is not a real ref type, to be changed
  struct storage_ref_type {
    const word_type* words_ref_;  ///< Words ref

    const rank_type* ranks_true_ref_;    ///< Ranks ref for 1 bits
    const size_type* selects_true_ref_;  ///< Selects ref for 1 bits

    const rank_type* ranks_false_ref_;    ///< Ranks ref for 0 bits
    const size_type* selects_false_ref_;  ///< Selects ref 0 bits
  };

  /**
   * @brief Device non-owning reference type of dynamic_bitset
   */
  class reference {
   public:
    /**
     * @brief Constructs a reference
     *
     * @param storage Struct with non-owning refs to bitset storage arrays
     */
    __host__ __device__ explicit constexpr reference(storage_ref_type storage) noexcept;

    /**
     * @brief Access value of a single bit
     *
     * @param key Position of bit
     *
     * @return Value of bit at position specified by key
     */
    [[nodiscard]] __device__ constexpr bool test(size_type key) const noexcept;

    /**
     * @brief Access a single word of internal storage
     *
     * @param word_id Index of word
     *
     * @return Word at position specified by index
     */
    [[nodiscard]] __device__ constexpr word_type word(size_type word_id) const noexcept;

    /**
     * @brief Find position of first set bit starting from a given position (inclusive)
     *
     * @param key Position of starting bit
     *
     * @return Index of next set bit
     */
    [[nodiscard]] __device__ size_type find_next(size_type key) const noexcept;

    /**
     * @brief Find number of set bits (rank) in all positions before the input position (exclusive)
     *
     * @param key Input bit position
     *
     * @return Rank of input position
     */
    [[nodiscard]] __device__ constexpr size_type rank(size_type key) const noexcept;

    /**
     * @brief Find position of Nth set (1) bit counting from start
     *
     * @param count Input N
     *
     * @return Position of Nth set bit
     */
    [[nodiscard]] __device__ constexpr size_type select(size_type count) const noexcept;

    /**
     * @brief Find position of Nth not-set (0) bit counting from start
     *
     * @param count Input N
     *
     * @return Position of Nth not-set bit
     */
    [[nodiscard]] __device__ constexpr size_type select_false(size_type count) const noexcept;

   private:
    /**
     * @brief Helper function for select operation that computes an initial rank estimate
     *
     * @param count Input count for which select operation is being performed
     * @param selects Selects array
     * @param ranks Ranks array
     *
     * @return index in ranks which corresponds to highest rank less than count (least upper bound)
     */
    template <typename SelectsRef, typename RanksRef>
    [[nodiscard]] __device__ constexpr size_type initial_rank_estimate(
      size_type count, const SelectsRef& selects, const RanksRef& ranks) const noexcept;

    /**
     * @brief Subtract rank estimate from input count and return an increment to word_id
     *
     * @tparam Rank type
     *
     * @param count Input count that will be updated
     * @param rank  Initial rank estimate for count
     *
     * @return Increment to word_id based on rank values
     */
    template <typename Rank>
    [[nodiscard]] __device__ constexpr size_type subtract_rank_from_count(size_type& count,
                                                                          Rank rank) const noexcept;

    /**
     * @brief Find position of Nth set bit in a 64-bit word
     *
     * @param N Input count
     *
     * @return Position of Nth set bit
     */
    [[nodiscard]] __device__ size_type select_bit_in_word(size_type N,
                                                          word_type word) const noexcept;

    storage_ref_type storage_;  ///< Non-owning storage
  };

  using ref_type = reference;  ///< Non-owning container ref type

  /**
   * @brief Gets non-owning device ref of the current object
   *
   * @return Device ref of the current `dynamic_bitset` object
   */
  [[nodiscard]] constexpr ref_type ref() const noexcept;

  /**
   * @brief Gets the number of bits dynamic_bitset holds
   *
   * @return Number of bits dynamic_bitset holds
   */
  [[nodiscard]] constexpr size_type size() const noexcept;

 private:
  /// Type of the allocator to (de)allocate ranks
  using rank_allocator_type =
    typename std::allocator_traits<Allocator>::template rebind_alloc<rank_type>;
  /// Type of the allocator to (de)allocate indices
  using size_allocator_type =
    typename std::allocator_traits<Allocator>::template rebind_alloc<size_type>;

  allocator_type allocator_;  ///< Words allocator
  size_type n_bits_;          ///< Number of bits dynamic_bitset currently holds
  bool is_built_;  ///< Flag indicating whether the rank and select indices are built or not

  /// Words vector that represents all bits
  thrust::device_vector<word_type, allocator_type> words_;
  /// Rank values for every 256-th bit (4-th word)
  thrust::device_vector<rank_type, rank_allocator_type> ranks_true_;
  /// Same as ranks_ but for `0` bits
  thrust::device_vector<rank_type, rank_allocator_type> ranks_false_;
  /// Block indices of (0, 256, 512...)th `1` bit
  thrust::device_vector<size_type, size_allocator_type> selects_true_;
  /// Same as selects_, but for `0` bits
  thrust::device_vector<size_type, size_allocator_type> selects_false_;

  /**
   * @brief Builds indexes for rank and select
   *
   * @param stream Stream to execute kernels
   */
  constexpr void build(cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief Populates rank and select indexes for true or false bits
   *
   * @param ranks Output array of ranks
   * @param selects Output array of selects
   * @param flip_bits If true, negate bits to construct indexes for false bits
   * @param stream Stream to execute kernels
   */
  constexpr void build_ranks_and_selects(
    thrust::device_vector<rank_type, rank_allocator_type>& ranks,
    thrust::device_vector<size_type, size_allocator_type>& selects,
    bool flip_bits,
    cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});
};

}  // namespace detail
}  // namespace experimental
}  // namespace cuco

#include <cuco/detail/trie/dynamic_bitset/dynamic_bitset.inl>
