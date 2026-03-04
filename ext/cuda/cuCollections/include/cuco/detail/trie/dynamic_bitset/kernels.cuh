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

#include <cuco/detail/trie/dynamic_bitset/dynamic_bitset.cuh>
#include <cuco/detail/utility/cuda.cuh>
#include <cuco/detail/utility/cuda.hpp>

#include <cuda/std/bit>

namespace cuco {
namespace experimental {
namespace detail {

CUCO_SUPPRESS_KERNEL_WARNINGS
/*
 * @brief Test bits for a range of keys
 *
 * @tparam BitsetRef Bitset reference type
 * @tparam KeyIt Device-accessible iterator whose `value_type` can be converted to bitset's
 * `size_type`
 * @tparam OutputIt Device-accessible iterator whose `value_type` can be constructed from boolean
 * type
 *
 * @param ref Bitset ref
 * @param keys Begin iterator to keys
 * @param outputs Begin iterator to outputs
 * @param num_keys Number of input keys
 */
template <typename BitsetRef, typename KeyIt, typename OutputIt>
CUCO_KERNEL void bitset_test_kernel(BitsetRef ref,
                                    KeyIt keys,
                                    OutputIt outputs,
                                    cuco::detail::index_type num_keys)
{
  auto key_id       = cuco::detail::global_thread_id();
  auto const stride = cuco::detail::grid_stride();

  while (key_id < num_keys) {
    outputs[key_id] = ref.test(keys[key_id]);
    key_id += stride;
  }
}

/*
 * @brief Gather rank values for a range of keys
 *
 * @tparam BitsetRef Bitset reference type
 * @tparam KeyIt Device-accessible iterator whose `value_type` can be converted to bitset's
 * `size_type`
 * @tparam OutputIt Device-accessible iterator whose `value_type` can be constructed from bitset's
 * `size_type`
 *
 * @param ref Bitset ref
 * @param keys Begin iterator to keys
 * @param outputs Begin iterator to outputs
 * @param num_keys Number of input keys
 */
template <typename BitsetRef, typename KeyIt, typename OutputIt>
CUCO_KERNEL void bitset_rank_kernel(BitsetRef ref,
                                    KeyIt keys,
                                    OutputIt outputs,
                                    cuco::detail::index_type num_keys)
{
  auto key_id       = cuco::detail::global_thread_id();
  auto const stride = cuco::detail::grid_stride();

  while (key_id < num_keys) {
    outputs[key_id] = ref.rank(keys[key_id]);
    key_id += stride;
  }
}

/*
 * @brief Gather select values for a range of keys
 *
 * @tparam BitsetRef Bitset reference type
 * @tparam KeyIt Device-accessible iterator whose `value_type` can be converted to bitset's
 * `size_type`
 * @tparam OutputIt Device-accessible iterator whose `value_type` can be constructed from bitset's
 * `size_type`
 *
 * @param ref Bitset ref
 * @param keys Begin iterator to keys
 * @param outputs Begin iterator to outputs
 * @param num_keys Number of input keys
 */
template <typename BitsetRef, typename KeyIt, typename OutputIt>
CUCO_KERNEL void bitset_select_kernel(BitsetRef ref,
                                      KeyIt keys,
                                      OutputIt outputs,
                                      cuco::detail::index_type num_keys)
{
  auto key_id       = cuco::detail::global_thread_id();
  auto const stride = cuco::detail::grid_stride();

  while (key_id < num_keys) {
    outputs[key_id] = ref.select(keys[key_id]);
    key_id += stride;
  }
}

/*
 * @brief Computes number of set or not-set bits in each word
 *
 * @tparam WordType Word type
 * @tparam SizeType Size type
 *
 * @param words Input array of words
 * @param bit_counts Output array of per-word bit counts
 * @param num_words Number of words
 * @param flip_bits Boolean to request negation of words before counting bits
 */
template <typename WordType, typename SizeType>
CUCO_KERNEL void bit_counts_kernel(WordType const* words,
                                   SizeType* bit_counts,
                                   cuco::detail::index_type num_words,
                                   bool flip_bits)
{
  auto word_id      = cuco::detail::global_thread_id();
  auto const stride = cuco::detail::grid_stride();

  while (word_id < num_words) {
    auto word           = words[word_id];
    bit_counts[word_id] = cuda::std::popcount(flip_bits ? ~word : word);
    word_id += stride;
  }
}

/*
 * @brief Compute rank values at block size intervals.
 *
 * ranks[i] = Number of set bits in [0, i) range
 * This kernel transforms prefix sum array of per-word bit counts
 * into base-delta encoding style of `rank` struct.
 * Since prefix sum is available, there are no dependencies across blocks.

 * @tparam SizeType Size type
 *
 * @param prefix_bit_counts Prefix sum array of per-word bit counts
 * @param ranks Output array of ranks
 * @param num_words Length of input array
 * @param num_blocks Length of ouput array
 * @param words_per_block Number of words in each block
 */
template <typename SizeType>
CUCO_KERNEL void encode_ranks_from_prefix_bit_counts(const SizeType* prefix_bit_counts,
                                                     rank* ranks,
                                                     SizeType num_words,
                                                     SizeType num_blocks,
                                                     SizeType words_per_block)
{
  auto rank_id      = cuco::detail::global_thread_id();
  auto const stride = cuco::detail::grid_stride();

  while (rank_id < num_blocks) {
    SizeType word_id = rank_id * words_per_block;

    // Set base value of rank
    auto& rank = ranks[rank_id];
    rank.set_base(prefix_bit_counts[word_id]);

    if (rank_id < num_blocks - 1) {
      // For each subsequent word in this block, compute deltas from base
      for (SizeType block_offset = 0; block_offset < words_per_block - 1; block_offset++) {
        auto delta = prefix_bit_counts[word_id + block_offset + 1] - prefix_bit_counts[word_id];
        rank.offsets_[block_offset] = delta;
      }
    }
    rank_id += stride;
  }
}

/*
 * @brief Compute select values at block size intervals.
 *
 * selects[i] = Position of (i+ 1)th set bit
 * This kernel check for blocks where prefix sum crosses a multiple of `bits_per_block`.
 * Such blocks are marked in the output boolean array
 *
 * @tparam SizeType Size type
 *
 * @param prefix_bit_counts Prefix sum array of per-word bit counts
 * @param selects_markers Ouput array indicating whether a block has selects entry or not
 * @param num_blocks Length of ouput array
 * @param words_per_block Number of words in each block
 * @param bits_per_block Number of bits in each block
 */
template <typename SizeType>
CUCO_KERNEL void mark_blocks_with_select_entries(SizeType const* prefix_bit_counts,
                                                 SizeType* select_markers,
                                                 SizeType num_blocks,
                                                 SizeType words_per_block,
                                                 SizeType bits_per_block)
{
  auto block_id     = cuco::detail::global_thread_id();
  auto const stride = cuco::detail::grid_stride();

  while (block_id < num_blocks) {
    if (block_id == 0) {  // Block 0 always has a selects entry
      select_markers[block_id] = 1;
      block_id += stride;
      continue;
    }

    select_markers[block_id] = 0;  // Always clear marker first
    SizeType word_id         = block_id * words_per_block;
    SizeType prev_count      = prefix_bit_counts[word_id];

    for (size_t block_offset = 1; block_offset <= words_per_block; block_offset++) {
      SizeType count = prefix_bit_counts[word_id + block_offset];

      // Selects entry is added when cumulative bitcount crosses a multiple of bits_per_block
      if ((prev_count - 1) / bits_per_block != (count - 1) / bits_per_block) {
        select_markers[block_id] = 1;
        break;
      }
      prev_count = count;
    }

    block_id += stride;
  }
}

}  // namespace detail
}  // namespace experimental
}  // namespace cuco
