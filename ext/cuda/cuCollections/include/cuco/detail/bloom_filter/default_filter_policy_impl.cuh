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

#include <cuco/detail/error.hpp>

#include <cuda/std/bit>
#include <cuda/std/limits>
#include <cuda/std/tuple>
#include <cuda/std/type_traits>

#include <cstdint>
#include <nv/target>

namespace cuco::detail {

template <class Hash, class Word, uint32_t WordsPerBlock>
class default_filter_policy_impl {
 public:
  using hasher             = Hash;
  using word_type          = Word;
  using hash_argument_type = typename hasher::argument_type;
  using hash_result_type   = decltype(std::declval<hasher>()(std::declval<hash_argument_type>()));

  static constexpr std::uint32_t words_per_block = WordsPerBlock;

 private:
  static constexpr std::uint32_t word_bits       = cuda::std::numeric_limits<word_type>::digits;
  static constexpr std::uint32_t bit_index_width = cuda::std::bit_width(word_bits - 1);

 public:
  __host__ __device__ explicit constexpr default_filter_policy_impl(uint32_t pattern_bits,
                                                                    Hash hash)
    : pattern_bits_{pattern_bits},
      min_bits_per_word_{pattern_bits_ / words_per_block},
      remainder_bits_{pattern_bits_ % words_per_block},
      hash_{hash}
  {
    NV_DISPATCH_TARGET(
      NV_IS_HOST,
      (  // This ensures each word in the block has at least one bit set; otherwise we would never
         // use some of the words
        constexpr uint32_t min_pattern_bits = words_per_block;

        // The maximum number of bits to be set for a key is capped by the total number of bits in
        // the filter block
        constexpr uint32_t max_pattern_bits = word_bits * words_per_block;

        constexpr uint32_t hash_bits = cuda::std::numeric_limits<hash_result_type>::digits;
        constexpr uint32_t max_pattern_bits_from_hash = hash_bits / bit_index_width;
        CUCO_EXPECTS(
          pattern_bits <= max_pattern_bits_from_hash,
          "`hash_result_type` too narrow to generate the requested number of `pattern_bits`");
        CUCO_EXPECTS(pattern_bits_ >= min_pattern_bits,
                     "`pattern_bits` must be at least `words_per_block`");
        CUCO_EXPECTS(pattern_bits_ <= max_pattern_bits,
                     "`pattern_bits` must be less than the total number of bits in a filter "
                     "block");))
    // TODO find a proper way to perform input checks/assertions on device without destroying the
    // context (e.g. __trap())
  }

  __device__ constexpr hash_result_type hash(hash_argument_type const& key) const
  {
    return hash_(key);
  }

  template <class Extent>
  __device__ constexpr auto block_index(hash_result_type hash, Extent num_blocks) const
  {
    return hash % num_blocks;
  }

  __device__ constexpr word_type word_pattern(hash_result_type hash, std::uint32_t word_index) const
  {
    word_type constexpr bit_index_mask = (word_type{1} << bit_index_width) - 1;

    auto const bits_so_far = min_bits_per_word_ * word_index +
                             (word_index < remainder_bits_ ? word_index : remainder_bits_);

    hash >>= bits_so_far * bit_index_width;

    word_type word              = 0;
    int32_t const bits_per_word = min_bits_per_word_ + (word_index < remainder_bits_ ? 1 : 0);

    for (int32_t bit = 0; bit < bits_per_word; ++bit) {
      word |= word_type{1} << (hash & bit_index_mask);
      hash >>= bit_index_width;
    }

    return word;
  }

 private:
  uint32_t pattern_bits_;
  uint32_t min_bits_per_word_;
  uint32_t remainder_bits_;
  hasher hash_;
};

}  // namespace cuco::detail