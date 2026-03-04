/*
 * Copyright (c) 2024, NVIDIA CORPORATION.
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

#include <cstdint>

namespace cuco {

template <class Hash, class Word, uint32_t WordsPerBlock>
__host__
  __device__ constexpr default_filter_policy<Hash, Word, WordsPerBlock>::default_filter_policy(
    uint32_t pattern_bits, Hash hash)
  : impl_{pattern_bits, hash}
{
}

template <class Hash, class Word, uint32_t WordsPerBlock>
__device__ constexpr typename default_filter_policy<Hash, Word, WordsPerBlock>::hash_result_type
default_filter_policy<Hash, Word, WordsPerBlock>::hash(
  typename default_filter_policy<Hash, Word, WordsPerBlock>::hash_argument_type const& key) const
{
  return impl_.hash(key);
}

template <class Hash, class Word, uint32_t WordsPerBlock>
template <class Extent>
__device__ constexpr auto default_filter_policy<Hash, Word, WordsPerBlock>::block_index(
  typename default_filter_policy<Hash, Word, WordsPerBlock>::hash_result_type hash,
  Extent num_blocks) const
{
  return impl_.block_index(hash, num_blocks);
}

template <class Hash, class Word, uint32_t WordsPerBlock>
__device__ constexpr typename default_filter_policy<Hash, Word, WordsPerBlock>::word_type
default_filter_policy<Hash, Word, WordsPerBlock>::word_pattern(
  default_filter_policy<Hash, Word, WordsPerBlock>::hash_result_type hash,
  std::uint32_t word_index) const
{
  return impl_.word_pattern(hash, word_index);
}

}  // namespace cuco