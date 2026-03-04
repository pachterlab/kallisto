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

#include <cuco/detail/trie/dynamic_bitset/kernels.cuh>
#include <cuco/detail/utility/cuda.hpp>
#include <cuco/detail/utils.hpp>

#include <cub/device/device_scan.cuh>
#include <cub/device/device_select.cuh>
#include <cuda/std/bit>
#include <thrust/device_vector.h>
#include <thrust/iterator/discard_iterator.h>

namespace cuco {
namespace experimental {
namespace detail {

template <class Allocator>
constexpr dynamic_bitset<Allocator>::dynamic_bitset(Allocator const& allocator)
  : allocator_{allocator},
    n_bits_{0},
    is_built_{false},
    words_{allocator},
    ranks_true_{allocator},
    ranks_false_{allocator},
    selects_true_{allocator},
    selects_false_{allocator}
{
}

template <class Allocator>
constexpr void dynamic_bitset<Allocator>::push_back(bool bit) noexcept
{
  if (n_bits_ % bits_per_block == 0) {
    words_.resize(words_.size() + words_per_block);  // Extend storage by one block
  }

  set(n_bits_++, bit);
}

template <class Allocator>
constexpr void dynamic_bitset<Allocator>::set(size_type index, bool bit) noexcept
{
  is_built_         = false;
  size_type word_id = index / bits_per_word;
  size_type bit_id  = index % bits_per_word;
  if (bit) {
    words_[word_id] |= 1UL << bit_id;
  } else {
    words_[word_id] &= ~(1UL << bit_id);
  }
}

template <class Allocator>
constexpr void dynamic_bitset<Allocator>::set_last(bool bit) noexcept
{
  set(n_bits_ - 1, bit);
}

template <class Allocator>
template <typename KeyIt, typename OutputIt>
constexpr void dynamic_bitset<Allocator>::test(KeyIt keys_begin,
                                               KeyIt keys_end,
                                               OutputIt outputs_begin,
                                               cuda::stream_ref stream) noexcept

{
  build();
  auto const num_keys = cuco::detail::distance(keys_begin, keys_end);
  if (num_keys == 0) { return; }

  auto const grid_size = cuco::detail::grid_size(num_keys);

  bitset_test_kernel<<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
    ref(), keys_begin, outputs_begin, num_keys);
}

template <class Allocator>
template <typename KeyIt, typename OutputIt>
constexpr void dynamic_bitset<Allocator>::rank(KeyIt keys_begin,
                                               KeyIt keys_end,
                                               OutputIt outputs_begin,
                                               cuda::stream_ref stream) noexcept
{
  build();
  auto const num_keys = cuco::detail::distance(keys_begin, keys_end);
  if (num_keys == 0) { return; }

  auto const grid_size = cuco::detail::grid_size(num_keys);

  bitset_rank_kernel<<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
    ref(), keys_begin, outputs_begin, num_keys);
}

template <class Allocator>
template <typename KeyIt, typename OutputIt>
constexpr void dynamic_bitset<Allocator>::select(KeyIt keys_begin,
                                                 KeyIt keys_end,
                                                 OutputIt outputs_begin,
                                                 cuda::stream_ref stream) noexcept

{
  build();
  auto const num_keys = cuco::detail::distance(keys_begin, keys_end);
  if (num_keys == 0) { return; }

  auto const grid_size = cuco::detail::grid_size(num_keys);

  bitset_select_kernel<<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
    ref(), keys_begin, outputs_begin, num_keys);
}

template <class Allocator>
constexpr void dynamic_bitset<Allocator>::build_ranks_and_selects(
  thrust::device_vector<rank_type, rank_allocator_type>& ranks,
  thrust::device_vector<size_type, size_allocator_type>& selects,
  bool flip_bits,
  cuda::stream_ref stream)
{
  if (n_bits_ == 0) { return; }

  // Step 1. Compute prefix sum of per-word bit counts
  // Population counts for each word
  size_type const num_words = words_.size();
  // Sized to have one extra entry for subsequent prefix sum
  auto const bit_counts_size = num_words + 1;

  thrust::device_vector<size_type, size_allocator_type> bit_counts(num_words + 1, this->allocator_);
  auto const bit_counts_begin = thrust::raw_pointer_cast(bit_counts.data());

  auto grid_size = cuco::detail::grid_size(num_words);
  bit_counts_kernel<<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
    thrust::raw_pointer_cast(words_.data()), bit_counts_begin, num_words, flip_bits);

  std::size_t temp_storage_bytes = 0;
  using temp_allocator_type =
    typename std::allocator_traits<allocator_type>::template rebind_alloc<char>;
  auto temp_allocator = temp_allocator_type{this->allocator_};

  CUCO_CUDA_TRY(cub::DeviceScan::ExclusiveSum(nullptr,
                                              temp_storage_bytes,
                                              bit_counts_begin,
                                              bit_counts_begin,
                                              bit_counts_size,
                                              stream.get()));

  // Allocate temporary storage
  auto d_temp_storage = temp_allocator.allocate(temp_storage_bytes);

  CUCO_CUDA_TRY(cub::DeviceScan::ExclusiveSum(thrust::raw_pointer_cast(d_temp_storage),
                                              temp_storage_bytes,
                                              bit_counts_begin,
                                              bit_counts_begin,
                                              bit_counts_size,
                                              stream.get()));

  temp_allocator.deallocate(d_temp_storage, temp_storage_bytes);

  // Step 2. Compute ranks
  auto const num_blocks = (num_words - 1) / words_per_block + 2;
  ranks.resize(num_blocks);

  grid_size = cuco::detail::grid_size(num_blocks);
  encode_ranks_from_prefix_bit_counts<<<grid_size,
                                        cuco::detail::default_block_size(),
                                        0,
                                        stream.get()>>>(bit_counts_begin,
                                                        thrust::raw_pointer_cast(ranks.data()),
                                                        num_words,
                                                        num_blocks,
                                                        words_per_block);

  // Step 3. Compute selects
  thrust::device_vector<size_type, size_allocator_type> select_markers(num_blocks,
                                                                       this->allocator_);
  auto const select_markers_begin = thrust::raw_pointer_cast(select_markers.data());

  mark_blocks_with_select_entries<<<grid_size,
                                    cuco::detail::default_block_size(),
                                    0,
                                    stream.get()>>>(
    bit_counts_begin, select_markers_begin, num_blocks, words_per_block, bits_per_block);

  auto d_sum = reinterpret_cast<size_type*>(
    thrust::raw_pointer_cast(temp_allocator.allocate(sizeof(size_type))));
  CUCO_CUDA_TRY(cub::DeviceReduce::Sum(
    nullptr, temp_storage_bytes, select_markers_begin, d_sum, num_blocks, stream.get()));

  d_temp_storage = temp_allocator.allocate(temp_storage_bytes);

  CUCO_CUDA_TRY(cub::DeviceReduce::Sum(thrust::raw_pointer_cast(d_temp_storage),
                                       temp_storage_bytes,
                                       select_markers_begin,
                                       d_sum,
                                       num_blocks,
                                       stream.get()));

  size_type num_selects{};
  CUCO_CUDA_TRY(
    cudaMemcpyAsync(&num_selects, d_sum, sizeof(size_type), cudaMemcpyDeviceToHost, stream.get()));
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
  stream.sync();
#else
  stream.wait();
#endif
  temp_allocator.deallocate(thrust::device_ptr<char>{reinterpret_cast<char*>(d_sum)},
                            sizeof(size_type));
  temp_allocator.deallocate(d_temp_storage, temp_storage_bytes);

  selects.resize(num_selects);

  auto const select_begin = thrust::raw_pointer_cast(selects.data());

  CUCO_CUDA_TRY(cub::DeviceSelect::Flagged(nullptr,
                                           temp_storage_bytes,
                                           thrust::make_counting_iterator(0UL),
                                           select_markers_begin,
                                           select_begin,
                                           thrust::make_discard_iterator(),
                                           num_blocks,
                                           stream.get()));

  d_temp_storage = temp_allocator.allocate(temp_storage_bytes);

  CUCO_CUDA_TRY(cub::DeviceSelect::Flagged(thrust::raw_pointer_cast(d_temp_storage),
                                           temp_storage_bytes,
                                           thrust::make_counting_iterator(0UL),
                                           select_markers_begin,
                                           select_begin,
                                           thrust::make_discard_iterator(),
                                           num_blocks,
                                           stream.get()));

  temp_allocator.deallocate(d_temp_storage, temp_storage_bytes);
}

template <class Allocator>
constexpr void dynamic_bitset<Allocator>::build(cuda::stream_ref stream) noexcept
{
  if (not is_built_) {
    build_ranks_and_selects(ranks_true_, selects_true_, false, stream);   // 1 bits
    build_ranks_and_selects(ranks_false_, selects_false_, true, stream);  // 0 bits
    is_built_ = true;
  }
}

template <class Allocator>
constexpr dynamic_bitset<Allocator>::ref_type dynamic_bitset<Allocator>::ref() const noexcept
{
  return ref_type{storage_ref_type{thrust::raw_pointer_cast(words_.data()),
                                   thrust::raw_pointer_cast(ranks_true_.data()),
                                   thrust::raw_pointer_cast(selects_true_.data()),
                                   thrust::raw_pointer_cast(ranks_false_.data()),
                                   thrust::raw_pointer_cast(selects_false_.data())}};
}

template <class Allocator>
constexpr dynamic_bitset<Allocator>::size_type dynamic_bitset<Allocator>::size() const noexcept
{
  return n_bits_;
}

// Device reference implementations

template <class Allocator>
__host__ __device__ constexpr dynamic_bitset<Allocator>::reference::reference(
  storage_ref_type storage) noexcept
  : storage_{storage}
{
}

template <class Allocator>
__device__ constexpr bool dynamic_bitset<Allocator>::reference::test(size_type key) const noexcept
{
  return (storage_.words_ref_[key / bits_per_word] >> (key % bits_per_word)) & 1UL;
}

template <class Allocator>
__device__ constexpr typename dynamic_bitset<Allocator>::word_type
dynamic_bitset<Allocator>::reference::word(size_type word_id) const noexcept
{
  return storage_.words_ref_[word_id];
}

template <class Allocator>
__device__ typename dynamic_bitset<Allocator>::size_type
dynamic_bitset<Allocator>::reference::find_next(size_type key) const noexcept
{
  size_type word_id = key / bits_per_word;
  size_type bit_id  = key % bits_per_word;
  word_type word    = storage_.words_ref_[word_id];
  word &= ~(0UL) << bit_id;
  while (word == 0) {
    word = storage_.words_ref_[++word_id];
  }
  return word_id * bits_per_word + __ffsll(word) - 1;  // cuda intrinsic
}

template <class Allocator>
__device__ constexpr typename dynamic_bitset<Allocator>::size_type
dynamic_bitset<Allocator>::reference::rank(size_type key) const noexcept
{
  size_type word_id   = key / bits_per_word;
  size_type bit_id    = key % bits_per_word;
  size_type rank_id   = word_id / words_per_block;
  size_type offset_id = word_id % words_per_block;

  auto rank   = storage_.ranks_true_ref_[rank_id];
  size_type n = rank.base();

  if (offset_id != 0) { n += rank.offsets_[offset_id - 1]; }

  n += cuda::std::popcount(storage_.words_ref_[word_id] & ((1UL << bit_id) - 1));

  return n;
}

template <class Allocator>
__device__ constexpr typename dynamic_bitset<Allocator>::size_type
dynamic_bitset<Allocator>::reference::select(size_type count) const noexcept
{
  auto rank_id = initial_rank_estimate(count, storage_.selects_true_ref_, storage_.ranks_true_ref_);
  auto rank    = storage_.ranks_true_ref_[rank_id];

  size_type word_id = rank_id * words_per_block;
  word_id += subtract_rank_from_count(count, rank);

  return word_id * bits_per_word + select_bit_in_word(count, storage_.words_ref_[word_id]);
}

template <class Allocator>
__device__ constexpr typename dynamic_bitset<Allocator>::size_type
dynamic_bitset<Allocator>::reference::select_false(size_type count) const noexcept
{
  auto rank_id =
    initial_rank_estimate(count, storage_.selects_false_ref_, storage_.ranks_false_ref_);
  auto rank = storage_.ranks_false_ref_[rank_id];

  size_type word_id = rank_id * words_per_block;
  word_id += subtract_rank_from_count(count, rank);

  return word_id * bits_per_word + select_bit_in_word(count, ~(storage_.words_ref_[word_id]));
}

template <class Allocator>
template <typename SelectsRef, typename RanksRef>
__device__ constexpr typename dynamic_bitset<Allocator>::size_type
dynamic_bitset<Allocator>::reference::initial_rank_estimate(size_type count,
                                                            SelectsRef const& selects,
                                                            RanksRef const& ranks) const noexcept
{
  size_type block_id = count / (bits_per_word * words_per_block);
  size_type begin    = selects[block_id];
  size_type end      = selects[block_id + 1] + 1UL;

  if (begin + 10 >= end) {  // Linear search
    while (count >= ranks[begin + 1].base()) {
      ++begin;
    }
  } else {  // Binary search
    while (begin + 1 < end) {
      size_type middle = (begin + end) / 2;
      if (count < ranks[middle].base()) {
        end = middle;
      } else {
        begin = middle;
      }
    }
  }
  return begin;
}

template <class Allocator>
template <typename Rank>
__device__ constexpr typename dynamic_bitset<Allocator>::size_type
dynamic_bitset<Allocator>::reference::subtract_rank_from_count(size_type& count,
                                                               Rank rank) const noexcept
{
  count -= rank.base();

  bool a0       = count >= rank.offsets_[0];
  bool a1       = count >= rank.offsets_[1];
  bool a2       = count >= rank.offsets_[2];
  size_type inc = a0 + a1 + a2;

  count -= (inc > 0) * rank.offsets_[inc - (inc > 0)];

  return inc;
}

template <class Allocator>
__device__ typename dynamic_bitset<Allocator>::size_type
dynamic_bitset<Allocator>::reference::select_bit_in_word(size_type N, word_type word) const noexcept
{
  for (size_type pos = 0; pos < N; pos++) {
    word &= word - 1;
  }
  return __ffsll(word & -word) - 1;  // cuda intrinsic
}
}  // namespace detail
}  // namespace experimental
}  // namespace cuco
