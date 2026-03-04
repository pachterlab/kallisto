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

#include <cuco/detail/bitwise_compare.cuh>
#include <cuco/detail/static_map/kernels.cuh>
#include <cuco/detail/utility/cuda.hpp>
#include <cuco/detail/utils.hpp>
#include <cuco/operator.hpp>
#include <cuco/static_map_ref.cuh>

#include <cuda/stream_ref>

#include <algorithm>
#include <cstddef>

namespace cuco {
namespace experimental {

template <typename Key,
          typename T,
          typename Extent,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename Allocator,
          typename Storage>
constexpr dynamic_map<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::
  dynamic_map(Extent initial_capacity,
              empty_key<Key> empty_key_sentinel,
              empty_value<T> empty_value_sentinel,
              KeyEqual const& pred,
              ProbingScheme const& probing_scheme,
              cuda_thread_scope<Scope> scope,
              Storage storage,
              Allocator const& alloc,
              cuda::stream_ref stream)
  : size_{0},
    capacity_{initial_capacity},
    min_insert_size_{static_cast<size_type>(1E4)},
    max_load_factor_{0.60},
    alloc_{alloc}
{
  submaps_.push_back(
    std::make_unique<
      cuco::static_map<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>>(
      initial_capacity,
      empty_key_sentinel,
      empty_value_sentinel,
      pred,
      probing_scheme,
      scope,
      storage,
      alloc,
      stream));
}

template <typename Key,
          typename T,
          typename Extent,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename Allocator,
          typename Storage>
template <typename InputIt>
void dynamic_map<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::insert(
  InputIt first, InputIt last, cuda::stream_ref stream)
{
  auto num_to_insert = cuco::detail::distance(first, last);
  this->reserve(size_ + num_to_insert, stream);

  std::size_t submap_idx = 0;
  while (num_to_insert > 0) {
    auto& cur = submaps_[submap_idx];

    auto capacity_remaining = max_load_factor_ * cur->capacity() - cur->size();
    // If we are tying to insert some of the remaining keys into this submap, we can insert
    // only if we meet the minimum insert size.
    if (capacity_remaining >= min_insert_size_) {
      auto const n = std::min(static_cast<detail::index_type>(capacity_remaining), num_to_insert);

      std::size_t h_num_successes = cur->insert(first, first + n, stream);

      size_ += h_num_successes;
      first += n;
      num_to_insert -= n;
    }
    submap_idx++;
  }
}

template <typename Key,
          typename T,
          typename Extent,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename Allocator,
          typename Storage>
void dynamic_map<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::reserve(
  size_type n, cuda::stream_ref stream)
{
  size_type num_elements_remaining = n;
  std::size_t submap_idx           = 0;
  while (num_elements_remaining > 0) {
    std::size_t submap_capacity;

    // if the submap already exists
    if (submap_idx < submaps_.size()) {
      submap_capacity = submaps_[submap_idx]->capacity();
    }
    // if the submap does not exist yet, create it
    else {
      empty_key<Key> empty_key_sentinel{submaps_.front()->empty_key_sentinel()};
      empty_value<T> empty_value_sentinel{submaps_.front()->empty_value_sentinel()};

      submap_capacity = capacity_;
      submaps_.push_back(std::make_unique<map_type>(submap_capacity,
                                                    empty_key_sentinel,
                                                    empty_value_sentinel,
                                                    KeyEqual{},
                                                    ProbingScheme{},
                                                    cuda_thread_scope<Scope>{},
                                                    Storage{},
                                                    alloc_,
                                                    stream));
      capacity_ *= 2;
    }

    num_elements_remaining -= max_load_factor_ * submap_capacity - min_insert_size_;
    submap_idx++;
  }
}

template <typename Key,
          typename T,
          typename Extent,
          cuda::thread_scope Scope,
          typename KeyEqual,
          typename ProbingScheme,
          typename Allocator,
          typename Storage>
template <typename InputIt, typename OutputIt>
void dynamic_map<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>::contains(
  InputIt first, InputIt last, OutputIt output_begin, cuda::stream_ref stream) const
{
  auto num_keys          = cuco::detail::distance(first, last);
  std::size_t traversed  = 0;
  std::size_t submap_idx = 0;
  while (num_keys > 0 && submap_idx < submaps_.size()) {
    const auto& cur       = submaps_[submap_idx];
    const size_t cur_size = cur->size();
    const size_t num_keys_to_process =
      std::min(static_cast<detail::index_type>(cur_size), num_keys);
    CUCO_CUDA_TRY(cudaStreamSynchronize(stream.get()));

    cur->contains(first, first + num_keys_to_process, output_begin + traversed, stream);

    traversed += num_keys_to_process;
    num_keys -= num_keys_to_process;
    submap_idx++;
    first += num_keys_to_process;
  }
}

}  // namespace experimental
}  // namespace cuco
