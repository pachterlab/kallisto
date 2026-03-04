/*
 * Copyright (c) 2020-2025, NVIDIA CORPORATION.
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

#include <cuco/detail/bitwise_compare.cuh>
#include <cuco/detail/utility/cuda.cuh>

#include <cub/block/block_reduce.cuh>
#include <cub/block/block_scan.cuh>
#include <cuda/std/atomic>

#include <cooperative_groups.h>

namespace cuco {
namespace detail {
namespace cg = cooperative_groups;

CUCO_SUPPRESS_KERNEL_WARNINGS

/**
 * @brief Inserts all key/value pairs in the range `[first, last)`.
 *
 * If multiple keys in `[first, last)` compare equal, it is unspecified which
 * element is inserted.
 *
 * @tparam block_size
 * @tparam pair_type Type of the pairs contained in the map
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `value_type`
 * @tparam viewT Type of the `static_map` device views
 * @tparam mutableViewT Type of the `static_map` device mutable views
 * @tparam atomicT Type of atomic storage
 * @tparam viewT Type of device view allowing access of hash map storage
 * @tparam Hash Unary callable type
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of key/value pairs
 * @param last End of the sequence of key/value pairs
 * @param submap_views Array of `static_map::device_view` objects used to
 * perform `contains` operations on each underlying `static_map`
 * @param submap_mutable_views Array of `static_map::device_mutable_view` objects
 * used to perform an `insert` into the target `static_map` submap
 * @param num_successes The number of successfully inserted key/value pairs
 * @param insert_idx The index of the submap we are inserting into
 * @param num_submaps The total number of submaps in the map
 * @param hash The unary function to apply to hash each key
 * @param key_equal The binary function used to compare two keys for equality
 */
template <uint32_t block_size,
          typename pair_type,
          typename InputIt,
          typename viewT,
          typename mutableViewT,
          typename atomicT,
          typename Hash,
          typename KeyEqual>
CUCO_KERNEL void insert(InputIt first,
                        InputIt last,
                        viewT* submap_views,
                        mutableViewT* submap_mutable_views,
                        atomicT* num_successes,
                        uint32_t insert_idx,
                        uint32_t num_submaps,
                        Hash hash,
                        KeyEqual key_equal)
{
  using BlockReduce = cub::BlockReduce<std::size_t, block_size>;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  std::size_t thread_num_successes = 0;

  auto tid = blockDim.x * blockIdx.x + threadIdx.x;

  while (first + tid < last) {
    pair_type insert_pair = *(first + tid);
    auto exists           = false;

    // manually check for duplicates in those submaps we are not inserting into
    for (auto i = 0; i < num_submaps; ++i) {
      if (i != insert_idx) {
        exists = submap_views[i].contains(insert_pair.first, hash, key_equal);
        if (exists) { break; }
      }
    }
    if (!exists) {
      if (submap_mutable_views[insert_idx].insert(insert_pair, hash, key_equal)) {
        thread_num_successes++;
      }
    }

    tid += gridDim.x * blockDim.x;
  }

  std::size_t const block_num_successes = BlockReduce(temp_storage).Sum(thread_num_successes);
  if (threadIdx.x == 0) {
    num_successes->fetch_add(block_num_successes, cuda::std::memory_order_relaxed);
  }
}

/**
 * @brief Inserts all key/value pairs in the range `[first, last)`.
 *
 * If multiple keys in `[first, last)` compare equal, it is unspecified which
 * element is inserted. Uses the CUDA Cooperative Groups API to leverage groups
 * of multiple threads to perform each key/value insertion. This provides a
 * significant boost in throughput compared to the non Cooperative Group
 * `insert` at moderate to high load factors.
 *
 * @tparam block_size
 * @tparam tile_size The number of threads in the Cooperative Groups used to perform
 * inserts
 * @tparam pair_type Type of the pairs contained in the map
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `value_type`
 * @tparam viewT Type of the `static_map` device views
 * @tparam mutableViewT Type of the `static_map` device mutable views
 * @tparam atomicT Type of atomic storage
 * @tparam viewT Type of device view allowing access of hash map storage
 * @tparam Hash Unary callable type
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of key/value pairs
 * @param last End of the sequence of key/value pairs
 * @param submap_views Array of `static_map::device_view` objects used to
 * perform `contains` operations on each underlying `static_map`
 * @param submap_mutable_views Array of `static_map::device_mutable_view` objects
 * used to perform an `insert` into the target `static_map` submap
 * @param submap_num_successes The number of successfully inserted key/value pairs for each submap
 * @param insert_idx The index of the submap we are inserting into
 * @param num_submaps The total number of submaps in the map
 * @param hash The unary function to apply to hash each key
 * @param key_equal The binary function used to compare two keys for equality
 */
template <uint32_t block_size,
          uint32_t tile_size,
          typename pair_type,
          typename InputIt,
          typename viewT,
          typename mutableViewT,
          typename atomicT,
          typename Hash,
          typename KeyEqual>
CUCO_KERNEL void insert(InputIt first,
                        InputIt last,
                        viewT* submap_views,
                        mutableViewT* submap_mutable_views,
                        atomicT** submap_num_successes,
                        uint32_t insert_idx,
                        uint32_t num_submaps,
                        Hash hash,
                        KeyEqual key_equal)
{
  using BlockReduce = cub::BlockReduce<std::size_t, block_size>;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  std::size_t thread_num_successes = 0;

  auto tile = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  auto tid  = blockDim.x * blockIdx.x + threadIdx.x;
  auto it   = first + tid / tile_size;

  while (it < last) {
    pair_type insert_pair = *it;
    auto exists           = false;

    // manually check for duplicates in those submaps we are not inserting into
    for (auto i = 0; i < num_submaps; ++i) {
      if (i != insert_idx) {
        exists = submap_views[i].contains(tile, insert_pair.first, hash, key_equal);
        if (exists) { break; }
      }
    }
    if (!exists) {
      if (submap_mutable_views[insert_idx].insert(tile, insert_pair, hash, key_equal) &&
          tile.thread_rank() == 0) {
        thread_num_successes++;
      }
    }

    it += (gridDim.x * blockDim.x) / tile_size;
  }

  std::size_t const block_num_successes = BlockReduce(temp_storage).Sum(thread_num_successes);
  if (threadIdx.x == 0) {
    submap_num_successes[insert_idx]->fetch_add(block_num_successes,
                                                cuda::std::memory_order_relaxed);
  }
}

/**
 * @brief Erases the key/value pairs corresponding to all keys in the range `[first, last)`.
 *
 * If the key `*(first + i)` exists in the map, its slot is erased and made available for future
   insertions.
 * Else, no effect.
 *
 * @tparam block_size The size of the thread block
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `key_type`
 * @tparam mutableViewT Type of device view allowing modification of hash map storage
 * @tparam atomicT Type of atomic storage
 * @tparam Hash Unary callable type
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of keys
 * @param last End of the sequence of keys
 * @param submap_mutable_views Array of `static_map::mutable_device_view` objects used to
 * perform `erase` operations on each underlying `static_map`
 * @param num_successes The number of successfully erased key/value pairs
 * @param submap_num_successes The number of successfully erased key/value pairs
 * in each submap
 * @param num_submaps The number of submaps in the map
 * @param hash The unary function to apply to hash each key
 * @param key_equal The binary function to compare two keys for equality
 */
template <uint32_t block_size,
          typename InputIt,
          typename mutableViewT,
          typename atomicT,
          typename Hash,
          typename KeyEqual>
CUCO_KERNEL void erase(InputIt first,
                       InputIt last,
                       mutableViewT* submap_mutable_views,
                       atomicT** submap_num_successes,
                       uint32_t num_submaps,
                       Hash hash,
                       KeyEqual key_equal)
{
  extern __shared__ unsigned long long submap_block_num_successes[];

  auto tid = block_size * blockIdx.x + threadIdx.x;
  auto it  = first + tid;

  for (auto i = threadIdx.x; i < num_submaps; i += block_size) {
    submap_block_num_successes[i] = 0;
  }
  __syncthreads();

  while (it < last) {
    for (auto i = 0; i < num_submaps; ++i) {
      if (submap_mutable_views[i].erase(*it, hash, key_equal)) {
        atomicAdd(&submap_block_num_successes[i], 1);
        break;
      }
    }
    it += gridDim.x * blockDim.x;
  }
  __syncthreads();

  for (auto i = 0; i < num_submaps; ++i) {
    if (threadIdx.x == 0) {
      submap_num_successes[i]->fetch_add(static_cast<std::size_t>(submap_block_num_successes[i]),
                                         cuda::std::memory_order_relaxed);
    }
  }
}

/**
 * @brief Erases the key/value pairs corresponding to all keys in the range `[first, last)`.
 *
 * If the key `*(first + i)` exists in the map, its slot is erased and made available for future
 * insertions.
 * Else, no effect.
 *
 * @tparam block_size The size of the thread block
 * @tparam tile_size The number of threads in the Cooperative Groups used to perform erase
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `key_type`
 * @tparam mutableViewT Type of device view allowing modification of hash map storage
 * @tparam atomicT Type of atomic storage
 * @tparam Hash Unary callable type
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of keys
 * @param last End of the sequence of keys
 * @param submap_mutable_views Array of `static_map::mutable_device_view` objects used to
 * perform `erase` operations on each underlying `static_map`
 * @param num_successes The number of successfully erased key/value pairs
 * @param submap_num_successes The number of successfully erased key/value pairs
 * in each submap
 * @param num_submaps The number of submaps in the map
 * @param hash The unary function to apply to hash each key
 * @param key_equal The binary function to compare two keys for equality
 */
template <uint32_t block_size,
          uint32_t tile_size,
          typename InputIt,
          typename mutableViewT,
          typename atomicT,
          typename Hash,
          typename KeyEqual>
CUCO_KERNEL void erase(InputIt first,
                       InputIt last,
                       mutableViewT* submap_mutable_views,
                       atomicT** submap_num_successes,
                       uint32_t num_submaps,
                       Hash hash,
                       KeyEqual key_equal)
{
  extern __shared__ unsigned long long submap_block_num_successes[];

  auto block = cg::this_thread_block();
  auto tile  = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  auto tid   = block_size * block.group_index().x + block.thread_rank();
  auto it    = first + tid / tile_size;

  for (auto i = threadIdx.x; i < num_submaps; i += block_size) {
    submap_block_num_successes[i] = 0;
  }
  block.sync();

  while (it < last) {
    auto erased = false;
    int i       = 0;
    for (i = 0; i < num_submaps; ++i) {
      erased = submap_mutable_views[i].erase(tile, *it, hash, key_equal);
      if (erased) { break; }
    }
    if (erased && tile.thread_rank() == 0) { atomicAdd(&submap_block_num_successes[i], 1); }
    it += (gridDim.x * blockDim.x) / tile_size;
  }
  block.sync();

  for (auto i = 0; i < num_submaps; ++i) {
    if (threadIdx.x == 0) {
      submap_num_successes[i]->fetch_add(static_cast<std::size_t>(submap_block_num_successes[i]),
                                         cuda::std::memory_order_relaxed);
    }
  }
}

/**
 * @brief Finds the values corresponding to all keys in the range `[first, last)`.
 *
 * If the key `*(first + i)` exists in the map, copies its associated value to `(output_begin + i)`.
 * Else, copies the empty value sentinel.
 *
 * @tparam block_size The number of threads in the thread block
 * @tparam Value The mapped value type for the map
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `key_type`
 * @tparam OutputIt Device accessible output iterator whose `value_type` is
 * convertible to the map's `mapped_type`
 * @tparam viewT Type of `static_map` device view
 * @tparam Hash Unary callable type
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of keys
 * @param last End of the sequence of keys
 * @param output_begin Beginning of the sequence of values retrieved for each key
 * @param submap_views Array of `static_map::device_view` objects used to
 * perform `find` operations on each underlying `static_map`
 * @param num_submaps The number of submaps in the map
 * @param hash The unary function to apply to hash each key
 * @param key_equal The binary function to compare two keys for equality
 */
template <uint32_t block_size,
          typename Value,
          typename InputIt,
          typename OutputIt,
          typename viewT,
          typename Hash,
          typename KeyEqual>
CUCO_KERNEL void find(InputIt first,
                      InputIt last,
                      OutputIt output_begin,
                      viewT* submap_views,
                      uint32_t num_submaps,
                      Hash hash,
                      KeyEqual key_equal)
{
  auto tid                  = blockDim.x * blockIdx.x + threadIdx.x;
  auto empty_value_sentinel = submap_views[0].get_empty_value_sentinel();
  __shared__ Value writeBuffer[block_size];

  while (first + tid < last) {
    auto key         = *(first + tid);
    auto found_value = empty_value_sentinel;
    for (auto i = 0; i < num_submaps; ++i) {
      auto submap_view = submap_views[i];
      auto found       = submap_view.find(key, hash, key_equal);
      if (found != submap_view.end()) {
        found_value = found->second.load(cuda::std::memory_order_relaxed);
        break;
      }
    }

    /*
     * The ld.relaxed.gpu instruction used in view.find causes L1 to
     * flush more frequently, causing increased sector stores from L2 to global memory.
     * By writing results to shared memory and then synchronizing before writing back
     * to global, we no longer rely on L1, preventing the increase in sector stores from
     * L2 to global and improving performance.
     */
    writeBuffer[threadIdx.x] = found_value;
    __syncthreads();
    *(output_begin + tid) = writeBuffer[threadIdx.x];
    tid += gridDim.x * blockDim.x;
  }
}

/**
 * @brief Finds the values corresponding to all keys in the range `[first, last)`.
 *
 * If the key `*(first + i)` exists in the map, copies its associated value to `(output_begin + i)`.
 * Else, copies the empty value sentinel. Uses the CUDA Cooperative Groups API to leverage groups
 * of multiple threads to find each key. This provides a significant boost in throughput compared
 * to the non Cooperative Group `find` at moderate to high load factors.
 *
 * @tparam block_size The number of threads in the thread block
 * @tparam tile_size The number of threads in the Cooperative Groups used to
 * perform find operations
 * @tparam Value The mapped value type for the map
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `key_type`
 * @tparam OutputIt Device accessible output iterator whose `value_type` is
 * convertible to the map's `mapped_type`
 * @tparam viewT Type of `static_map` device view
 * @tparam Hash Unary callable type
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of keys
 * @param last End of the sequence of keys
 * @param output_begin Beginning of the sequence of values retrieved for each key
 * @param submap_views Array of `static_map::device_view` objects used to
 * perform `find` operations on each underlying `static_map`
 * @param num_submaps The number of submaps in the map
 * @param hash The unary function to apply to hash each key
 * @param key_equal The binary function to compare two keys for equality
 */
template <uint32_t block_size,
          uint32_t tile_size,
          typename Value,
          typename InputIt,
          typename OutputIt,
          typename viewT,
          typename Hash,
          typename KeyEqual>
CUCO_KERNEL void find(InputIt first,
                      InputIt last,
                      OutputIt output_begin,
                      viewT* submap_views,
                      uint32_t num_submaps,
                      Hash hash,
                      KeyEqual key_equal)
{
  auto tile    = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  auto tid     = blockDim.x * blockIdx.x + threadIdx.x;
  auto key_idx = tid / tile_size;
  auto empty_value_sentinel = submap_views[0].get_empty_value_sentinel();
  __shared__ Value writeBuffer[block_size];

  while (first + key_idx < last) {
    auto key         = *(first + key_idx);
    auto found_value = empty_value_sentinel;
    for (auto i = 0; i < num_submaps; ++i) {
      auto submap_view = submap_views[i];
      auto found       = submap_view.find(tile, key, hash, key_equal);
      if (found != submap_view.end()) {
        found_value = found->second.load(cuda::std::memory_order_relaxed);
        break;
      }
    }

    /*
     * The ld.relaxed.gpu instruction used in view.find causes L1 to
     * flush more frequently, causing increased sector stores from L2 to global memory.
     * By writing results to shared memory and then synchronizing before writing back
     * to global, we no longer rely on L1, preventing the increase in sector stores from
     * L2 to global and improving performance.
     */
    if (tile.thread_rank() == 0) { writeBuffer[threadIdx.x / tile_size] = found_value; }
    __syncthreads();
    if (tile.thread_rank() == 0) {
      *(output_begin + key_idx) = writeBuffer[threadIdx.x / tile_size];
    }
    key_idx += (gridDim.x * blockDim.x) / tile_size;
  }
}

/**
 * @brief Retrieves all of the keys and their associated values.
 *
 * The order in which keys are returned is implementation defined and not guaranteed to be
 * consistent between subsequent calls to `retrieve_all`.
 *
 * Behavior is undefined if the range beginning at `keys_out` or `values_out` is less than
 * `get_size()`
 *
 * @tparam block_size The number of threads in the thread block
 * @tparam KeyOutputIt Device accessible output iterator whose `value_type` is
 * convertible to the map's `key_type`
 * @tparam ValueOutputIt Device accessible output iterator whose `value_type` is
 * convertible to the map's `mapped_type`
 * @tparam viewT Type of `static_map` device view
 * @tparam AtomicT Atomic counter type
 *
 * @param keys_out Beginning output iterator for keys
 * @param values_out Beginning output iterator for values
 * @param submap_views Array of `static_map::device_view` objects used to
 * perform `retrieve_all` operations on each underlying `static_map`
 * @param num_submaps The number of submaps in the map
 * @param capacity The total number of slots of all submaps
 * @param d_num_out Pointer to the device memory location where the number of keys/vals retrieved
 * are stored
 * @param cap_prefix_sum Array of prefix sums of the number of slots in each submap
 * @return Pair of iterators indicating the last elements in the output
 */
template <uint32_t block_size,
          typename KeyOutputIt,
          typename ValueOutputIt,
          typename viewT,
          typename AtomicT>
CUCO_KERNEL void retrieve_all(KeyOutputIt keys_out,
                              ValueOutputIt values_out,
                              viewT* submap_views,
                              uint32_t num_submaps,
                              uint64_t capacity,
                              AtomicT* d_num_out,
                              size_t* cap_prefix_sum)
{
  using BlockScan = cub::BlockScan<unsigned int, block_size>;

  __shared__ typename BlockScan::TempStorage scan_temp_storage;
  __shared__ unsigned int block_base;

  auto tid                       = blockDim.x * blockIdx.x + threadIdx.x;
  auto const empty_key_sentinel  = submap_views[0].get_empty_key_sentinel();
  auto const erased_key_sentinel = submap_views[0].get_erased_key_sentinel();

  while ((tid - threadIdx.x) < capacity) {
    uint32_t submap_idx    = 0;
    uint32_t submap_offset = tid;
    while (tid >= cap_prefix_sum[submap_idx] && submap_idx < num_submaps) {
      ++submap_idx;
    }
    if (submap_idx > 0) { submap_offset = tid - cap_prefix_sum[submap_idx - 1]; }

    auto const& current_slot = submap_views[submap_idx].get_slots()[submap_offset];
    auto const existing_key  = current_slot.first.load(cuda::std::memory_order_relaxed);

    bool const is_filled = not(cuco::detail::bitwise_compare(existing_key, empty_key_sentinel) or
                               cuco::detail::bitwise_compare(existing_key, erased_key_sentinel));

    unsigned int local_idx   = 0;
    unsigned int block_valid = 0;
    BlockScan(scan_temp_storage).ExclusiveSum(is_filled ? 1u : 0u, local_idx, block_valid);

    if (threadIdx.x == 0) {
      block_base = d_num_out->fetch_add(block_valid, cuda::memory_order_relaxed);
    }
    __syncthreads();

    if (is_filled) {
      auto const value                 = current_slot.second.load(cuda::std::memory_order_relaxed);
      keys_out[block_base + local_idx] = existing_key;
      values_out[block_base + local_idx] = value;
    }
    tid += gridDim.x * blockDim.x;
  }
}

/**
 * @brief Indicates whether the keys in the range `[first, last)` are contained in the map.
 *
 * Writes a `bool` to `(output + i)` indicating if the key `*(first + i)` exists in the map.
 *
 * @tparam block_size The number of threads in the thread block
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `key_type`
 * @tparam OutputIt Device accessible output iterator whose `value_type` is
 * convertible to the map's `mapped_type`
 * @tparam viewT Type of `static_map` device view
 * @tparam Hash Unary callable type
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of keys
 * @param last End of the sequence of keys
 * @param output_begin Beginning of the sequence of booleans for the presence of each key
 * @param submap_views Array of `static_map::device_view` objects used to
 * perform `contains` operations on each underlying `static_map`
 * @param num_submaps The number of submaps in the map
 * @param hash The unary function to apply to hash each key
 * @param key_equal The binary function to compare two keys for equality
 */
template <uint32_t block_size,
          typename InputIt,
          typename OutputIt,
          typename viewT,
          typename Hash,
          typename KeyEqual>
CUCO_KERNEL void contains(InputIt first,
                          InputIt last,
                          OutputIt output_begin,
                          viewT* submap_views,
                          uint32_t num_submaps,
                          Hash hash,
                          KeyEqual key_equal)
{
  auto tid = blockDim.x * blockIdx.x + threadIdx.x;
  __shared__ bool writeBuffer[block_size];

  while (first + tid < last) {
    auto key   = *(first + tid);
    auto found = false;
    for (auto i = 0; i < num_submaps; ++i) {
      found = submap_views[i].contains(key, hash, key_equal);
      if (found) { break; }
    }

    /*
     * The ld.relaxed.gpu instruction used in view.find causes L1 to
     * flush more frequently, causing increased sector stores from L2 to global memory.
     * By writing results to shared memory and then synchronizing before writing back
     * to global, we no longer rely on L1, preventing the increase in sector stores from
     * L2 to global and improving performance.
     */
    writeBuffer[threadIdx.x] = found;
    __syncthreads();
    *(output_begin + tid) = writeBuffer[threadIdx.x];
    tid += gridDim.x * blockDim.x;
  }
}

/**
 * @brief Indicates whether the keys in the range `[first, last)` are contained in the map.
 *
 * Writes a `bool` to `(output + i)` indicating if the key `*(first + i)` exists in the map.
 * Uses the CUDA Cooperative Groups API to leverage groups of multiple threads to perform the
 * contains operation for each key. This provides a significant boost in throughput compared
 * to the non Cooperative Group `contains` at moderate to high load factors.
 *
 * @tparam block_size The number of threads in the thread block
 * @tparam tile_size The number of threads in the Cooperative Groups used to
 * perform find operations
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `key_type`
 * @tparam OutputIt Device accessible output iterator whose `value_type` is
 * convertible to the map's `mapped_type`
 * @tparam viewT Type of `static_map` device view
 * @tparam Hash Unary callable type
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of keys
 * @param last End of the sequence of keys
 * @param output_begin Beginning of the sequence of booleans for the presence of each key
 * @param submap_views Array of `static_map::device_view` objects used to
 * perform `contains` operations on each underlying `static_map`
 * @param num_submaps The number of submaps in the map
 * @param hash The unary function to apply to hash each key
 * @param key_equal The binary function to compare two keys for equality
 */
template <uint32_t block_size,
          uint32_t tile_size,
          typename InputIt,
          typename OutputIt,
          typename viewT,
          typename Hash,
          typename KeyEqual>
CUCO_KERNEL void contains(InputIt first,
                          InputIt last,
                          OutputIt output_begin,
                          viewT* submap_views,
                          uint32_t num_submaps,
                          Hash hash,
                          KeyEqual key_equal)
{
  auto tile    = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  auto tid     = blockDim.x * blockIdx.x + threadIdx.x;
  auto key_idx = tid / tile_size;
  __shared__ bool writeBuffer[block_size];

  while (first + key_idx < last) {
    auto key   = *(first + key_idx);
    auto found = false;
    for (auto i = 0; i < num_submaps; ++i) {
      found = submap_views[i].contains(tile, key, hash, key_equal);
      if (found) { break; }
    }

    /*
     * The ld.relaxed.gpu instruction used in view.find causes L1 to
     * flush more frequently, causing increased sector stores from L2 to global memory.
     * By writing results to shared memory and then synchronizing before writing back
     * to global, we no longer rely on L1, preventing the increase in sector stores from
     * L2 to global and improving performance.
     */
    if (tile.thread_rank() == 0) { writeBuffer[threadIdx.x / tile_size] = found; }
    __syncthreads();
    if (tile.thread_rank() == 0) {
      *(output_begin + key_idx) = writeBuffer[threadIdx.x / tile_size];
    }
    key_idx += (gridDim.x * blockDim.x) / tile_size;
  }
}
}  // namespace detail
}  // namespace cuco
