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
#pragma once

#include <cuco/detail/utility/cuda.cuh>
#include <cuco/pair.cuh>

#include <cub/block/block_reduce.cuh>
#include <cuda/std/atomic>
#include <cuda/std/iterator>
#include <thrust/type_traits/is_contiguous_iterator.h>

namespace cuco {
namespace detail {
namespace cg = cooperative_groups;

CUCO_SUPPRESS_KERNEL_WARNINGS
/**
 * @brief Initializes each slot in the flat `slots` storage to contain `k` and `v`.
 *
 * Each space in `slots` that can hold a key value pair is initialized to a
 * `pair_atomic_type` containing the key `k` and the value `v`.
 *
 * @tparam atomic_key_type Type of the `Key` atomic container
 * @tparam atomic_mapped_type Type of the `Value` atomic container
 * @tparam Key key type
 * @tparam Value value type
 * @tparam pair_atomic_type key/value pair type
 *
 * @param slots Pointer to flat storage for the map's key/value pairs
 * @param k Key to which all keys in `slots` are initialized
 * @param v Value to which all values in `slots` are initialized
 * @param size Size of the storage pointed to by `slots`
 */
template <typename atomic_key_type,
          typename atomic_mapped_type,
          typename Key,
          typename Value,
          typename pair_atomic_type>
CUCO_KERNEL void initialize(pair_atomic_type* const slots, Key k, Value v, int64_t size)
{
  int64_t const loop_stride = gridDim.x * blockDim.x;
  int64_t idx               = threadIdx.x + blockIdx.x * blockDim.x;
  while (idx < size) {
    new (&slots[idx].first) atomic_key_type{k};
    new (&slots[idx].second) atomic_mapped_type{v};
    idx += loop_stride;
  }
}

/**
 * @brief Inserts all key/value pairs in the range `[first, last)`.
 *
 * Uses the CUDA Cooperative Groups API to leverage groups of multiple threads to perform each
 * key/value insertion. This provides a significant boost in throughput compared to the non
 * Cooperative Group `insert` at moderate to high load factors.
 *
 * @tparam block_size The size of the thread block
 * @tparam tile_size The number of threads in the Cooperative Groups used to perform
 * inserts
 * @tparam InputIt Device accessible random access input iterator where
 * `std::is_convertible<std::iterator_traits<InputIt>::value_type,
 * static_multimap<K, V>::value_type>` is `true`
 * @tparam viewT Type of device view allowing access of hash map storage
 *
 * @param first Beginning of the sequence of key/value pairs
 * @param n Number of key/value pairs to insert
 * @param view Mutable device view used to access the hash map's slot storage
 */
template <uint32_t block_size, uint32_t tile_size, typename InputIt, typename viewT>
CUCO_KERNEL void insert(InputIt first, int64_t n, viewT view)
{
  auto tile = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  int64_t const loop_stride = gridDim.x * block_size / tile_size;
  int64_t idx               = (block_size * blockIdx.x + threadIdx.x) / tile_size;

  while (idx < n) {
    // force conversion to value_type
    typename viewT::value_type const insert_pair{*(first + idx)};
    view.insert(tile, insert_pair);
    idx += loop_stride;
  }
}

/**
 * @brief Inserts key/value pairs in the range `[first, first + n)` if `pred` of the
 * corresponding stencil returns true.
 *
 * The key/value pair `*(first + i)` is inserted if `pred( *(stencil + i) )` returns true.
 *
 * Uses the CUDA Cooperative Groups API to leverage groups of multiple threads to perform each
 * key/value insertion. This provides a significant boost in throughput compared to the non
 * Cooperative Group `insert` at moderate to high load factors.
 *
 * @tparam block_size The size of the thread block
 * @tparam tile_size The number of threads in the Cooperative Groups used to perform
 * inserts
 * @tparam InputIt Device accessible random access input iterator where
 * `std::is_convertible<std::iterator_traits<InputIt>::value_type,
 * static_multimap<K, V>::value_type>` is `true`
 * @tparam StencilIt Device accessible random access iterator whose value_type is
 * convertible to Predicate's argument type
 * @tparam viewT Type of device view allowing access of hash map storage
 * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
 * argument type is convertible from `std::iterator_traits<StencilIt>::value_type`.
 *
 * @param first Beginning of the sequence of key/value pairs
 * @param s Beginning of the stencil sequence
 * @param n Number of elements to insert
 * @param view Mutable device view used to access the hash map's slot storage
 * @param pred Predicate to test on every element in the range `[s, s + n)`
 */
template <uint32_t block_size,
          uint32_t tile_size,
          typename InputIt,
          typename StencilIt,
          typename viewT,
          typename Predicate>
CUCO_KERNEL void insert_if_n(InputIt first, StencilIt s, int64_t n, viewT view, Predicate pred)
{
  auto tile = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  int64_t const loop_stride = gridDim.x * block_size / tile_size;
  int64_t idx               = (block_size * blockIdx.x + threadIdx.x) / tile_size;

  while (idx < n) {
    if (pred(*(s + idx))) {
      typename viewT::value_type const insert_pair{*(first + idx)};
      // force conversion to value_type
      view.insert(tile, insert_pair);
    }
    idx += loop_stride;
  }
}

/**
 * @brief Indicates whether the elements in the range `[first, last)` are contained in the map.
 *
 * Stores `true` or `false` to `(output + i)` indicating if the element `*(first + i)` exists in the
 * map.
 *
 * Uses the CUDA Cooperative Groups API to leverage groups of multiple threads to perform the
 * contains operation for each element. This provides a significant boost in throughput compared
 * to the non Cooperative Group `contains` at moderate to high load factors.
 *
 * @tparam is_pair_contains `true` if it's a `pair_contains` implementation
 * @tparam block_size The size of the thread block
 * @tparam tile_size The number of threads in the Cooperative Groups
 * @tparam InputIt Device accessible input iterator
 * @tparam OutputIt Device accessible output iterator assignable from `bool`
 * @tparam viewT Type of device view allowing access of hash map storage
 * @tparam Equal Binary callable type
 *
 * @param first Beginning of the sequence of elements
 * @param n Number of elements to query
 * @param output_begin Beginning of the sequence of booleans for the presence of each element
 * @param view Device view used to access the hash map's slot storage
 * @param equal The binary function to compare input element and slot content for equality
 */
template <bool is_pair_contains,
          uint32_t block_size,
          uint32_t tile_size,
          typename InputIt,
          typename OutputIt,
          typename viewT,
          typename Equal>
CUCO_KERNEL void contains(InputIt first, int64_t n, OutputIt output_begin, viewT view, Equal equal)
{
  auto tile = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  int64_t const loop_stride = gridDim.x * block_size / tile_size;
  int64_t idx               = (block_size * blockIdx.x + threadIdx.x) / tile_size;
  __shared__ bool writeBuffer[block_size / tile_size];

  while (idx < n) {
    typename cuda::std::iterator_traits<InputIt>::value_type element = *(first + idx);
    auto found                                                       = [&]() {
      if constexpr (is_pair_contains) { return view.pair_contains(tile, element, equal); }
      if constexpr (not is_pair_contains) { return view.contains(tile, element, equal); }
    }();

    /*
     * The ld.relaxed.gpu instruction used in view.find causes L1 to
     * flush more frequently, causing increased sector stores from L2 to global memory.
     * By writing results to shared memory and then synchronizing before writing back
     * to global, we no longer rely on L1, preventing the increase in sector stores from
     * L2 to global and improving performance.
     */
    if (tile.thread_rank() == 0) { writeBuffer[threadIdx.x / tile_size] = found; }
    __syncthreads();
    if (tile.thread_rank() == 0) { *(output_begin + idx) = writeBuffer[threadIdx.x / tile_size]; }
    idx += loop_stride;
  }
}

/**
 * @brief Counts the occurrences of keys in `[first, last)` contained in the multimap.
 *
 * For each key, `k = *(first + i)`, counts all matching keys, `k'`, as determined by `key_equal(k,
 * k')` and stores the sum of all matches for all keys to `num_matches`. If `k` does not have any
 * matches, it contributes 1 to the final sum only if `is_outer` is true.
 *
 * @tparam block_size The size of the thread block
 * @tparam tile_size The number of threads in the Cooperative Groups used to perform counts
 * @tparam uses_vector_load Boolean flag indicating whether vector loads are used or not
 * @tparam is_outer Boolean flag indicating whether non-matches are counted
 * @tparam InputIt Device accessible input iterator whose `value_type` is convertible to the map's
 * `key_type`
 * @tparam atomicT Type of atomic storage
 * @tparam viewT Type of device view allowing access of hash map storage
 * @tparam KeyEqual Binary callable
 *
 * @param first Beginning of the sequence of keys to count
 * @param n Number of the keys to query
 * @param num_matches The number of all the matches for a sequence of keys
 * @param view Device view used to access the hash map's slot storage
 * @param key_equal Binary function to compare two keys for equality
 */
template <uint32_t block_size,
          uint32_t tile_size,
          bool is_outer,
          typename InputIt,
          typename atomicT,
          typename viewT,
          typename KeyEqual>
CUCO_KERNEL void count(
  InputIt first, int64_t n, atomicT* num_matches, viewT view, KeyEqual key_equal)
{
  auto tile = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  int64_t const loop_stride = gridDim.x * block_size / tile_size;
  int64_t idx               = (block_size * blockIdx.x + threadIdx.x) / tile_size;

  typedef cub::BlockReduce<std::size_t, block_size> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  std::size_t thread_num_matches = 0;

  while (idx < n) {
    auto key = *(first + idx);
    if constexpr (is_outer) {
      thread_num_matches += view.count_outer(tile, key, key_equal);
    } else {
      thread_num_matches += view.count(tile, key, key_equal);
    }
    idx += loop_stride;
  }

  // compute number of successfully inserted elements for each block
  // and atomically add to the grand total
  std::size_t block_num_matches = BlockReduce(temp_storage).Sum(thread_num_matches);
  if (threadIdx.x == 0) {
    num_matches->fetch_add(block_num_matches, cuda::std::memory_order_relaxed);
  }
}

/**
 * @brief Counts the occurrences of key/value pairs in `[first, last)` contained in the multimap.
 *
 * For pair, `p = *(first + i)`, counts all matching pairs, `p'`, as determined by `pair_equal(p,
 * p')` and stores the sum of all matches for all pairs to `num_matches`. If `p` does not have any
 * matches, it contributes 1 to the final sum only if `is_outer` is true.
 *
 * @tparam block_size The size of the thread block
 * @tparam tile_size The number of threads in the Cooperative Groups used to perform counts
 * @tparam is_outer Boolean flag indicating whether non-matches are counted
 * @tparam InputIt Device accessible random access input iterator where
 * `std::is_convertible<std::iterator_traits<InputIt>::value_type,
 * static_multimap<K, V>::value_type>` is `true`
 * @tparam atomicT Type of atomic storage
 * @tparam viewT Type of device view allowing access of hash map storage
 * @tparam PairEqual Binary callable
 *
 * @param first Beginning of the sequence of pairs to count
 * @param n Number of the pairs to query
 * @param num_matches The number of all the matches for a sequence of pairs
 * @param view Device view used to access the hash map's slot storage
 * @param pair_equal Binary function to compare two pairs for equality
 */
template <uint32_t block_size,
          uint32_t tile_size,
          bool is_outer,
          typename InputIt,
          typename atomicT,
          typename viewT,
          typename PairEqual>
CUCO_KERNEL void pair_count(
  InputIt first, int64_t n, atomicT* num_matches, viewT view, PairEqual pair_equal)
{
  auto tile = cg::tiled_partition<tile_size, cg::thread_block>(cg::this_thread_block());
  int64_t const loop_stride = gridDim.x * block_size / tile_size;
  int64_t idx               = (block_size * blockIdx.x + threadIdx.x) / tile_size;

  typedef cub::BlockReduce<std::size_t, block_size> BlockReduce;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  std::size_t thread_num_matches = 0;

  while (idx < n) {
    typename viewT::value_type const pair = *(first + idx);
    if constexpr (is_outer) {
      thread_num_matches += view.pair_count_outer(tile, pair, pair_equal);
    } else {
      thread_num_matches += view.pair_count(tile, pair, pair_equal);
    }
    idx += loop_stride;
  }

  // compute number of successfully inserted elements for each block
  // and atomically add to the grand total
  std::size_t block_num_matches = BlockReduce(temp_storage).Sum(thread_num_matches);
  if (threadIdx.x == 0) {
    num_matches->fetch_add(block_num_matches, cuda::std::memory_order_relaxed);
  }
}

/**
 * @brief Retrieves all the values corresponding to all keys in the range `[first, last)`.
 *
 * For key `k = *(first + i)` existing in the map, copies `k` and all associated values to
 * unspecified locations in `[output_begin, output_end)`. If `k` does not have any matches, copies
 * `k` and `empty_value_sentinel()` into the output only if `is_outer` is true.
 *
 * Behavior is undefined if the total number of matching keys exceeds `std::distance(output_begin,
 * output_begin + *num_matches - 1)`. Use `count()` to determine the size of the output range.
 *
 * @tparam block_size The size of the thread block
 * @tparam flushing_cg_size The size of the CG used to flush output buffers
 * @tparam probing_cg_size The size of the CG for parallel retrievals
 * @tparam buffer_size Size of the output buffer
 * @tparam is_outer Boolean flag indicating whether non-matches are included in the output
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the map's `key_type`
 * @tparam OutputIt Device accessible output iterator whose `value_type` is
 * constructible from the map's `value_type`
 * @tparam atomicT Type of atomic storage
 * @tparam viewT Type of device view allowing access of hash map storage
 * @tparam KeyEqual Binary callable type
 *
 * @param first Beginning of the sequence of keys
 * @param n Number of the keys to query
 * @param output_begin Beginning of the sequence of values retrieved for each key
 * @param num_matches Size of the output sequence
 * @param view Device view used to access the hash map's slot storage
 * @param key_equal The binary function to compare two keys for equality
 */
template <uint32_t block_size,
          uint32_t flushing_cg_size,
          uint32_t probing_cg_size,
          uint32_t buffer_size,
          bool is_outer,
          typename InputIt,
          typename OutputIt,
          typename atomicT,
          typename viewT,
          typename KeyEqual>
CUCO_KERNEL void retrieve(InputIt first,
                          int64_t n,
                          OutputIt output_begin,
                          atomicT* num_matches,
                          viewT view,
                          KeyEqual key_equal)
{
  using pair_type = typename viewT::value_type;

  constexpr uint32_t num_flushing_cgs = block_size / flushing_cg_size;
  const uint32_t flushing_cg_id       = threadIdx.x / flushing_cg_size;

  auto flushing_cg =
    cg::tiled_partition<flushing_cg_size, cg::thread_block>(cg::this_thread_block());
  auto probing_cg = cg::tiled_partition<probing_cg_size, cg::thread_block>(cg::this_thread_block());
  int64_t const loop_stride = gridDim.x * block_size / probing_cg_size;
  int64_t idx               = (block_size * blockIdx.x + threadIdx.x) / probing_cg_size;

  __shared__ pair_type output_buffer[num_flushing_cgs][buffer_size];
  // TODO: replace this with shared memory cuda::atomic variables once the dynamiic initialization
  // warning issue is solved __shared__ atomicT counter[num_flushing_cgs][buffer_size];
  __shared__ uint32_t flushing_cg_counter[num_flushing_cgs];

  if (flushing_cg.thread_rank() == 0) { flushing_cg_counter[flushing_cg_id] = 0; }

  flushing_cg.sync();

  while (flushing_cg.any(idx < n)) {
    bool active_flag        = idx < n;
    auto active_flushing_cg = cg::binary_partition<flushing_cg_size>(flushing_cg, active_flag);

    if (active_flag) {
      auto key = *(first + idx);
      if constexpr (is_outer) {
        view.template retrieve_outer<buffer_size>(active_flushing_cg,
                                                  probing_cg,
                                                  key,
                                                  &flushing_cg_counter[flushing_cg_id],
                                                  output_buffer[flushing_cg_id],
                                                  num_matches,
                                                  output_begin,
                                                  key_equal);
      } else {
        view.template retrieve<buffer_size>(active_flushing_cg,
                                            probing_cg,
                                            key,
                                            &flushing_cg_counter[flushing_cg_id],
                                            output_buffer[flushing_cg_id],
                                            num_matches,
                                            output_begin,
                                            key_equal);
      }
    }
    idx += loop_stride;
  }

  flushing_cg.sync();
  // Final flush of output buffer
  if (flushing_cg_counter[flushing_cg_id] > 0) {
    view.flush_output_buffer(flushing_cg,
                             flushing_cg_counter[flushing_cg_id],
                             output_buffer[flushing_cg_id],
                             num_matches,
                             output_begin);
  }
}

/**
 * @brief Retrieves all pairs matching the input probe pair in the range `[first, last)`.
 *
 * If pair_equal(*(first + i), slot[j]) returns true, then *(first+i) is stored to unspecified
 * locations in `probe_output_begin`, and slot[j] is stored to unspecified locations in
 * `contained_output_begin`. If the given pair has no matches in the map, copies *(first + i) in
 * `probe_output_begin` and a pair of `empty_key_sentinel` and `empty_value_sentinel` in
 * `contained_output_begin` only when `is_outer` is `true`.
 *
 * Behavior is undefined if the total number of matching pairs exceeds `std::distance(output_begin,
 * output_begin + *num_matches - 1)`. Use `pair_count()` to determine the size of the output range.
 *
 * @tparam block_size The size of the thread block
 * @tparam flushing_cg_size The size of the CG used to flush output buffers
 * @tparam probing_cg_size The size of the CG for parallel retrievals
 * @tparam buffer_size Size of the output buffer
 * @tparam is_outer Boolean flag indicating whether non-matches are included in the output
 * @tparam InputIt Device accessible random access input iterator where
 * `std::is_convertible<std::iterator_traits<InputIt>::value_type,
 * static_multimap<K, V>::value_type>` is `true`
 * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
 * `InputIt`s `value_type`.
 * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
 * the map's `value_type`.
 * @tparam atomicT Type of atomic storage
 * @tparam viewT Type of device view allowing access of hash map storage
 * @tparam PairEqual Binary callable type
 *
 * @param first Beginning of the sequence of keys
 * @param n Number of keys to query
 * @param probe_output_begin Beginning of the sequence of the matched probe pairs
 * @param contained_output_begin Beginning of the sequence of the matched contained pairs
 * @param num_matches Size of the output sequence
 * @param view Device view used to access the hash map's slot storage
 * @param pair_equal The binary function to compare two pairs for equality
 */
template <uint32_t block_size,
          uint32_t flushing_cg_size,
          uint32_t probing_cg_size,
          uint32_t buffer_size,
          bool is_outer,
          typename InputIt,
          typename OutputIt1,
          typename OutputIt2,
          typename atomicT,
          typename viewT,
          typename PairEqual>
CUCO_KERNEL void pair_retrieve(InputIt first,
                               int64_t n,
                               OutputIt1 probe_output_begin,
                               OutputIt2 contained_output_begin,
                               atomicT* num_matches,
                               viewT view,
                               PairEqual pair_equal)
{
  using pair_type = typename viewT::value_type;

  constexpr uint32_t num_flushing_cgs = block_size / flushing_cg_size;
  const uint32_t flushing_cg_id       = threadIdx.x / flushing_cg_size;

  auto flushing_cg =
    cg::tiled_partition<flushing_cg_size, cg::thread_block>(cg::this_thread_block());
  auto probing_cg = cg::tiled_partition<probing_cg_size, cg::thread_block>(cg::this_thread_block());
  int64_t const loop_stride = gridDim.x * block_size / probing_cg_size;
  int64_t idx               = (block_size * blockIdx.x + threadIdx.x) / probing_cg_size;

  __shared__ pair_type probe_output_buffer[num_flushing_cgs][buffer_size];
  __shared__ pair_type contained_output_buffer[num_flushing_cgs][buffer_size];
  // TODO: replace this with shared memory cuda::atomic variables once the dynamiic initialization
  // warning issue is solved __shared__ atomicT counter[num_flushing_cgs][buffer_size];
  __shared__ uint32_t flushing_cg_counter[num_flushing_cgs];

  if (flushing_cg.thread_rank() == 0) { flushing_cg_counter[flushing_cg_id] = 0; }

  flushing_cg.sync();

  while (flushing_cg.any(idx < n)) {
    bool active_flag        = idx < n;
    auto active_flushing_cg = cg::binary_partition<flushing_cg_size>(flushing_cg, active_flag);

    if (active_flag) {
      pair_type pair = *(first + idx);
      if constexpr (is_outer) {
        view.pair_retrieve_outer<buffer_size>(active_flushing_cg,
                                              probing_cg,
                                              pair,
                                              &flushing_cg_counter[flushing_cg_id],
                                              probe_output_buffer[flushing_cg_id],
                                              contained_output_buffer[flushing_cg_id],
                                              num_matches,
                                              probe_output_begin,
                                              contained_output_begin,
                                              pair_equal);
      } else {
        view.pair_retrieve<buffer_size>(active_flushing_cg,
                                        probing_cg,
                                        pair,
                                        &flushing_cg_counter[flushing_cg_id],
                                        probe_output_buffer[flushing_cg_id],
                                        contained_output_buffer[flushing_cg_id],
                                        num_matches,
                                        probe_output_begin,
                                        contained_output_begin,
                                        pair_equal);
      }
    }
    idx += loop_stride;
  }

  flushing_cg.sync();
  // Final flush of output buffer
  if (flushing_cg_counter[flushing_cg_id] > 0) {
    view.flush_output_buffer(flushing_cg,
                             flushing_cg_counter[flushing_cg_id],
                             probe_output_buffer[flushing_cg_id],
                             contained_output_buffer[flushing_cg_id],
                             num_matches,
                             probe_output_begin,
                             contained_output_begin);
  }
}
}  // namespace detail
}  // namespace cuco
