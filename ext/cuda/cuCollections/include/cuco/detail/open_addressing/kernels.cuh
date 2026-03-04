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

#include <cuco/detail/utility/cuda.cuh>

#include <cub/block/block_reduce.cuh>
#include <cuda/atomic>
#include <cuda/functional>
#include <cuda/std/iterator>
#include <cuda/std/type_traits>

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

namespace cuco::detail::open_addressing_ns {
CUCO_SUPPRESS_KERNEL_WARNINGS

/**
 * @brief Inserts all elements in the range `[first, first + n)` and returns the number of
 * successful insertions if `pred` of the corresponding stencil returns true.
 *
 * @note If multiple elements in `[first, first + n)` compare equal, it is unspecified which element
 * is inserted.
 * @note The key `*(first + i)` is inserted if `pred( *(stencil + i) )` returns true.
 *
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the `value_type` of the data structure
 * @tparam StencilIt Device accessible random access iterator whose value_type is
 * convertible to Predicate's argument type
 * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool`
 * and argument type is convertible from `std::iterator_traits<StencilIt>::value_type`
 * @tparam AtomicT Atomic counter type
 * @tparam Ref Type of non-owning device container ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param stencil Beginning of the stencil sequence
 * @param pred Predicate to test on every element in the range `[stencil, stencil + n)`
 * @param num_successes Number of successful inserted elements
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <int CGSize,
          int BlockSize,
          typename InputIt,
          typename StencilIt,
          typename Predicate,
          typename AtomicT,
          typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void insert_if_n(InputIt first,
                                                          cuco::detail::index_type n,
                                                          StencilIt stencil,
                                                          Predicate pred,
                                                          AtomicT* num_successes,
                                                          Ref ref)
{
  using BlockReduce = cub::BlockReduce<typename Ref::size_type, BlockSize>;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  typename Ref::size_type thread_num_successes = 0;

  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  while (idx < n) {
    if (pred(*(stencil + idx))) {
      typename cuda::std::iterator_traits<InputIt>::value_type const insert_element{*(first + idx)};
      if constexpr (CGSize == 1) {
        if (ref.insert(insert_element)) { thread_num_successes++; };
      } else {
        auto const tile =
          cooperative_groups::tiled_partition<CGSize, cooperative_groups::thread_block>(
            cooperative_groups::this_thread_block());
        if (ref.insert(tile, insert_element) && tile.thread_rank() == 0) { thread_num_successes++; }
      }
    }
    idx += loop_stride;
  }

  // compute number of successfully inserted elements for each block
  // and atomically add to the grand total
  auto const block_num_successes = BlockReduce(temp_storage).Sum(thread_num_successes);
  if (threadIdx.x == 0) {
    num_successes->fetch_add(block_num_successes, cuda::std::memory_order_relaxed);
  }
}

/**
 * @brief Inserts all elements in the range `[first, first + n)` if `pred` of the corresponding
 * stencil returns true.
 *
 * @note If multiple elements in `[first, first + n)` compare equal, it is unspecified which element
 * is inserted.
 * @note The key `*(first + i)` is inserted if `pred( *(stencil + i) )` returns true.
 *
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the `value_type` of the data structure
 * @tparam StencilIt Device accessible random access iterator whose value_type is
 * convertible to Predicate's argument type
 * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool`
 * and argument type is convertible from `std::iterator_traits<StencilIt>::value_type`
 * @tparam Ref Type of non-owning device ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param stencil Beginning of the stencil sequence
 * @param pred Predicate to test on every element in the range `[stencil, stencil + n)`
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <int CGSize,
          int BlockSize,
          typename InputIt,
          typename StencilIt,
          typename Predicate,
          typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void insert_if_n(
  InputIt first, cuco::detail::index_type n, StencilIt stencil, Predicate pred, Ref ref)
{
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  while (idx < n) {
    if (pred(*(stencil + idx))) {
      typename cuda::std::iterator_traits<InputIt>::value_type const insert_element{*(first + idx)};
      if constexpr (CGSize == 1) {
        ref.insert(insert_element);
      } else {
        auto const tile =
          cooperative_groups::tiled_partition<CGSize, cooperative_groups::thread_block>(
            cooperative_groups::this_thread_block());
        ref.insert(tile, insert_element);
      }
    }
    idx += loop_stride;
  }
}

/**
 * @brief Asynchronously erases keys in the range `[first, first + n)`.
 *
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the `value_type` of the data structure
 * @tparam Ref Type of non-owning device ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <int CGSize, int BlockSize, typename InputIt, typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void erase(InputIt first,
                                                    cuco::detail::index_type n,
                                                    Ref ref)
{
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  while (idx < n) {
    typename cuda::std::iterator_traits<InputIt>::value_type const erase_element{*(first + idx)};
    if constexpr (CGSize == 1) {
      ref.erase(erase_element);
    } else {
      auto const tile =
        cooperative_groups::tiled_partition<CGSize, cooperative_groups::thread_block>(
          cooperative_groups::this_thread_block());
      ref.erase(tile, erase_element);
    }
    idx += loop_stride;
  }
}

/**
 * @brief For each key in the range [first, first + n), applies the function object `callback_op` to
 * the copy of all corresponding matches found in the container.
 *
 * @note The return value of `callback_op`, if any, is ignored.
 *
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the `key_type` of the data structure
 * @tparam CallbackOp Type of unary callback function object
 * @tparam Ref Type of non-owning device ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param callback_op Function to call on every matched slot found in the container
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <int CGSize, int BlockSize, typename InputIt, typename CallbackOp, typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void for_each_n(InputIt first,
                                                         cuco::detail::index_type n,
                                                         CallbackOp callback_op,
                                                         Ref ref)
{
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  while (idx < n) {
    typename cuda::std::iterator_traits<InputIt>::value_type const key{*(first + idx)};
    if constexpr (CGSize == 1) {
      ref.for_each(key, callback_op);
    } else {
      auto const tile =
        cooperative_groups::tiled_partition<CGSize, cooperative_groups::thread_block>(
          cooperative_groups::this_thread_block());
      ref.for_each(tile, key, callback_op);
    }
    idx += loop_stride;
  }
}

/**
 * @brief Indicates whether the keys in the range `[first, first + n)` are contained in the data
 * structure if `pred` of the corresponding stencil returns true.
 *
 * @note If `pred( *(stencil + i) )` is true, stores `true` or `false` to `(output_begin + i)`
 * indicating if the key `*(first + i)` is present in the container. If `pred( *(stencil + i) )` is
 * false, stores false to `(output_begin + i)`.
 *
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize The size of the thread block
 * @tparam InputIt Device accessible input iterator
 * @tparam StencilIt Device accessible random access iterator whose value_type is
 * convertible to Predicate's argument type
 * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool`
 * and argument type is convertible from `std::iterator_traits<StencilIt>::value_type`
 * @tparam OutputIt Device accessible output iterator assignable from `bool`
 * @tparam Ref Type of non-owning device ref allowing access to storage
 *
 * @param first Beginning of the sequence of keys
 * @param n Number of keys
 * @param stencil Beginning of the stencil sequence
 * @param pred Predicate to test on every element in the range `[stencil, stencil + n)`
 * @param output_begin Beginning of the sequence of booleans for the presence of each key
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <int CGSize,
          int BlockSize,
          typename InputIt,
          typename StencilIt,
          typename Predicate,
          typename OutputIt,
          typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void contains_if_n(InputIt first,
                                                            cuco::detail::index_type n,
                                                            StencilIt stencil,
                                                            Predicate pred,
                                                            OutputIt output_begin,
                                                            Ref ref)
{
  namespace cg = cooperative_groups;

  auto const block       = cg::this_thread_block();
  auto const thread_idx  = block.thread_rank();
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  __shared__ bool output_buffer[BlockSize / CGSize];

  while ((idx - thread_idx / CGSize) < n) {  // the whole thread block falls into the same iteration
    if constexpr (CGSize == 1) {
      if (idx < n) {
        typename cuda::std::iterator_traits<InputIt>::value_type const key = *(first + idx);
        /*
         * The ld.relaxed.gpu instruction causes L1 to flush more frequently, causing increased
         * sector stores from L2 to global memory. By writing results to shared memory and then
         * synchronizing before writing back to global, we no longer rely on L1, preventing the
         * increase in sector stores from L2 to global and improving performance.
         */
        output_buffer[thread_idx] = pred(*(stencil + idx)) ? ref.contains(key) : false;
      }
      block.sync();
      if (idx < n) { *(output_begin + idx) = output_buffer[thread_idx]; }
    } else {
      auto const tile = cg::tiled_partition<CGSize, cg::thread_block>(block);
      if (idx < n) {
        typename cuda::std::iterator_traits<InputIt>::value_type const key = *(first + idx);
        auto const found = pred(*(stencil + idx)) ? ref.contains(tile, key) : false;
        if (tile.thread_rank() == 0) { *(output_begin + idx) = found; }
      }
    }
    idx += loop_stride;
  }
}

/**
 * @brief Helper to determine the buffer type for the find kernel
 *
 * @tparam Container Container type
 */
template <typename Container, typename = void>
struct find_buffer {
  using type = typename Container::key_type;  ///< Buffer type
};

/**
 * @brief Helper to determine the buffer type for the find kernel
 *
 * @note Specialization if `mapped_type` exists
 *
 * @tparam Container Container type
 */
template <typename Container>
struct find_buffer<Container, cuda::std::void_t<typename Container::mapped_type>> {
  using type = typename Container::mapped_type;  ///< Buffer type
};

/**
 * @brief Finds the equivalent container elements of all keys in the range `[first, first + n)`.
 *
 * @note If the key `*(first + i)` has a match in the container, copies the match to `(output_begin
 * + i)`. Else, copies the empty sentinel. Uses the CUDA Cooperative Groups API to leverage groups
 * of multiple threads to find each key. This provides a significant boost in throughput compared to
 * the non Cooperative Group `find` at moderate to high load factors.
 *
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize The size of the thread block
 * @tparam InputIt Device accessible input iterator
 * @tparam StencilIt Device accessible random access iterator whose value_type is
 * convertible to Predicate's argument type
 * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool`
 * and argument type is convertible from `std::iterator_traits<StencilIt>::value_type`
 * @tparam OutputIt Device accessible output iterator
 * @tparam Ref Type of non-owning device ref allowing access to storage
 *
 * @param first Beginning of the sequence of keys
 * @param n Number of keys to query
 * @param stencil Beginning of the stencil sequence
 * @param pred Predicate to test on every element in the range `[stencil, stencil + n)`
 * @param output_begin Beginning of the sequence of matched payloads retrieved for each key
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <int CGSize,
          int BlockSize,
          typename InputIt,
          typename StencilIt,
          typename Predicate,
          typename OutputIt,
          typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void find_if_n(InputIt first,
                                                        cuco::detail::index_type n,
                                                        StencilIt stencil,
                                                        Predicate pred,
                                                        OutputIt output_begin,
                                                        Ref ref)
{
  namespace cg = cooperative_groups;

  auto const block       = cg::this_thread_block();
  auto const thread_idx  = block.thread_rank();
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  using output_type = typename find_buffer<Ref>::type;
  __shared__ output_type output_buffer[BlockSize / CGSize];

  auto constexpr has_payload =
    not cuda::std::is_same_v<typename Ref::key_type, typename Ref::value_type>;

  auto const sentinel = [&]() {
    if constexpr (has_payload) {
      return ref.empty_value_sentinel();
    } else {
      return ref.empty_key_sentinel();
    }
  }();

  auto output = cuda::proclaim_return_type<output_type>([&] __device__(auto found) {
    if constexpr (has_payload) {
      return found == ref.end() ? sentinel : found->second;
    } else {
      return found == ref.end() ? sentinel : *found;
    }
  });

  while ((idx - thread_idx / CGSize) < n) {  // the whole thread block falls into the same iteration
    if constexpr (CGSize == 1) {
      if (idx < n) {
        typename cuda::std::iterator_traits<InputIt>::value_type const key = *(first + idx);
        auto const found                                                   = ref.find(key);
        /*
         * The ld.relaxed.gpu instruction causes L1 to flush more frequently, causing increased
         * sector stores from L2 to global memory. By writing results to shared memory and then
         * synchronizing before writing back to global, we no longer rely on L1, preventing the
         * increase in sector stores from L2 to global and improving performance.
         */
        output_buffer[thread_idx] = pred(*(stencil + idx)) ? output(found) : sentinel;
      }
      block.sync();
      if (idx < n) { *(output_begin + idx) = output_buffer[thread_idx]; }
    } else {
      auto const tile = cg::tiled_partition<CGSize, cg::thread_block>(block);
      if (idx < n) {
        typename cuda::std::iterator_traits<InputIt>::value_type const key = *(first + idx);
        auto const found                                                   = ref.find(tile, key);

        if (tile.thread_rank() == 0) {
          *(output_begin + idx) = pred(*(stencil + idx)) ? output(found) : sentinel;
        }
      }
    }
    idx += loop_stride;
  }
}

/**
 * @brief Inserts all elements in the range `[first, last)`.
 *
 * @note: For a given element `*(first + i)`, if the container doesn't already contain an element
 * with an equivalent key, inserts the element at a location pointed by `iter` and writes
 * `iter` to `found_begin + i` and writes `true` to `inserted_begin + i`. Otherwise, finds the
 * location of the equivalent element, `iter` and writes `iter` to `found_begin + i` and writes
 * `false` to `inserted_begin + i`.
 *
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the `value_type` of the data structure
 * @tparam FoundIt Device accessible random access output iterator whose `value_type`
 * is constructible from `map::iterator` type
 * @tparam InsertedIt Device accessible random access output iterator whose `value_type`
 * is constructible from `bool`
 * @tparam Ref Type of non-owning device container ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param found_begin Beginning of the sequence of elements found for each key
 * @param inserted_begin Beginning of the sequence of booleans for the presence of each key
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <int CGSize,
          int BlockSize,
          typename InputIt,
          typename FoundIt,
          typename InsertedIt,
          typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void insert_and_find(InputIt first,
                                                              cuco::detail::index_type n,
                                                              FoundIt found_begin,
                                                              InsertedIt inserted_begin,
                                                              Ref ref)
{
  namespace cg = cooperative_groups;

  auto const block       = cg::this_thread_block();
  auto const thread_idx  = block.thread_rank();
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  using output_type = typename find_buffer<Ref>::type;

  auto constexpr has_payload =
    not cuda::std::is_same_v<typename Ref::key_type, typename Ref::value_type>;

  auto output = cuda::proclaim_return_type<output_type>([&] __device__(auto found) {
    if constexpr (has_payload) {
      return found->second;
    } else {
      return *found;
    }
  });

  __shared__ output_type output_location_buffer[BlockSize / CGSize];
  __shared__ bool output_inserted_buffer[BlockSize / CGSize];

  while ((idx - thread_idx / CGSize) < n) {  // the whole thread block falls into the same iteration
    if constexpr (CGSize == 1) {
      if (idx < n) {
        typename cuda::std::iterator_traits<InputIt>::value_type const insert_element{
          *(first + idx)};
        auto const [iter, inserted] = ref.insert_and_find(insert_element);
        /*
         * The ld.relaxed.gpu instruction causes L1 to flush more frequently, causing increased
         * sector stores from L2 to global memory. By writing results to shared memory and then
         * synchronizing before writing back to global, we no longer rely on L1, preventing the
         * increase in sector stores from L2 to global and improving performance.
         */
        output_location_buffer[thread_idx] = output(iter);
        output_inserted_buffer[thread_idx] = inserted;
      }
      block.sync();
      if (idx < n) {
        *(found_begin + idx)    = output_location_buffer[thread_idx];
        *(inserted_begin + idx) = output_inserted_buffer[thread_idx];
      }
    } else {
      auto const tile = cg::tiled_partition<CGSize, cg::thread_block>(cg::this_thread_block());
      if (idx < n) {
        typename cuda::std::iterator_traits<InputIt>::value_type const insert_element{
          *(first + idx)};
        auto const [iter, inserted] = ref.insert_and_find(tile, insert_element);
        if (tile.thread_rank() == 0) {
          *(found_begin + idx)    = output(iter);
          *(inserted_begin + idx) = inserted;
        }
      }
    }
    idx += loop_stride;
  }
}

/**
 * @brief Counts the occurrences of keys in `[first, last)` contained in the container
 *
 * @tparam IsOuter Flag indicating whether it's an outer count or not
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam InputIt Device accessible input iterator
 * @tparam AtomicT Atomic counter type
 * @tparam Ref Type of non-owning device container ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param count Number of matches
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <bool IsOuter, int CGSize, int BlockSize, typename InputIt, typename AtomicT, typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void count(InputIt first,
                                                    cuco::detail::index_type n,
                                                    AtomicT* count,
                                                    Ref ref)
{
  using size_type = typename Ref::size_type;

  size_type constexpr outer_min_count = 1;

  using BlockReduce = cub::BlockReduce<size_type, BlockSize>;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  size_type thread_count = 0;

  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  while (idx < n) {
    typename cuda::std::iterator_traits<InputIt>::value_type const key = *(first + idx);
    if constexpr (CGSize == 1) {
      if constexpr (IsOuter) {
        thread_count += max(ref.count(key), outer_min_count);
      } else {
        thread_count += ref.count(key);
      }
    } else {
      auto const tile =
        cooperative_groups::tiled_partition<CGSize, cooperative_groups::thread_block>(
          cooperative_groups::this_thread_block());
      if constexpr (IsOuter) {
        auto temp_count = ref.count(tile, key);
        if (tile.all(temp_count == 0) and tile.thread_rank() == 0) { ++temp_count; }
        thread_count += temp_count;
      } else {
        thread_count += ref.count(tile, key);
      }
    }
    idx += loop_stride;
  }

  auto const block_count = BlockReduce(temp_storage).Sum(thread_count);
  if (threadIdx.x == 0) { count->fetch_add(block_count, cuda::std::memory_order_relaxed); }
}

/**
 * @brief Counts the occurrences of each key in `[first, last)` contained in the container
 * and stores the counts in the output array.
 *
 * @tparam IsOuter Flag indicating whether it's an outer count or not
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam InputIt Device accessible input iterator
 * @tparam OutputIt Device accessible output iterator
 * @tparam Ref Type of non-owning device container ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param output_begin Beginning of the sequence of counts for each key
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <bool IsOuter,
          int CGSize,
          int BlockSize,
          typename InputIt,
          typename OutputIt,
          typename Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void count_each(InputIt first,
                                                         cuco::detail::index_type n,
                                                         OutputIt output_begin,
                                                         Ref ref)
{
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  using size_type                     = typename Ref::size_type;
  size_type constexpr outer_min_count = 1;

  while (idx < n) {
    typename cuda::std::iterator_traits<InputIt>::value_type const key = *(first + idx);
    if constexpr (CGSize == 1) {
      if constexpr (IsOuter) {
        *(output_begin + idx) = max(ref.count(key), size_type{outer_min_count});
      } else {
        *(output_begin + idx) = ref.count(key);
      }
    } else {
      auto const tile =
        cooperative_groups::tiled_partition<CGSize, cooperative_groups::thread_block>(
          cooperative_groups::this_thread_block());
      if constexpr (IsOuter) {
        auto temp_count = ref.count(tile, key);
        if (tile.all(temp_count == 0) and tile.thread_rank() == 0) { ++temp_count; }
        auto const count =
          cooperative_groups::reduce(tile, temp_count, cooperative_groups::plus<size_type>());
        if (tile.thread_rank() == 0) { *(output_begin + idx) = count; }
      } else {
        auto const count = cooperative_groups::reduce(
          tile, ref.count(tile, key), cooperative_groups::plus<size_type>());
        if (tile.thread_rank() == 0) { *(output_begin + idx) = count; }
      }
    }
    idx += loop_stride;
  }
}

/**
 * @brief Retrieves the equivalent container elements of all keys in the range `[input_probe,
 * input_probe + n)`.
 *
 * If key `k = *(input_probe + i)` has one or more matches in the container, copies `k` to
 * `output_probe` and associated slot contents to `output_match`, respectively. The output order is
 * unspecified.
 *
 * @tparam IsOuter Flag indicating whether it's an outer count or not
 * @tparam block_size The size of the thread block
 * @tparam InputProbeIt Device accessible input iterator
 * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
 * convertible to the `InputProbeIt`'s `value_type`
 * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
 * convertible to the container's `value_type`
 * @tparam AtomicCounter Integral atomic type that follows the same semantics as
 * `cuda::(std::)atomic(_ref)`
 * @tparam Ref Type of non-owning device ref allowing access to storage
 *
 * @param input_probe Beginning of the sequence of input keys
 * @param n Number of the keys to query
 * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
 * `output_match`
 * @param output_match Beginning of the sequence of matching elements
 * @param atomic_counter Pointer to an atomic object of integral type that is used to count the
 * number of output elements
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <bool IsOuter,
          int BlockSize,
          class InputProbeIt,
          class OutputProbeIt,
          class OutputMatchIt,
          class AtomicCounter,
          class Ref>
CUCO_KERNEL void retrieve(InputProbeIt input_probe,
                          cuco::detail::index_type n,
                          OutputProbeIt output_probe,
                          OutputMatchIt output_match,
                          AtomicCounter* atomic_counter,
                          Ref ref)
{
  namespace cg = cooperative_groups;

  auto const block              = cg::this_thread_block();
  auto constexpr tiles_in_block = BlockSize / Ref::cg_size;

  auto const block_begin_offset = block.group_index().x * tiles_in_block;
  auto const block_end_offset =
    min(n, static_cast<cuco::detail::index_type>(block_begin_offset + tiles_in_block));

  if (block_begin_offset < block_end_offset) {
    if constexpr (IsOuter) {
      ref.template retrieve_outer<BlockSize>(block,
                                             input_probe + block_begin_offset,
                                             input_probe + block_end_offset,
                                             output_probe,
                                             output_match,
                                             *atomic_counter);
    } else {
      ref.template retrieve<BlockSize>(block,
                                       input_probe + block_begin_offset,
                                       input_probe + block_end_offset,
                                       output_probe,
                                       output_match,
                                       *atomic_counter);
    }
  }
}

/**
 * @brief Calculates the number of filled slots for the given bucket storage.
 *
 * @tparam BlockSize Number of threads in each block
 * @tparam StorageRef Type of non-owning ref allowing access to storage
 * @tparam Predicate Type of predicate indicating if the given slot is filled
 * @tparam AtomicT Atomic counter type
 *
 * @param storage Non-owning device ref used to access the slot storage
 * @param is_filled Predicate indicating if the given slot is filled
 * @param count Number of filled slots
 */
template <int BlockSize, typename StorageRef, typename Predicate, typename AtomicT>
CUCO_KERNEL __launch_bounds__(BlockSize) void size(StorageRef storage,
                                                   Predicate is_filled,
                                                   AtomicT* count)
{
  using size_type = typename StorageRef::size_type;

  auto const loop_stride = cuco::detail::grid_stride();
  auto idx               = cuco::detail::global_thread_id();

  size_type thread_count = 0;
  auto const n           = storage.capacity();

  while (idx < n) {
    thread_count += static_cast<size_type>(is_filled(*(storage.data() + idx)));

    idx += loop_stride;
  }

  using BlockReduce = cub::BlockReduce<size_type, BlockSize>;
  __shared__ typename BlockReduce::TempStorage temp_storage;
  auto const block_count = BlockReduce(temp_storage).Sum(thread_count);
  if (threadIdx.x == 0) { count->fetch_add(block_count, cuda::std::memory_order_relaxed); }
}

template <int BlockSize, typename ContainerRef, typename Predicate>
CUCO_KERNEL __launch_bounds__(BlockSize) void rehash(
  typename ContainerRef::storage_ref_type storage_ref,
  ContainerRef container_ref,
  Predicate is_filled)
{
  namespace cg = cooperative_groups;

  __shared__ typename ContainerRef::value_type buffer[BlockSize];
  __shared__ unsigned int buffer_size;

  auto constexpr cg_size = ContainerRef::cg_size;
  auto const block       = cg::this_thread_block();
  auto const tile        = cg::tiled_partition<cg_size, cg::thread_block>(block);

  auto const thread_rank         = block.thread_rank();
  auto constexpr tiles_per_block = BlockSize / cg_size;  // tile.meta_group_size() but constexpr
  auto const tile_rank           = tile.meta_group_rank();
  auto const loop_stride         = cuco::detail::grid_stride();
  auto idx                       = cuco::detail::global_thread_id();
  auto const n                   = storage_ref.num_buckets();

  while (idx - thread_rank < n) {
    if (thread_rank == 0) { buffer_size = 0; }
    block.sync();

    // gather values in shmem buffer
    if (idx < n) {
      auto const bucket = storage_ref[idx];

      for (auto const& slot : bucket) {
        if (is_filled(slot)) { buffer[atomicAdd_block(&buffer_size, 1)] = slot; }
      }
    }
    block.sync();

    auto const local_buffer_size = buffer_size;

    // insert from shmem buffer into the container
    for (auto tidx = tile_rank; tidx < local_buffer_size; tidx += tiles_per_block) {
      container_ref.insert(tile, buffer[tidx]);
    }
    block.sync();

    idx += loop_stride;
  }
}
}  // namespace cuco::detail::open_addressing_ns
