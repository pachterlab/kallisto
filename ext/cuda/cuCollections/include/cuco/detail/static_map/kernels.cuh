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

#include <cuco/detail/bitwise_compare.cuh>
#include <cuco/detail/utility/cuda.cuh>
#include <cuco/operator.hpp>
#include <cuco/pair.cuh>
#include <cuco/types.cuh>

#include <cub/block/block_reduce.cuh>
#include <cuda/atomic>
#include <cuda/std/cstdint>
#include <cuda/std/iterator>

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

namespace cuco::detail::static_map_ns {
CUCO_SUPPRESS_KERNEL_WARNINGS

// TODO user insert_or_assign internally
/**
 * @brief For any key-value pair `{k, v}` in the range `[first, first + n)`, if a key equivalent to
 * `k` already exists in the container, assigns `v` to the mapped_type corresponding to the key `k`.
 * If the key does not exist, inserts the pair as if by insert.
 *
 * @note If multiple elements in `[first, first + n)` compare equal, it is unspecified which element
 * is inserted.
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
CUCO_KERNEL __launch_bounds__(BlockSize) void insert_or_assign(InputIt first,
                                                               cuco::detail::index_type n,
                                                               Ref ref)
{
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  while (idx < n) {
    typename cuda::std::iterator_traits<InputIt>::value_type const insert_pair = *(first + idx);
    if constexpr (CGSize == 1) {
      ref.insert_or_assign(insert_pair);
    } else {
      auto const tile =
        cooperative_groups::tiled_partition<CGSize, cooperative_groups::thread_block>(
          cooperative_groups::this_thread_block());
      ref.insert_or_assign(tile, insert_pair);
    }
    idx += loop_stride;
  }
}

/**
 * @brief For any key-value pair `{k, v}` in the range `[first, first + n)`, if a key equivalent to
 * `k` already exists in the container, then binary operation is applied using `op` callable object
 * on the existing value at slot and the element to insert. If the key does not exist, inserts the
 * pair as if by insert.
 *
 * @note Callable object to perform binary operation should be able to invoke as
 * Op(cuda::atomic_ref<T,Scope>, T>)
 * @note If `HasInit` is `true` and if `init == empty_sentinel_value`, we directly
 * `apply` the `op` instead of atomic store and then waiting for the payload to get materalized.
 * This has potential speedups when insert strategy is not `packed_cas`.
 *
 * @tparam HasInit Boolean to dispatch based on init parameter
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the `value_type` of the data structure
 * @tparam Init Type of init value convertible to payload type
 * @tparam Op Callable type used to peform `apply` operation.
 * @tparam Ref Type of non-owning device ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param init The init value of the op
 * @param op Callable object to perform apply operation.
 * @param ref Non-owning container device ref used to access the slot storage
 */
template <bool HasInit,
          int CGSize,
          int BlockSize,
          typename InputIt,
          typename Init,
          typename Op,
          typename Ref>
__global__ void insert_or_apply(
  InputIt first, cuco::detail::index_type n, [[maybe_unused]] Init init, Op op, Ref ref)
{
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  while (idx < n) {
    using value_type             = typename cuda::std::iterator_traits<InputIt>::value_type;
    value_type const insert_pair = *(first + idx);
    if constexpr (CGSize == 1) {
      if constexpr (HasInit) {
        ref.insert_or_apply(insert_pair, init, op);
      } else {
        ref.insert_or_apply(insert_pair, op);
      }
    } else {
      auto const tile =
        cooperative_groups::tiled_partition<CGSize, cooperative_groups::thread_block>(
          cooperative_groups::this_thread_block());
      if constexpr (HasInit) {
        ref.insert_or_apply(tile, insert_pair, init, op);
      } else {
        ref.insert_or_apply(tile, insert_pair, op);
      }
    }
    idx += loop_stride;
  }
}

/**
 * @brief For any key-value pair `{k, v}` in the range `[first, first + n)`, if a key equivalent to
 * `k` already exists in the container, then binary operation is applied using `op` callable object
 * on the existing value at slot and the element to insert. If the key does not exist, inserts the
 * pair as if by insert.
 *
 * @note Callable object to perform binary operation should be able to invoke as
 * Op(cuda::atomic_ref<T,Scope>, T>)
 * @note If `HasInit` is `true` and if `init == empty_sentinel_value`, we directly
 * `apply` the `op` instead of atomic store and then waiting for the payload to get materalized.
 * This has potential speedups when insert strategy is not `packed_cas`.
 *
 * @tparam HasInit Boolean to dispatch based on init parameter
 * @tparam CGSize Number of threads in each CG
 * @tparam BlockSize Number of threads in each block
 * @tparam SharedMapRefType The Shared Memory Map Ref Type
 * @tparam InputIt Device accessible input iterator whose `value_type` is
 * convertible to the `value_type` of the data structure
 * @tparam Init Type of init value convertible to payload type
 * @tparam Op Callable type used to peform `apply` operation.
 * @tparam Ref Type of non-owning device ref allowing access to storage
 *
 * @param first Beginning of the sequence of input elements
 * @param n Number of input elements
 * @param init The init value of the op
 * @param op Callable object to perform apply operation.
 * @param ref Non-owning container device ref used to access the slot storage
 * @param bucket_extent Bucket Extent used for shared memory map slot storage
 */
template <bool HasInit,
          int CGSize,
          int BlockSize,
          class SharedMapRefType,
          class InputIt,
          class Init,
          class Op,
          class Ref>
CUCO_KERNEL __launch_bounds__(BlockSize) void insert_or_apply_shmem(
  InputIt first,
  cuco::detail::index_type n,
  [[maybe_unused]] Init init,
  Op op,
  Ref ref,
  typename SharedMapRefType::extent_type bucket_extent)
{
  static_assert(CGSize == 1, "use shared_memory kernel only if cg_size == 1");
  namespace cg     = cooperative_groups;
  using Key        = typename Ref::key_type;
  using Value      = typename Ref::mapped_type;
  using value_type = typename cuda::std::iterator_traits<InputIt>::value_type;

  auto const block       = cg::this_thread_block();
  auto const thread_idx  = block.thread_rank();
  auto const loop_stride = cuco::detail::grid_stride() / CGSize;
  auto idx               = cuco::detail::global_thread_id() / CGSize;

  auto warp                  = cg::tiled_partition<32, cg::thread_block>(block);
  auto const warp_thread_idx = warp.thread_rank();

  // Shared map initialization
  __shared__ typename SharedMapRefType::value_type slots[bucket_extent.value()];
  using storage_ref_type = typename SharedMapRefType::storage_ref_type;
  auto storage           = storage_ref_type(bucket_extent, slots);
  auto const num_buckets = storage.num_buckets();

  using atomic_type = cuda::atomic<cuda::std::int32_t, cuda::thread_scope_block>;
  __shared__ atomic_type block_cardinality;
  if (thread_idx == 0) { new (&block_cardinality) atomic_type{}; }
  block.sync();

  auto shared_map     = SharedMapRefType{cuco::empty_key<Key>{ref.empty_key_sentinel()},
                                     cuco::empty_value<Value>{ref.empty_value_sentinel()},
                                     ref.key_eq(),
                                     ref.probing_scheme(),
                                         {},
                                     storage};
  auto shared_map_ref = shared_map.rebind_operators(cuco::op::insert_or_apply);
  shared_map_ref.initialize(block);
  block.sync();

  while ((idx - thread_idx / CGSize) < n) {
    cuda::std::int32_t inserted         = 0;
    cuda::std::int32_t warp_cardinality = 0;
    // insert-or-apply into the shared map first
    if (idx < n) {
      value_type const insert_pair = *(first + idx);
      if constexpr (HasInit) {
        inserted = shared_map_ref.insert_or_apply(insert_pair, init, op);
      } else {
        inserted = shared_map_ref.insert_or_apply(insert_pair, op);
      }
    }
    if (idx - warp_thread_idx < n) {  // all threads in warp particpate
      warp_cardinality = cg::reduce(warp, inserted, cg::plus<cuda::std::int32_t>());
    }
    if (warp_thread_idx == 0) {
      block_cardinality.fetch_add(warp_cardinality, cuda::memory_order_relaxed);
    }
    block.sync();
    if (block_cardinality > BlockSize) { break; }
    idx += loop_stride;
  }

  // insert-or-apply from shared map to global map
  auto bucket_idx = thread_idx;
  while (bucket_idx < num_buckets) {
    auto const slot = storage[bucket_idx][0];
    if (not cuco::detail::bitwise_compare(slot.first, ref.empty_key_sentinel())) {
      if constexpr (HasInit) {
        ref.insert_or_apply(slot, init, op);
      } else {
        ref.insert_or_apply(slot, op);
      }
    }
    bucket_idx += BlockSize;
  }

  // insert-or-apply into global map for the remaining elements whose block_cardinality
  // exceeds the cardinality threshold.
  if (block_cardinality > BlockSize) {
    idx += loop_stride;
    while (idx < n) {
      value_type const insert_pair = *(first + idx);
      if constexpr (HasInit) {
        ref.insert_or_apply(insert_pair, init, op);
      } else {
        ref.insert_or_apply(insert_pair, op);
      }
      idx += loop_stride;
    }
  }
}
}  // namespace cuco::detail::static_map_ns
