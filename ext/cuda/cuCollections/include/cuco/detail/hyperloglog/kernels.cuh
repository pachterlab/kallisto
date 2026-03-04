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

#include <cuco/detail/utility/cuda.cuh>
#include <cuco/utility/cuda_thread_scope.cuh>

#include <cuda/std/array>
#include <cuda/std/span>

#include <cooperative_groups.h>

#include <cstddef>

namespace cuco::hyperloglog_ns::detail {
CUCO_SUPPRESS_KERNEL_WARNINGS

template <class RefType>
CUCO_KERNEL void clear(RefType ref)
{
  auto const block = cooperative_groups::this_thread_block();
  if (block.group_index().x == 0) { ref.clear(block); }
}

template <int32_t VectorSize, class RefType>
CUCO_KERNEL void add_shmem_vectorized(typename RefType::value_type const* first,
                                      cuco::detail::index_type n,
                                      RefType ref)
{
  using value_type     = typename RefType::value_type;
  using vector_type    = cuda::std::array<value_type, VectorSize>;
  using local_ref_type = typename RefType::template with_scope<cuda::thread_scope_block>;

  // Base address of dynamic shared memory is guaranteed to be aligned to at least 16 bytes which is
  // sufficient for this purpose
  extern __shared__ cuda::std::byte local_sketch[];

  auto const loop_stride = cuco::detail::grid_stride();
  auto idx               = cuco::detail::global_thread_id();
  auto const grid        = cooperative_groups::this_grid();
  auto const block       = cooperative_groups::this_thread_block();

  local_ref_type local_ref(cuda::std::span{local_sketch, ref.sketch_bytes()}, {});
  local_ref.clear(block);
  block.sync();

  // each thread processes VectorSize-many items per iteration
  vector_type vec;
  while (idx < n / VectorSize) {
    vec = *reinterpret_cast<vector_type*>(
      __builtin_assume_aligned(first + idx * VectorSize, sizeof(vector_type)));
    for (auto const& i : vec) {
      local_ref.add(i);
    }
    idx += loop_stride;
  }
  // a single thread processes the remaining items
#if defined(CUCO_HAS_CG_INVOKE_ONE)
  cooperative_groups::invoke_one(grid, [&]() {
    auto const remainder = n % VectorSize;
    for (int i = 0; i < remainder; ++i) {
      local_ref.add(*(first + n - i - 1));
    }
  });
#else
  if (grid.thread_rank() == 0) {
    auto const remainder = n % VectorSize;
    for (int i = 0; i < remainder; ++i) {
      local_ref.add(*(first + n - i - 1));
    }
  }
#endif
  block.sync();

  ref.merge(block, local_ref);
}

template <class InputIt, class RefType>
CUCO_KERNEL void add_shmem(InputIt first, cuco::detail::index_type n, RefType ref)
{
  using local_ref_type = typename RefType::template with_scope<cuda::thread_scope_block>;

  // TODO assert alignment
  extern __shared__ cuda::std::byte local_sketch[];

  auto const loop_stride = cuco::detail::grid_stride();
  auto idx               = cuco::detail::global_thread_id();
  auto const block       = cooperative_groups::this_thread_block();

  local_ref_type local_ref(cuda::std::span{local_sketch, ref.sketch_bytes()}, {});
  local_ref.clear(block);
  block.sync();

  while (idx < n) {
    local_ref.add(*(first + idx));
    idx += loop_stride;
  }
  block.sync();

  ref.merge(block, local_ref);
}

template <class InputIt, class RefType>
CUCO_KERNEL void add_gmem(InputIt first, cuco::detail::index_type n, RefType ref)
{
  auto const loop_stride = cuco::detail::grid_stride();
  auto idx               = cuco::detail::global_thread_id();

  while (idx < n) {
    ref.add(*(first + idx));
    idx += loop_stride;
  }
}

template <class OtherRefType, class RefType>
CUCO_KERNEL void merge(OtherRefType other_ref, RefType ref)
{
  auto const block = cooperative_groups::this_thread_block();
  if (block.group_index().x == 0) { ref.merge(block, other_ref); }
}

// TODO this kernel currently isn't being used
template <class RefType>
CUCO_KERNEL void estimate(std::size_t* cardinality, RefType ref)
{
  auto const block = cooperative_groups::this_thread_block();
  if (block.group_index().x == 0) {
    auto const estimate = ref.estimate(block);
    if (block.thread_rank() == 0) { *cardinality = estimate; }
  }
}
}  // namespace cuco::hyperloglog_ns::detail
