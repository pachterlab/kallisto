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
#include <cuco/hyperloglog.cuh>

#include <cuda/std/cstddef>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>

#include <iostream>

/**
 * @file device_ref_example.cu
 * @brief Demonstrates usage of `cuco::hyperloglog` device-side APIs.
 *
 * This example demonstrates how the non-owning reference type `cuco::hyperloglog_ref`
 * can be used to implement a custom kernel that fuses the cardinality estimation step with any
 * other workload that traverses the input data.
 */

template <class RefType, class InputIt>
__global__ void fused_kernel(RefType ref, InputIt first, std::size_t n)
{
  // Transform the reference type (with device scope) to a reference type with block scope
  using local_ref_type = typename RefType::template with_scope<cuda::thread_scope_block>;

  // Shared memory storage for the block-local estimator
  extern __shared__ cuda::std::byte local_sketch[];

  // The following check is optional since the base address of dynamic shared memory is guaranteed
  // to meet the alignment requirements
  /*
  auto const alignment =
    1ull << cuda::std::countr_zero(reinterpret_cast<std::uintptr_t>(local_sketch));
  assert(alignment >= local_ref_type::sketch_alignment());
  */

  auto const loop_stride = gridDim.x * blockDim.x;
  auto idx               = blockDim.x * blockIdx.x + threadIdx.x;
  auto const block       = cooperative_groups::this_thread_block();

  // Create the local estimator with the shared memory storage
  local_ref_type local_ref(cuda::std::span{local_sketch, ref.sketch_bytes()});

  // Initialize the local estimator
  local_ref.clear(block);
  block.sync();

  while (idx < n) {
    auto const& item = *(first + idx);

    // Add each item to the local estimator
    local_ref.add(item);

    /*
    Here we can add some custom workload that takes the input `item`.

    The idea is that cardinality estimation can be fused with any other workload that
    traverses the data. Since `local_ref.add` can run close to the SOL of the DRAM bandwidth, we get
    the estimate "for free" while performing other computations over the data.
    */

    idx += loop_stride;
  }
  block.sync();

  // We can also compute the local estimate on the device
  // auto const local_estimate = local_ref.estimate(block);
  if (block.thread_rank() == 0) {
    // The local estimate should approximately be `num_items`/`gridDim.x`
    // printf("Estimate for block %d = %llu\n", blockIdx.x, local_estimate);
  }

  // In the end, we merge the shared memory estimator into the global estimator which gives us the
  // final result
  ref.merge(block, local_ref);
}

template <typename Ref, typename InputIt, typename OutputIt>
__global__ void device_estimate_kernel(cuco::sketch_size_kb sketch_size_kb,
                                       InputIt in,
                                       size_t n,
                                       OutputIt out)
{
  extern __shared__ cuda::std::byte local_sketch[];

  auto const block = cooperative_groups::this_thread_block();

  // only a single block computes the estimate
  if (block.group_index().x == 0) {
    Ref estimator(cuda::std::span(local_sketch, Ref::sketch_bytes(sketch_size_kb)));

    estimator.clear(block);
    block.sync();

    for (int i = block.thread_rank(); i < n; i += block.num_threads()) {
      estimator.add(*(in + i));
    }
    block.sync();
    // we can compute the final estimate on the device and return the result to the host
    auto const estimate = estimator.estimate(block);

    if (block.thread_rank() == 0) { *out = estimate; }
  }
}

int main(void)
{
  using T                         = int;
  using estimator_type            = cuco::hyperloglog<T>;
  constexpr std::size_t num_items = 1ull << 28;  // 1GB
  auto const sketch_size_kb       = 32_KB;

  thrust::device_vector<T> items(num_items);

  // Generate `num_items` distinct items
  thrust::sequence(items.begin(), items.end(), 0);

  // Initialize the estimator
  estimator_type estimator(sketch_size_kb);

  // Add all items to the estimator
  estimator.add(items.begin(), items.end());

  // Calculate the cardinality estimate from the bulk operation
  std::size_t const estimated_cardinality_bulk = estimator.estimate();

  // Clear the estimator so it can be reused
  estimator.clear();

  // Number of dynamic shared memory bytes required to store a CTA-local sketch
  auto const sketch_bytes = estimator.sketch_bytes();

  // Call the custom kernel and pass a non-owning reference to the estimator to the GPU
  fused_kernel<<<10, 512, sketch_bytes>>>(estimator.ref(), items.begin(), num_items);

  // Calculate the cardinality estimate from the custom kernel
  std::size_t const estimated_cardinality_custom = estimator.estimate();

  thrust::device_vector<std::size_t> device_estimate(1);
  device_estimate_kernel<typename estimator_type::ref_type<cuda::thread_scope_block>>
    <<<1, 512, sketch_bytes>>>(sketch_size_kb, items.begin(), num_items, device_estimate.begin());

  std::size_t const estimated_cardinality_device = device_estimate[0];

  if (estimated_cardinality_custom == estimated_cardinality_bulk and
      estimated_cardinality_device == estimated_cardinality_bulk) {
    std::cout << "Success! Cardinality estimates are identical" << std::endl;
  }

  return 0;
}
