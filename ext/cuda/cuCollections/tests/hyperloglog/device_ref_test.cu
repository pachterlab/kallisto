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

#include <test_utils.hpp>

#include <cuco/hash_functions.cuh>
#include <cuco/hyperloglog.cuh>

#include <cuda/std/cstddef>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cmath>
#include <cstdint>

template <typename Ref, typename InputIt, typename OutputIt>
__global__ void estimate_kernel(cuco::sketch_size_kb sketch_size_kb,
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
    auto const estimate = estimator.estimate(block);
    if (block.thread_rank() == 0) { *out = estimate; }
  }
}

TEMPLATE_TEST_CASE_SIG("hyperloglog: device ref",
                       "",
                       ((typename T, typename Hash), T, Hash),
                       (int32_t, cuco::xxhash_64<int32_t>),
                       (int64_t, cuco::xxhash_64<int64_t>),
                       (__int128_t, cuco::xxhash_64<__int128_t>))
{
  using estimator_type = cuco::hyperloglog<T, cuda::thread_scope_device, Hash>;

  auto num_items_pow2 = GENERATE(25, 26, 28);
  auto hll_precision  = GENERATE(8, 10, 12, 13);
  auto sketch_size_kb = 4 * (1ull << hll_precision) / 1024;
  INFO("hll_precision=" << hll_precision);
  INFO("sketch_size_kb=" << sketch_size_kb);
  INFO("num_items=2^" << num_items_pow2);
  auto num_items = 1ull << num_items_pow2;

  thrust::device_vector<T> items(num_items);

  // Generate `num_items` distinct items
  thrust::sequence(items.begin(), items.end(), 0);

  // Initialize the estimator
  estimator_type estimator{cuco::sketch_size_kb(sketch_size_kb)};

  // Add all items to the estimator
  estimator.add(items.begin(), items.end());

  auto const host_estimate = estimator.estimate();

  thrust::device_vector<std::size_t> device_estimate(1);
  estimate_kernel<typename estimator_type::template ref_type<cuda::thread_scope_block>>
    <<<1, 512, estimator.sketch_bytes()>>>(
      cuco::sketch_size_kb(sketch_size_kb), items.begin(), num_items, device_estimate.begin());

  REQUIRE(device_estimate[0] == host_estimate);
}
