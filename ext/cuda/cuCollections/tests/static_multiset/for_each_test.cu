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

#include <cuco/detail/utility/cuda.hpp>
#include <cuco/static_multiset.cuh>

#include <cuda/atomic>
#include <cuda/functional>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>

#include <catch2/catch_template_test_macros.hpp>

#include <cstddef>

template <class Ref, class InputIt, class AtomicErrorCounter>
CUCO_KERNEL void for_each_check_scalar(Ref ref,
                                       InputIt first,
                                       std::size_t n,
                                       std::size_t multiplicity,
                                       AtomicErrorCounter* error_counter)
{
  static_assert(Ref::cg_size == 1, "Scalar test must have cg_size==1");
  auto const loop_stride = cuco::detail::grid_stride();
  auto idx               = cuco::detail::global_thread_id();

  while (idx < n) {
    auto const& key     = *(first + idx);
    std::size_t matches = 0;
    ref.for_each(key, [&] __device__(auto const slot) {
      if (ref.key_eq()(key, slot)) { matches++; }
    });
    if (matches != multiplicity) { error_counter->fetch_add(1, cuda::memory_order_relaxed); }
    idx += loop_stride;
  }
}

template <bool Synced, class Ref, class InputIt, class AtomicErrorCounter>
CUCO_KERNEL void for_each_check_cooperative(Ref ref,
                                            InputIt first,
                                            std::size_t n,
                                            std::size_t multiplicity,
                                            AtomicErrorCounter* error_counter)
{
  auto const loop_stride = cuco::detail::grid_stride() / Ref::cg_size;
  auto idx               = cuco::detail::global_thread_id() / Ref::cg_size;
  ;

  while (idx < n) {
    auto const tile =
      cooperative_groups::tiled_partition<Ref::cg_size, cooperative_groups::thread_block>(
        cooperative_groups::this_thread_block());
    auto const& key            = *(first + idx);
    std::size_t thread_matches = 0;
    if constexpr (Synced) {
      ref.for_each(
        tile,
        key,
        [&] __device__(auto const slot) {
          if (ref.key_eq()(key, slot)) { thread_matches++; }
        },
        [] __device__(auto group) { group.sync(); });
    } else {
      ref.for_each(tile, key, [&] __device__(auto const slot) {
        if (ref.key_eq()(key, slot)) { thread_matches++; }
      });
    }
    auto const tile_matches =
      cooperative_groups::reduce(tile, thread_matches, cooperative_groups::plus<std::size_t>());
    if (tile_matches != multiplicity and tile.thread_rank() == 0) {
      error_counter->fetch_add(1, cuda::memory_order_relaxed);
    }
    idx += loop_stride;
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_multiset for_each tests",
  "",
  ((typename Key, cuco::test::probe_sequence Probe, int CGSize), Key, Probe, CGSize),
  (int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, cuco::test::probe_sequence::linear_probing, 2))
{
  constexpr size_t num_unique_keys{400};
  constexpr size_t key_multiplicity{5};
  constexpr size_t num_keys{num_unique_keys * key_multiplicity};

  using probe = std::conditional_t<Probe == cuco::test::probe_sequence::linear_probing,
                                   cuco::linear_probing<CGSize, cuco::default_hash_function<Key>>,
                                   cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>>;

  auto set =
    cuco::static_multiset{num_keys, cuco::empty_key<Key>{-1}, {}, probe{}, {}, cuco::storage<2>{}};

  auto unique_keys_begin  = thrust::counting_iterator<Key>(0);
  auto gen_duplicate_keys = cuda::proclaim_return_type<Key>(
    [] __device__(auto const& k) { return static_cast<Key>(k % num_unique_keys); });
  auto keys_begin = thrust::make_transform_iterator(unique_keys_begin, gen_duplicate_keys);

  set.insert(keys_begin, keys_begin + num_keys);

  using error_counter_type = cuda::atomic<std::size_t, cuda::thread_scope_system>;
  error_counter_type* error_counter;
  CUCO_CUDA_TRY(cudaMallocHost(&error_counter, sizeof(error_counter_type)));
  new (error_counter) error_counter_type{0};

  auto const grid_size  = cuco::detail::grid_size(num_unique_keys, CGSize);
  auto const block_size = cuco::detail::default_block_size();

  // test scalar for_each
  if constexpr (CGSize == 1) {
    for_each_check_scalar<<<grid_size, block_size>>>(
      set.ref(cuco::for_each), unique_keys_begin, num_unique_keys, key_multiplicity, error_counter);
    CUCO_CUDA_TRY(cudaDeviceSynchronize());
    REQUIRE(error_counter->load() == 0);
    error_counter->store(0);
  }

  // test CG for_each
  for_each_check_cooperative<false><<<grid_size, block_size>>>(
    set.ref(cuco::for_each), unique_keys_begin, num_unique_keys, key_multiplicity, error_counter);
  CUCO_CUDA_TRY(cudaDeviceSynchronize());
  REQUIRE(error_counter->load() == 0);
  error_counter->store(0);

  // test synchronized CG for_each
  for_each_check_cooperative<true><<<grid_size, block_size>>>(
    set.ref(cuco::for_each), unique_keys_begin, num_unique_keys, key_multiplicity, error_counter);
  CUCO_CUDA_TRY(cudaDeviceSynchronize());
  REQUIRE(error_counter->load() == 0);

  CUCO_CUDA_TRY(cudaFreeHost(error_counter));
}