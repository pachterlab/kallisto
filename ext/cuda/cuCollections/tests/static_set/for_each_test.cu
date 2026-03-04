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

#include <test_utils.hpp>

#include <cuco/static_set.cuh>

#include <cuda/atomic>
#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = std::size_t;

template <typename Set>
void test_for_each(Set& set, size_type num_keys)
{
  using Key = typename Set::key_type;

  REQUIRE(num_keys % 2 == 0);

  cuda::stream_ref stream{cudaStream_t{nullptr}};

  // Insert keys
  auto keys_begin = thrust::make_transform_iterator(
    thrust::counting_iterator<size_type>{0}, cuda::proclaim_return_type<Key>([] __device__(auto i) {
      // generates a sequence of 1, 2, 1, 2, ...
      return static_cast<Key>(i);
    }));
  set.insert(keys_begin, keys_begin + num_keys, stream);

  using Allocator = cuco::cuda_allocator<cuda::atomic<size_type, cuda::thread_scope_device>>;
  cuco::detail::counter_storage<size_type, cuda::thread_scope_device, Allocator> counter_storage(
    Allocator{}, stream);
  counter_storage.reset(stream);

  // count the sum of all even keys
  set.for_each(
    [counter = counter_storage.data()] __device__(auto const slot) {
      if (slot % 2 == 0) { counter->fetch_add(slot, cuda::memory_order_relaxed); }
    },
    stream);
  REQUIRE(counter_storage.load_to_host(stream) == 249'500);

  counter_storage.reset(stream);

  // count the sum of all odd keys
  set.for_each(
    thrust::counting_iterator<size_type>(0),
    thrust::counting_iterator<size_type>(2 * num_keys),  // test for false-positives
    [counter = counter_storage.data()] __device__(auto const slot) {
      if (!(slot % 2 == 0)) { counter->fetch_add(slot, cuda::memory_order_relaxed); }
    },
    stream);
  REQUIRE(counter_storage.load_to_host(stream) == 250'000);
}

TEMPLATE_TEST_CASE_SIG(
  "static_set for_each tests",
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
  constexpr size_type num_keys{1'000};
  using probe = std::conditional_t<
    Probe == cuco::test::probe_sequence::linear_probing,
    cuco::linear_probing<CGSize, cuco::murmurhash3_32<Key>>,
    cuco::double_hashing<CGSize, cuco::murmurhash3_32<Key>, cuco::murmurhash3_32<Key>>>;

  using set_t = cuco::static_set<Key,
                                 cuco::extent<size_type>,
                                 cuda::thread_scope_device,
                                 cuda::std::equal_to<Key>,
                                 probe,
                                 cuco::cuda_allocator<cuda::std::byte>,
                                 cuco::storage<2>>;

  auto set = set_t{num_keys, cuco::empty_key<Key>{-1}};
  test_for_each(set, num_keys);
}
