/*
 * Copyright (c) 2022, Jonas Hahnfeld, CERN.
 * Copyright (c) 2022-2025, NVIDIA CORPORATION.
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

#include <cuco/static_map.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = std::size_t;

TEMPLATE_TEST_CASE_SIG(
  "static_map insert_and_find tests",
  "",
  ((typename Key, typename Value, cuco::test::probe_sequence Probe, int CGSize),
   Key,
   Value,
   Probe,
   CGSize),
  (int32_t, int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int32_t, int64_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, int64_t, cuco::test::probe_sequence::linear_probing, 2))
{
#if !defined(CUCO_HAS_INDEPENDENT_THREADS)
  if constexpr (cuco::detail::is_packable<cuco::pair<Key, Value>>())
#endif
  {
    using probe = std::conditional_t<
      Probe == cuco::test::probe_sequence::linear_probing,
      cuco::linear_probing<CGSize, cuco::murmurhash3_32<Key>>,
      cuco::double_hashing<CGSize, cuco::murmurhash3_32<Key>, cuco::murmurhash3_32<Key>>>;

    constexpr size_type num_keys{400};

    auto map = cuco::static_map<Key,
                                Value,
                                cuco::extent<size_type>,
                                cuda::thread_scope_device,
                                cuda::std::equal_to<Key>,
                                probe,
                                cuco::cuda_allocator<cuda::std::byte>,
                                cuco::storage<2>>{
      num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{-1}};

    auto pairs_begin = thrust::make_transform_iterator(
      thrust::counting_iterator<size_type>(0),
      cuda::proclaim_return_type<cuco::pair<Key, Value>>(
        [] __device__(auto i) { return cuco::pair<Key, Value>{i, 1}; }));

    thrust::device_vector<size_type> found1(num_keys);
    thrust::device_vector<size_type> found2(num_keys);

    thrust::device_vector<bool> inserted(num_keys);

    // insert first time, fills inserted with true
    map.insert_and_find(pairs_begin, pairs_begin + num_keys, found1.begin(), inserted.begin());
    REQUIRE(cuco::test::all_of(inserted.begin(), inserted.end(), cuda::std::identity{}));

    // insert second time, fills inserted with false as keys already in map
    map.insert_and_find(pairs_begin, pairs_begin + num_keys, found2.begin(), inserted.begin());
    REQUIRE(cuco::test::none_of(inserted.begin(), inserted.end(), cuda::std::identity{}));

    // both found1 and found2 should be same, as keys will be referring to same slot
    REQUIRE(
      cuco::test::equal(found1.begin(), found1.end(), found2.begin(), cuda::std::equal_to<Key>{}));
  }
}
