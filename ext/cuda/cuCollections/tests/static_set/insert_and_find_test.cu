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

#include <test_utils.hpp>

#include <cuco/static_set.cuh>

#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>

#include <catch2/catch_template_test_macros.hpp>

template <typename Set>
void test_insert_and_find(Set& set, std::size_t num_keys)
{
  using Key                     = typename Set::key_type;
  static auto constexpr cg_size = Set::cg_size;

  auto const keys_begin = thrust::counting_iterator<Key>(0);
  auto const keys_end   = thrust::counting_iterator<Key>(num_keys);

  thrust::device_vector<Key> iters1(num_keys);
  thrust::device_vector<int> iters2(num_keys);

  thrust::device_vector<bool> inserted(num_keys);

  // insert first time, fills inserted with true
  set.insert_and_find(keys_begin, keys_end, iters1.begin(), inserted.begin());
  REQUIRE(cuco::test::all_of(inserted.begin(), inserted.end(), cuda::std::identity{}));

  // insert second time, fills inserted with false as keys already in set
  set.insert_and_find(keys_begin, keys_end, iters2.begin(), inserted.begin());
  REQUIRE(cuco::test::none_of(inserted.begin(), inserted.end(), cuda::std::identity{}));

  // both iters1 and iters2 should be same, as keys will be referring to same slot
  REQUIRE(
    cuco::test::equal(iters1.begin(), iters1.end(), iters2.begin(), cuda::std::equal_to<Key>{}));
}

TEMPLATE_TEST_CASE_SIG(
  "static_set Insert and find",
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
  constexpr std::size_t num_keys{400};

  using probe = std::conditional_t<Probe == cuco::test::probe_sequence::linear_probing,
                                   cuco::linear_probing<CGSize, cuco::default_hash_function<Key>>,
                                   cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>>;

  auto set =
    cuco::static_set{num_keys, cuco::empty_key<Key>{-1}, {}, probe{}, {}, cuco::storage<2>{}};

  test_insert_and_find(set, num_keys);
}
