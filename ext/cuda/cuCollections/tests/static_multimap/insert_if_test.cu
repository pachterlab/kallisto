/*
 * Copyright (c) 2021-2024, NVIDIA CORPORATION.
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

#include <cuco/static_multimap.cuh>

#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

#include <catch2/catch_template_test_macros.hpp>

template <typename Map>
void test_insert_if(Map& map, std::size_t size)
{
  using Key   = typename Map::key_type;
  using Value = typename Map::mapped_type;

  // 50% insertion
  auto const pred       = [] __device__(Key k) { return k % 2 == 0; };
  auto const keys_begin = thrust::counting_iterator<Key>{0};

  SECTION("Count of n / 2 insertions should be n / 2.")
  {
    auto const pairs_begin = thrust::make_transform_iterator(
      keys_begin, cuda::proclaim_return_type<cuco::pair<Key, Value>>([] __device__(auto i) {
        return cuco::pair<Key, Value>{i, i};
      }));

    auto const num = map.insert_if(pairs_begin, pairs_begin + size, keys_begin, pred);
    REQUIRE(num * 2 == size);

    auto const count = map.count(keys_begin, keys_begin + size);
    REQUIRE(count * 2 == size);
  }

  SECTION("Inserting the same element n / 2 times should return n / 2.")
  {
    auto const pairs_begin = thrust::constant_iterator<cuco::pair<Key, Value>>{{1, 1}};

    auto const num = map.insert_if(pairs_begin, pairs_begin + size, keys_begin, pred);
    REQUIRE(num * 2 == size);

    auto const count = map.count(keys_begin, keys_begin + size);
    REQUIRE(count * 2 == size);
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_multimap insert_if",
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
  constexpr std::size_t num_keys{1'000};

  using extent_type = cuco::extent<std::size_t>;
  using probe       = std::conditional_t<
          Probe == cuco::test::probe_sequence::linear_probing,
          cuco::linear_probing<CGSize, cuco::murmurhash3_32<Key>>,
          cuco::double_hashing<CGSize, cuco::murmurhash3_32<Key>, cuco::murmurhash3_32<Key>>>;

  auto map = cuco::experimental::static_multimap<Key,
                                                 Value,
                                                 extent_type,
                                                 cuda::thread_scope_device,
                                                 cuda::std::equal_to<Key>,
                                                 probe,
                                                 cuco::cuda_allocator<cuda::std::byte>,
                                                 cuco::storage<2>>{
    num_keys * 2, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{-1}};

  test_insert_if(map, num_keys);
}
