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

#include <cuco/hash_functions.cuh>
#include <cuco/static_map.cuh>

#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = std::size_t;

template <typename Key, typename Hash>
void test_hash_function()
{
  using Value = int64_t;

  constexpr size_type num_keys{400};

  auto map = cuco::static_map<Key,
                              Value,
                              cuco::extent<size_type>,
                              cuda::thread_scope_device,
                              cuda::std::equal_to<Key>,
                              cuco::linear_probing<1, Hash>,
                              cuco::cuda_allocator<cuda::std::byte>,
                              cuco::storage<2>>{
    num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{-1}};

  auto keys_begin = thrust::counting_iterator<Key>(1);

  auto pairs_begin = thrust::make_transform_iterator(
    keys_begin, cuda::proclaim_return_type<cuco::pair<Key, Value>>([] __device__(auto i) {
      return cuco::pair<Key, Value>(i, i);
    }));

  thrust::device_vector<bool> d_keys_exist(num_keys);

  map.insert(pairs_begin, pairs_begin + num_keys);

  REQUIRE(map.size() == num_keys);

  map.contains(keys_begin, keys_begin + num_keys, d_keys_exist.begin());

  REQUIRE(cuco::test::all_of(d_keys_exist.begin(), d_keys_exist.end(), cuda::std::identity{}));
}

TEMPLATE_TEST_CASE_SIG("static_map hash tests", "", ((typename Key)), (int32_t), (int64_t))
{
  test_hash_function<Key, cuco::murmurhash3_32<Key>>();
  test_hash_function<Key, cuco::murmurhash3_x64_128<Key>>();
  test_hash_function<Key, cuco::xxhash_32<Key>>();
  test_hash_function<Key, cuco::xxhash_64<Key>>();
}
