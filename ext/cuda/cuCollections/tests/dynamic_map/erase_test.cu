/*
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

#include <cuco/dynamic_map.cuh>

#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>

TEMPLATE_TEST_CASE_SIG("dynamic_map erase tests",
                       "",
                       ((typename Key, typename Value), Key, Value),
                       (int32_t, int32_t),
                       (int32_t, int64_t),
                       (int64_t, int32_t),
                       (int64_t, int64_t))
{
  constexpr std::size_t num_keys = 1'000'000;
  cuco::dynamic_map<Key, Value> map{num_keys * 2,
                                    cuco::empty_key<Key>{-1},
                                    cuco::empty_value<Value>{-1},
                                    cuco::erased_key<Key>{-2}};

  SECTION("Check single submap insert/erase")
  {
    thrust::device_vector<Key> d_keys(num_keys);
    thrust::device_vector<Value> d_values(num_keys);
    thrust::device_vector<bool> d_keys_exist(num_keys);

    thrust::sequence(thrust::device, d_keys.begin(), d_keys.end(), 1);
    thrust::sequence(thrust::device, d_values.begin(), d_values.end(), 1);

    auto pairs_begin =
      thrust::make_zip_iterator(cuda::std::tuple{d_keys.begin(), d_values.begin()});

    map.insert(pairs_begin, pairs_begin + num_keys);

    REQUIRE(map.get_size() == num_keys);

    map.erase(d_keys.begin(), d_keys.end());

    // delete decreases count correctly
    REQUIRE(map.get_size() == 0);

    map.contains(d_keys.begin(), d_keys.end(), d_keys_exist.begin());

    // keys were actaully deleted
    REQUIRE(cuco::test::none_of(d_keys_exist.begin(), d_keys_exist.end(), cuda::std::identity{}));

    // ensures that map is reusing deleted slots
    map.insert(pairs_begin, pairs_begin + num_keys);

    REQUIRE(map.get_size() == num_keys);

    map.contains(d_keys.begin(), d_keys.end(), d_keys_exist.begin());

    REQUIRE(cuco::test::all_of(d_keys_exist.begin(), d_keys_exist.end(), cuda::std::identity{}));

    // erase can act selectively
    map.erase(d_keys.begin(), d_keys.begin() + num_keys / 2);
    map.contains(d_keys.begin(), d_keys.end(), d_keys_exist.begin());

    REQUIRE(cuco::test::none_of(
      d_keys_exist.begin(), d_keys_exist.begin() + num_keys / 2, cuda::std::identity{}));

    REQUIRE(cuco::test::all_of(
      d_keys_exist.begin() + num_keys / 2, d_keys_exist.end(), cuda::std::identity{}));

    // clear map
    map.erase(d_keys.begin() + num_keys / 2, d_keys.end());
  }

  SECTION("Check multiple submaps insert/erase")
  {
    constexpr std::size_t num = 4 * num_keys;

    thrust::device_vector<Key> d_keys(num);
    thrust::device_vector<Value> d_values(num);
    thrust::device_vector<bool> d_keys_exist(num);

    thrust::sequence(thrust::device, d_keys.begin(), d_keys.end(), 1);
    thrust::sequence(thrust::device, d_values.begin(), d_values.end(), 1);

    auto pairs_begin =
      thrust::make_zip_iterator(cuda::std::tuple{d_keys.begin(), d_values.begin()});

    map.insert(pairs_begin, pairs_begin + num);

    // map should resize twice if the erased slots are successfully reused
    REQUIRE(map.get_capacity() == 2 * num);
    // check that keys can be successfully deleted from only the first and second submaps
    map.erase(d_keys.begin(), d_keys.begin() + 2 * num_keys);
    map.contains(d_keys.begin(), d_keys.end(), d_keys_exist.begin());

    REQUIRE(cuco::test::none_of(
      d_keys_exist.begin(), d_keys_exist.begin() + 2 * num_keys, cuda::std::identity{}));

    REQUIRE(cuco::test::all_of(
      d_keys_exist.begin() + 2 * num_keys, d_keys_exist.end(), cuda::std::identity{}));

    REQUIRE(map.get_size() == 2 * num_keys);
    // check that keys can be successfully deleted from all submaps (some will be unsuccessful
    // erases)
    map.erase(d_keys.begin(), d_keys.end());

    map.contains(d_keys.begin(), d_keys.end(), d_keys_exist.begin());

    REQUIRE(cuco::test::none_of(d_keys_exist.begin(), d_keys_exist.end(), cuda::std::identity{}));

    REQUIRE(map.get_size() == 0);
  }
}
