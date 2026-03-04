/*
 * Copyright (c) 2025, NVIDIA CORPORATION.
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

#include <cuda/functional>
#include <cuda/std/functional>
#include <cuda/std/tuple>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>

TEMPLATE_TEST_CASE_SIG("dynamic_map find tests",
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

  SECTION("Check single submap insert/find")
  {
    thrust::device_vector<Key> d_keys(num_keys);
    thrust::device_vector<Value> d_values(num_keys);
    thrust::device_vector<Value> d_found_values(num_keys);

    thrust::sequence(thrust::device, d_keys.begin(), d_keys.end(), 1);
    thrust::sequence(thrust::device, d_values.begin(), d_values.end(), 1);

    auto pairs_begin =
      thrust::make_zip_iterator(cuda::std::tuple{d_keys.begin(), d_values.begin()});

    map.insert(pairs_begin, pairs_begin + num_keys);

    REQUIRE(map.get_size() == num_keys);

    // Find all inserted keys
    map.find(d_keys.begin(), d_keys.end(), d_found_values.begin());

    // Verify that all keys were found with correct values
    auto zip_equal = cuda::proclaim_return_type<bool>(
      [] __device__(auto const& p) { return cuda::std::get<0>(p) == cuda::std::get<1>(p); });
    auto zip =
      thrust::make_zip_iterator(cuda::std::tuple{d_values.begin(), d_found_values.begin()});
    REQUIRE(cuco::test::all_of(zip, zip + num_keys, zip_equal));

    // Test finding non-existent keys
    thrust::device_vector<Key> d_nonexistent_keys(100);
    thrust::device_vector<Value> d_nonexistent_values(100);

    thrust::sequence(thrust::device,
                     d_nonexistent_keys.begin(),
                     d_nonexistent_keys.end(),
                     static_cast<Key>(num_keys + 1));

    map.find(d_nonexistent_keys.begin(), d_nonexistent_keys.end(), d_nonexistent_values.begin());

    // Verify that non-existent keys return empty value sentinel
    auto empty_zip = thrust::make_zip_iterator(
      cuda::std::tuple{d_nonexistent_values.begin(),
                       thrust::constant_iterator<Value>{cuco::empty_value<Value>{-1}.value}});
    REQUIRE(cuco::test::all_of(empty_zip, empty_zip + 100, zip_equal));

    // Test finding a mix of existing and non-existing keys
    thrust::device_vector<Key> d_mixed_keys(200);
    thrust::device_vector<Value> d_mixed_values(200);

    // First half: existing keys
    thrust::copy(d_keys.begin(), d_keys.begin() + 100, d_mixed_keys.begin());
    // Second half: non-existing keys
    thrust::sequence(thrust::device,
                     d_mixed_keys.begin() + 100,
                     d_mixed_keys.end(),
                     static_cast<Key>(num_keys + 1));

    map.find(d_mixed_keys.begin(), d_mixed_keys.end(), d_mixed_values.begin());

    // Verify first half found correct values
    auto first_half_zip =
      thrust::make_zip_iterator(cuda::std::tuple{d_values.begin(), d_mixed_values.begin()});
    REQUIRE(cuco::test::all_of(first_half_zip, first_half_zip + 100, zip_equal));

    // Verify second half returned empty value sentinel
    auto second_half_empty_zip = thrust::make_zip_iterator(
      cuda::std::tuple{d_mixed_values.begin() + 100,
                       thrust::constant_iterator<Value>{cuco::empty_value<Value>{-1}.value}});
    REQUIRE(cuco::test::all_of(second_half_empty_zip, second_half_empty_zip + 100, zip_equal));
  }

  SECTION("Check find after erase")
  {
    thrust::device_vector<Key> d_keys(num_keys);
    thrust::device_vector<Value> d_values(num_keys);
    thrust::device_vector<Value> d_found_values(num_keys);

    thrust::sequence(thrust::device, d_keys.begin(), d_keys.end(), 1);
    thrust::sequence(thrust::device, d_values.begin(), d_values.end(), 1);

    auto pairs_begin =
      thrust::make_zip_iterator(cuda::std::tuple{d_keys.begin(), d_values.begin()});

    map.insert(pairs_begin, pairs_begin + num_keys);

    REQUIRE(map.get_size() == num_keys);

    // Find all keys before erase
    map.find(d_keys.begin(), d_keys.end(), d_found_values.begin());

    auto zip_equal = cuda::proclaim_return_type<bool>(
      [] __device__(auto const& p) { return cuda::std::get<0>(p) == cuda::std::get<1>(p); });
    auto zip =
      thrust::make_zip_iterator(cuda::std::tuple{d_values.begin(), d_found_values.begin()});
    REQUIRE(cuco::test::all_of(zip, zip + num_keys, zip_equal));

    // Erase first half of keys
    map.erase(d_keys.begin(), d_keys.begin() + num_keys / 2);

    REQUIRE(map.get_size() == num_keys / 2);

    // Find all keys after erase
    map.find(d_keys.begin(), d_keys.end(), d_found_values.begin());

    // First half should return empty value sentinel (erased)
    auto first_half_empty_zip = thrust::make_zip_iterator(
      cuda::std::tuple{d_found_values.begin(),
                       thrust::constant_iterator<Value>{cuco::empty_value<Value>{-1}.value}});
    REQUIRE(
      cuco::test::all_of(first_half_empty_zip, first_half_empty_zip + num_keys / 2, zip_equal));

    // Second half should return correct values (not erased)
    auto second_half_zip = thrust::make_zip_iterator(
      cuda::std::tuple{d_values.begin() + num_keys / 2, d_found_values.begin() + num_keys / 2});
    REQUIRE(cuco::test::all_of(second_half_zip, second_half_zip + num_keys / 2, zip_equal));
  }

  cuco::dynamic_map<Key, Value> indentity_hash_map{
    num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{-1}, cuco::erased_key<Key>{-2}};

  SECTION("Check find in a all erased submap")
  {
    constexpr float default_load_factor     = 0.60;
    constexpr std::size_t first_insert_size = num_keys * default_load_factor;

    thrust::device_vector<Key> d_keys(num_keys);
    thrust::device_vector<Value> d_values(num_keys);
    thrust::device_vector<Value> d_found_values(num_keys);

    thrust::sequence(thrust::device, d_keys.begin(), d_keys.end(), 1);
    thrust::sequence(thrust::device, d_values.begin(), d_values.end(), 1);

    auto pairs_begin =
      thrust::make_zip_iterator(cuda::std::tuple{d_keys.begin(), d_values.begin()});

    // To construct a map with all erased keys, we can't insert all at once
    indentity_hash_map.insert(
      pairs_begin, pairs_begin + first_insert_size, cuco::identity_hash<Key>());
    REQUIRE(indentity_hash_map.get_size() == first_insert_size);

    indentity_hash_map.erase(
      d_keys.begin(), d_keys.begin() + first_insert_size, cuco::identity_hash<Key>());
    REQUIRE(indentity_hash_map.get_size() == 0);

    indentity_hash_map.insert(
      pairs_begin + first_insert_size, pairs_begin + num_keys, cuco::identity_hash<Key>());
    REQUIRE(indentity_hash_map.get_size() == num_keys - first_insert_size);

    indentity_hash_map.erase(
      d_keys.begin() + first_insert_size, d_keys.end(), cuco::identity_hash<Key>());
    REQUIRE(indentity_hash_map.get_size() == 0);

    // we've construct a dynamic_map with one submap whose keys are all erased keys (-2 in this
    // case) this find would run forever if we don't check whether we have iterated all keys in a
    // submap
    indentity_hash_map.find(
      d_keys.begin(), d_keys.end(), d_found_values.begin(), cuco::identity_hash<Key>());
    // all d_found_values should be empty value sentinel (-1 in this case)
    auto empty_zip = thrust::make_zip_iterator(
      cuda::std::tuple{d_found_values.begin(),
                       thrust::constant_iterator<Value>{cuco::empty_value<Value>{-1}.value}});
    auto zip_equal = cuda::proclaim_return_type<bool>(
      [] __device__(auto const& p) { return cuda::std::get<0>(p) == cuda::std::get<1>(p); });
    REQUIRE(cuco::test::all_of(empty_zip, empty_zip + num_keys, zip_equal));
  }
}