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

#include <cuco/static_set.cuh>

#include <cuda/functional>
#include <cuda/std/tuple>
#include <thrust/device_vector.h>
#include <thrust/distance.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sort.h>
#include <thrust/transform.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = int32_t;

int32_t constexpr SENTINEL = -1;

template <typename Set>
void test_unique_sequence(Set& set, size_type num_keys)
{
  using Key = typename Set::key_type;

  auto keys_begin = thrust::counting_iterator<Key>{0};
  thrust::device_vector<bool> d_contained(num_keys);

  auto zip_equal = cuda::proclaim_return_type<bool>(
    [] __device__(auto const& p) { return cuda::std::get<0>(p) == cuda::std::get<1>(p); });
  auto is_even =
    cuda::proclaim_return_type<bool>([] __device__(auto const& i) { return i % 2 == 0; });

  SECTION("Non-inserted keys should not be contained.")
  {
    REQUIRE(set.size() == 0);

    set.contains(keys_begin, keys_begin + num_keys, d_contained.begin());
    REQUIRE(cuco::test::none_of(d_contained.begin(), d_contained.end(), cuda::std::identity{}));
  }

  SECTION("Non-inserted keys have no matches")
  {
    thrust::device_vector<Key> d_results(num_keys);

    set.find(keys_begin, keys_begin + num_keys, d_results.begin());
    auto zip = thrust::make_zip_iterator(cuda::std::tuple{
      d_results.begin(), thrust::constant_iterator<Key>{set.empty_key_sentinel()}});

    REQUIRE(cuco::test::all_of(zip, zip + num_keys, zip_equal));
  }

  SECTION("All conditionally inserted keys should be contained")
  {
    auto const inserted = set.insert_if(
      keys_begin, keys_begin + num_keys, thrust::counting_iterator<std::size_t>(0), is_even);
    REQUIRE(inserted == num_keys / 2);
    REQUIRE(set.size() == num_keys / 2);

    set.contains(keys_begin, keys_begin + num_keys, d_contained.begin());
    REQUIRE(cuco::test::equal(
      d_contained.begin(),
      d_contained.end(),
      thrust::counting_iterator<std::size_t>(0),
      cuda::proclaim_return_type<bool>([] __device__(auto const& idx_contained, auto const& idx) {
        return ((idx % 2) == 0) == idx_contained;
      })));
  }

  set.clear();
  auto const inserted = set.insert(keys_begin, keys_begin + num_keys);
  REQUIRE(inserted == num_keys);
  REQUIRE(set.size() == num_keys);

  SECTION("All inserted keys should be contained.")
  {
    set.contains(keys_begin, keys_begin + num_keys, d_contained.begin());
    REQUIRE(cuco::test::all_of(d_contained.begin(), d_contained.end(), cuda::std::identity{}));
  }

  SECTION("Conditional contains should return true on even inputs.")
  {
    set.contains_if(keys_begin,
                    keys_begin + num_keys,
                    thrust::counting_iterator<std::size_t>(0),
                    is_even,
                    d_contained.begin());
    auto gold_iter =
      thrust::make_transform_iterator(thrust::counting_iterator<std::size_t>(0), is_even);
    auto zip = thrust::make_zip_iterator(cuda::std::tuple{d_contained.begin(), gold_iter});
    REQUIRE(cuco::test::all_of(zip, zip + num_keys, zip_equal));
  }

  SECTION("All inserted keys should be correctly recovered during find")
  {
    thrust::device_vector<Key> d_results(num_keys);

    set.find(keys_begin, keys_begin + num_keys, d_results.begin());
    auto zip = thrust::make_zip_iterator(cuda::std::tuple{d_results.begin(), keys_begin});

    REQUIRE(cuco::test::all_of(zip, zip + num_keys, zip_equal));
  }

  SECTION("Conditional find should return valid values on even inputs.")
  {
    auto found_results = thrust::device_vector<Key>(num_keys);
    auto gold_fn       = cuda::proclaim_return_type<Key>(
      [] __device__(auto const& i) { return i % 2 == 0 ? static_cast<Key>(i) : Key{SENTINEL}; });

    set.find_if(keys_begin,
                keys_begin + num_keys,
                thrust::counting_iterator<std::size_t>{0},
                is_even,
                found_results.begin());

    REQUIRE(cuco::test::equal(
      found_results.begin(),
      found_results.end(),
      thrust::make_transform_iterator(thrust::counting_iterator<Key>{0}, gold_fn),
      cuda::proclaim_return_type<bool>(
        [] __device__(auto const& found, auto const& gold) { return found == gold; })));
  }
}

TEMPLATE_TEST_CASE_SIG(
  "static_set unique sequence tests",
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
  constexpr size_type num_keys{400};
  using probe = std::conditional_t<Probe == cuco::test::probe_sequence::linear_probing,
                                   cuco::linear_probing<CGSize, cuco::default_hash_function<Key>>,
                                   cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>>;

  constexpr size_type gold_capacity = [&]() {
    if constexpr (cuco::is_double_hashing<probe>::value) {
      return (CGSize == 1) ? 422   // 211 x 1 x 2
                           : 412;  // 103 x 2 x 2
    } else {
      return 400;
    }
  }();

  auto set =
    cuco::static_set{num_keys, cuco::empty_key<Key>{SENTINEL}, {}, probe{}, {}, cuco::storage<2>{}};

  REQUIRE(set.capacity() == gold_capacity);

  test_unique_sequence(set, num_keys);
}
