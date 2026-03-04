/*
 * Copyright (c) 2020-2025, NVIDIA CORPORATION.
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
#include <cuda/std/functional>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/transform.h>

#include <catch2/catch_template_test_macros.hpp>

#include <cstddef>
#include <cstdint>
#include <tuple>

// User-defined key type
template <typename T>
struct key_pair_type {
  T a;
  T b;

  __host__ __device__ key_pair_type() {}
  __host__ __device__ key_pair_type(T x) : a{x}, b{x} {}

  // Device equality operator is mandatory due to libcudacxx bug:
  // https://github.com/NVIDIA/libcudacxx/issues/223
  __device__ bool operator==(key_pair_type const& other) const
  {
    return a == other.a and b == other.b;
  }
};

// User-defined key type
template <typename T>
struct large_key_type {
  T a;
  T b;
  T c;

  __host__ __device__ large_key_type() {}
  __host__ __device__ large_key_type(T x) : a{x}, b{x}, c{x} {}

  // Device equality operator is mandatory due to libcudacxx bug:
  // https://github.com/NVIDIA/libcudacxx/issues/223
  __device__ bool operator==(large_key_type const& other) const
  {
    return a == other.a and b == other.b and c == other.c;
  }
};

// User-defined value type
template <typename T>
struct value_pair_type {
  T f;
  T s;

  __host__ __device__ value_pair_type() {}
  __host__ __device__ value_pair_type(T x) : f{x}, s{x} {}

  __device__ bool operator==(value_pair_type const& other) const
  {
    return f == other.f and s == other.s;
  }
};

// User-defined device hasher
struct hash_custom_key {
  template <typename custom_type>
  __device__ uint32_t operator()(custom_type k)
  {
    return thrust::raw_reference_cast(k).a;
  };
};

// User-defined device key equality
struct custom_key_equals {
  template <typename lhs_type, typename rhs_type>
  __device__ bool operator()(lhs_type lhs, rhs_type rhs)
  {
    return lhs == static_cast<lhs_type>(rhs);
  }
};

TEMPLATE_TEST_CASE_SIG("static_map custom key and value type tests",
                       "",
                       ((typename Key, typename Value), Key, Value),
#if defined(CUCO_HAS_INDEPENDENT_THREADS)  // Key type larger than 8B only supported for sm_70 and
                                           // up
                       (key_pair_type<int64_t>, value_pair_type<int32_t>),
                       (key_pair_type<int64_t>, value_pair_type<int64_t>),
                       (large_key_type<int32_t>, value_pair_type<int32_t>),
#endif
                       (key_pair_type<int32_t>, value_pair_type<int32_t>))
{
  auto const sentinel_key   = Key{-1};
  auto const sentinel_value = Value{-1};

  constexpr std::size_t num      = 100;
  constexpr std::size_t capacity = num * 2;
  cuco::legacy::static_map<Key, Value> map{
    capacity, cuco::empty_key<Key>{sentinel_key}, cuco::empty_value<Value>{sentinel_value}};

  auto insert_keys = thrust::make_transform_iterator(
    thrust::counting_iterator<int>(0),
    cuda::proclaim_return_type<Key>([] __device__(auto i) { return Key{i}; }));

  auto insert_values = thrust::make_transform_iterator(
    thrust::counting_iterator<int>(0),
    cuda::proclaim_return_type<Value>([] __device__(auto i) { return Value{i}; }));

  auto insert_pairs = thrust::make_transform_iterator(
    thrust::make_counting_iterator<int>(0),
    cuda::proclaim_return_type<cuco::pair<Key, Value>>(
      [] __device__(auto i) { return cuco::pair<Key, Value>(i, i); }));

  SECTION("All inserted keys-value pairs should be correctly recovered during find")
  {
    thrust::device_vector<Value> found_values(num);
    map.insert(insert_pairs, insert_pairs + num, hash_custom_key{}, custom_key_equals{});

    REQUIRE(num == map.get_size());

    map.find(
      insert_keys, insert_keys + num, found_values.begin(), hash_custom_key{}, custom_key_equals{});

    REQUIRE(cuco::test::equal(insert_values,
                              insert_values + num,
                              found_values.begin(),
                              cuda::proclaim_return_type<bool>([] __device__(Value lhs, Value rhs) {
                                return cuda::std::tie(lhs.f, lhs.s) == cuda::std::tie(rhs.f, rhs.s);
                              })));
  }

  SECTION("All inserted keys-value pairs should be contained")
  {
    thrust::device_vector<bool> contained(num);
    map.insert(insert_pairs, insert_pairs + num, hash_custom_key{}, custom_key_equals{});
    map.contains(
      insert_keys, insert_keys + num, contained.begin(), hash_custom_key{}, custom_key_equals{});
    REQUIRE(cuco::test::all_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  SECTION("All conditionally inserted keys-value pairs should be contained")
  {
    thrust::device_vector<bool> contained(num);
    map.insert_if(
      insert_pairs,
      insert_pairs + num,
      thrust::counting_iterator<int>(0),
      cuda::proclaim_return_type<bool>([] __device__(auto const& key) { return (key % 2) == 0; }),
      hash_custom_key{},
      custom_key_equals{});

    REQUIRE(num / 2 == map.get_size());

    map.contains(
      insert_keys, insert_keys + num, contained.begin(), hash_custom_key{}, custom_key_equals{});

    REQUIRE(cuco::test::equal(
      contained.begin(),
      contained.end(),
      thrust::counting_iterator<int>(0),
      cuda::proclaim_return_type<bool>([] __device__(auto const& idx_contained, auto const& idx) {
        return ((idx % 2) == 0) == idx_contained;
      })));
  }

  SECTION("Non-inserted keys-value pairs should not be contained")
  {
    thrust::device_vector<bool> contained(num);
    map.contains(
      insert_keys, insert_keys + num, contained.begin(), hash_custom_key{}, custom_key_equals{});
    REQUIRE(cuco::test::none_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  SECTION("All inserted keys-value pairs should be contained")
  {
    thrust::device_vector<bool> contained(num);
    map.insert(insert_pairs, insert_pairs + num, hash_custom_key{}, custom_key_equals{});
    auto view = map.get_device_view();
    REQUIRE(cuco::test::all_of(
      insert_pairs,
      insert_pairs + num,
      cuda::proclaim_return_type<bool>([view] __device__(cuco::pair<Key, Value> const& pair) {
        return view.contains(pair.first, hash_custom_key{}, custom_key_equals{});
      })));
  }

  SECTION("Inserting unique keys should return insert success.")
  {
    auto m_view = map.get_device_mutable_view();
    REQUIRE(cuco::test::all_of(insert_pairs,
                               insert_pairs + num,
                               cuda::proclaim_return_type<bool>(
                                 [m_view] __device__(cuco::pair<Key, Value> const& pair) mutable {
                                   return m_view.insert(
                                     pair, hash_custom_key{}, custom_key_equals{});
                                 })));
  }

  SECTION("Cannot find any key in an empty hash map")
  {
    SECTION("non-const view")
    {
      auto view = map.get_device_view();
      REQUIRE(cuco::test::all_of(
        insert_pairs,
        insert_pairs + num,
        cuda::proclaim_return_type<bool>(
          [view] __device__(cuco::pair<Key, Value> const& pair) mutable {
            return view.find(pair.first, hash_custom_key{}, custom_key_equals{}) == view.end();
          })));
    }

    SECTION("const view")
    {
      auto const view = map.get_device_view();
      REQUIRE(cuco::test::all_of(
        insert_pairs,
        insert_pairs + num,
        cuda::proclaim_return_type<bool>([view] __device__(cuco::pair<Key, Value> const& pair) {
          return view.find(pair.first, hash_custom_key{}, custom_key_equals{}) == view.end();
        })));
    }
  }
}
