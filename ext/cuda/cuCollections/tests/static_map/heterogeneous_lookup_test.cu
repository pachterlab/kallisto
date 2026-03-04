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

#include <cuco/static_map.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/transform.h>

#include <catch2/catch_template_test_macros.hpp>

#include <tuple>

// insert key type
template <typename T>
struct key_pair {
  T a;
  T b;

  __host__ __device__ key_pair() {}
  __host__ __device__ key_pair(T x) : a{x}, b{x} {}

  // Device equality operator is mandatory due to libcudacxx bug:
  // https://github.com/NVIDIA/libcudacxx/issues/223
  __device__ bool operator==(key_pair const& other) const { return a == other.a and b == other.b; }

  __device__ explicit operator T() const noexcept { return a; }
};

// probe key type
template <typename T>
struct key_triplet {
  T a;
  T b;
  T c;

  __host__ __device__ key_triplet() {}
  __host__ __device__ key_triplet(T x) : a{x}, b{x}, c{x} {}

  // Device equality operator is mandatory due to libcudacxx bug:
  // https://github.com/NVIDIA/libcudacxx/issues/223
  __device__ bool operator==(key_triplet const& other) const
  {
    return a == other.a and b == other.b and c == other.c;
  }
};

// User-defined device hasher
struct custom_hasher {
  template <typename CustomKey>
  __device__ uint32_t operator()(CustomKey const& k) const
  {
    return k.a;
  };
};

// User-defined device key equality, Slot key always on the right-hand side
struct custom_key_equal {
  template <typename InputKey, typename SlotKey>
  __device__ bool operator()(InputKey const& lhs, SlotKey const& rhs) const
  {
    return lhs.a == rhs;
  }
};

TEMPLATE_TEST_CASE_SIG("static_map heterogeneous lookup tests",
                       "",
                       ((typename T, int CGSize), T, CGSize),
#if defined(CUCO_HAS_INDEPENDENT_THREADS)  // Key type larger than 8B only supported for sm_70 and
                                           // up
                       (int64_t, 1),
                       (int64_t, 2),
#endif

                       (int32_t, 1),
                       (int32_t, 2))
{
  using Key        = T;
  using Value      = T;
  using InsertKey  = key_pair<T>;
  using ProbeKey   = key_triplet<T>;
  using probe_type = cuco::double_hashing<CGSize, custom_hasher, custom_hasher>;

  auto const sentinel_key   = Key{-1};
  auto const sentinel_value = Value{-1};

  constexpr std::size_t num      = 100;
  constexpr std::size_t capacity = num * 2;
  auto const probe               = probe_type{custom_hasher{}, custom_hasher{}};

  auto my_map = cuco::static_map{capacity,
                                 cuco::empty_key<Key>{sentinel_key},
                                 cuco::empty_value{sentinel_value},
                                 custom_key_equal{},
                                 probe};

  auto insert_pairs = thrust::make_transform_iterator(
    thrust::counting_iterator<int>(0),
    cuda::proclaim_return_type<cuco::pair<InsertKey, Value>>(
      [] __device__(auto i) { return cuco::pair<InsertKey, Value>(i, i); }));
  auto probe_keys = thrust::make_transform_iterator(
    thrust::counting_iterator<int>(0),
    cuda::proclaim_return_type<ProbeKey>([] __device__(auto i) { return ProbeKey{i}; }));

  SECTION("All inserted keys-value pairs should be contained")
  {
    thrust::device_vector<bool> contained(num);
    my_map.insert(insert_pairs, insert_pairs + num);
    my_map.contains(probe_keys, probe_keys + num, contained.begin());
    REQUIRE(cuco::test::all_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  SECTION("Non-inserted keys-value pairs should not be contained")
  {
    thrust::device_vector<bool> contained(num);
    my_map.contains(probe_keys, probe_keys + num, contained.begin());
    REQUIRE(cuco::test::none_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }
}
