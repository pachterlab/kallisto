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

#include <cuco/static_multimap.cuh>

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
    return thrust::raw_reference_cast(k).a;
  };
};

// User-defined device key equality
struct custom_key_equal {
  template <typename LHS, typename RHS>
  __device__ bool operator()(LHS const& lhs, RHS const& rhs)
  {
    return thrust::raw_reference_cast(lhs).a == thrust::raw_reference_cast(rhs).a;
  }
};

TEMPLATE_TEST_CASE("static_multimap heterogeneous lookup tests",
                   "",
#if defined(CUCO_HAS_INDEPENDENT_THREADS)  // Key type larger than 8B only supported for sm_70 and
                                           // up
                   int64_t,
#endif
                   int32_t)
{
  using Key      = key_pair<TestType>;
  using Value    = TestType;
  using ProbeKey = key_triplet<TestType>;

  auto const sentinel_key   = Key{-1};
  auto const sentinel_value = Value{-1};

  constexpr std::size_t num      = 100;
  constexpr std::size_t capacity = num * 2;
  cuco::static_multimap<Key,
                        Value,
                        cuda::thread_scope_device,
                        cuco::cuda_allocator<char>,
                        cuco::legacy::linear_probing<1, custom_hasher>>
    map{capacity, cuco::empty_key<Key>{sentinel_key}, cuco::empty_value<Value>{sentinel_value}};

  auto insert_pairs = thrust::make_transform_iterator(
    thrust::counting_iterator<int>(0),
    cuda::proclaim_return_type<cuco::pair<Key, Value>>(
      [] __device__(auto i) { return cuco::pair<Key, Value>(i, i); }));
  auto probe_keys = thrust::make_transform_iterator(
    thrust::counting_iterator<int>(0),
    cuda::proclaim_return_type<ProbeKey>([] __device__(auto i) { return ProbeKey(i); }));

  SECTION("All inserted keys-value pairs should be contained")
  {
    thrust::device_vector<bool> contained(num);
    map.insert(insert_pairs, insert_pairs + num);
    map.contains(probe_keys, probe_keys + num, contained.begin(), custom_key_equal{});
    REQUIRE(cuco::test::all_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }

  SECTION("Non-inserted keys-value pairs should not be contained")
  {
    thrust::device_vector<bool> contained(num);
    map.contains(probe_keys, probe_keys + num, contained.begin(), custom_key_equal{});
    REQUIRE(cuco::test::none_of(contained.begin(), contained.end(), cuda::std::identity{}));
  }
}
