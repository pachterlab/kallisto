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

#include <cuco/static_map.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/logical.h>
#include <thrust/transform.h>

// User-defined key type
struct custom_key_type {
  int32_t a;
  int32_t b;

  __host__ __device__ custom_key_type() {}
  __host__ __device__ custom_key_type(int32_t x) : a{x}, b{x} {}
};

// User-defined value type
struct custom_value_type {
  int32_t f;
  int32_t s;

  __host__ __device__ custom_value_type() {}
  __host__ __device__ custom_value_type(int32_t x) : f{x}, s{x} {}
};

// User-defined device hash callable
struct custom_hash {
  __device__ uint32_t operator()(custom_key_type const& k) const noexcept { return k.a; };
};

// User-defined device key equal callable
struct custom_key_equal {
  __device__ bool operator()(custom_key_type const& lhs, custom_key_type const& rhs) const noexcept
  {
    return lhs.a == rhs.a;
  }
};

int main(void)
{
  constexpr std::size_t num_pairs = 80'000;

  // Set emtpy sentinels
  auto const empty_key_sentinel   = custom_key_type{-1};
  auto const empty_value_sentinel = custom_value_type{-1};

  // Create an iterator of input key/value pairs
  auto pairs_begin = thrust::make_transform_iterator(
    thrust::make_counting_iterator<int32_t>(0),
    cuda::proclaim_return_type<cuco::pair<custom_key_type, custom_value_type>>(
      [] __device__(auto i) { return cuco::pair{custom_key_type{i}, custom_value_type{i}}; }));

  // Construct a map with 100,000 slots using the given empty key/value sentinels. Note the
  // capacity is chosen knowing we will insert 80,000 keys, for an load factor of 80%.
  auto map = cuco::static_map{cuco::extent<std::size_t, 100'000>{},
                              cuco::empty_key{empty_key_sentinel},
                              cuco::empty_value{empty_value_sentinel},
                              custom_key_equal{},
                              cuco::linear_probing<1, custom_hash>{}};

  // Inserts 80,000 pairs into the map by using the custom hasher and custom equality callable
  map.insert(pairs_begin, pairs_begin + num_pairs);

  // Reproduce inserted keys
  auto insert_keys =
    thrust::make_transform_iterator(thrust::make_counting_iterator<int32_t>(0),
                                    cuda::proclaim_return_type<custom_key_type>(
                                      [] __device__(auto i) { return custom_key_type{i}; }));

  thrust::device_vector<bool> contained(num_pairs);

  // Determine if all the inserted keys can be found by using the same hasher and equality
  // function as `insert`. If a key `insert_keys[i]` doesn't exist, `contained[i] == false`.
  map.contains(insert_keys, insert_keys + num_pairs, contained.begin());
  // This will fail due to inconsistent hash and key equal.
  // map.contains(insert_keys, insert_keys + num_pairs, contained.begin());

  // All inserted keys are contained
  auto const all_contained =
    thrust::all_of(contained.begin(), contained.end(), [] __device__(auto const& b) { return b; });
  if (all_contained) { std::cout << "Success! Found all values.\n"; }

  return 0;
}
