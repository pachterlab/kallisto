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

#include <cuco/static_map.cuh>

#include <cuda/functional>
#include <cuda/std/tuple>
#include <thrust/detail/raw_reference_cast.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <iostream>

/**
 * @file heterogeneous_lookup_example.cu
 *
 * @brief Demonstrates how to perform heterogeneous lookups with `cuco::static_map`.
 *
 * This example demonstrates heterogeneous lookup, which allows you to perform lookups with a key
 * type that is different from the container's key type, without having to first construct an
 * object of the container's key type.
 *
 * In many workflows the format of the keys used when inserting into a hash table differs from the
 * format that is available at query time. This example stores keys as `cuco::pair<int, int>`
 * representing `(sensor_id, channel)` but performs lookups directly using 3-element tuples
 * `(sensor_id, channel, timestamp)` without needing to construct intermediate `cuco::pair` objects.
 *
 * Heterogeneous lookup is enabled by custom hash and equality functors that can operate on
 * "compatible" key types. The functors only consider the first two elements, allowing both
 * the stored `cuco::pair` and query `tuple` types to interoperate transparently and efficiently.
 */

using stored_key = cuco::pair<int, int>;  // Key type used for insertion: (sensor_id, channel)
using probe_key =
  cuda::std::tuple<int, int, int>;  // Key type used for querying: (sensor_id, channel, timestamp)
using value_type = float;

// Declare that value_type is bitwise comparable since float doesn't have unique object
// representations
CUCO_DECLARE_BITWISE_COMPARABLE(value_type);

// Heterogeneous hasher that can hash both cuco::pair and tuple types without conversion.
// The template allows it to accept any key type and extract the first two elements.
struct heterogeneous_hasher {
  template <typename Key>
  __device__ std::size_t operator()(Key const& key) const
  {
    auto const& ref  = thrust::raw_reference_cast(key);
    auto const major = cuda::std::get<0>(ref);  // Works for both pair.first and get<0>(tuple)
    auto const minor = cuda::std::get<1>(ref);  // Works for both pair.second and get<1>(tuple)
    return static_cast<std::size_t>(major * 131 + minor);
  }
};

// Heterogeneous equality functor that can compare cuco::pair and tuple types without conversion.
// The template allows it to accept any combination of key types and compare their first two
// elements.
struct heterogeneous_key_equal {
  template <typename LHS, typename RHS>
  __device__ bool operator()(LHS const& lhs, RHS const& rhs) const
  {
    auto const& left  = thrust::raw_reference_cast(lhs);
    auto const& right = thrust::raw_reference_cast(rhs);
    return (cuda::std::get<0>(left) == cuda::std::get<0>(right)) and  // Compare first elements
           (cuda::std::get<1>(left) == cuda::std::get<1>(right));     // Compare second elements
  }
};

int main()
{
  constexpr std::size_t num_entries = 4;
  auto constexpr empty_key          = stored_key{-1, -1};
  auto constexpr empty_value        = value_type{-1.0f};

  // Allocate a map with ~50% load factor.
  auto map =
    cuco::static_map{cuco::extent<std::size_t>{num_entries * 2},
                     cuco::empty_key{empty_key},
                     cuco::empty_value{empty_value},
                     heterogeneous_key_equal{},
                     cuco::linear_probing<1, heterogeneous_hasher>{heterogeneous_hasher{}}};

  thrust::device_vector<stored_key> d_keys = {
    stored_key{101, 3},
    stored_key{104, 8},
    stored_key{215, 1},
    stored_key{305, 0},
  };
  thrust::device_vector<value_type> d_values = {36.5f, 41.2f, 27.1f, 33.8f};

  auto pairs_begin = thrust::make_transform_iterator(
    thrust::counting_iterator{0},
    cuda::proclaim_return_type<cuco::pair<stored_key, value_type>>(
      [keys = d_keys.begin(), values = d_values.begin()] __device__(int i) {
        return cuco::pair<stored_key, value_type>{keys[i], values[i]};
      }));

  map.insert(pairs_begin, pairs_begin + num_entries);

  // Query using 3-element tuples that include an additional timestamp field.
  // The heterogeneous hash and equality functors only consider the first two components
  // (sensor_id, channel) when comparing against the stored cuco::pair keys.
  thrust::device_vector<probe_key> d_queries{
    probe_key{101, 3, 1210},  // present in the map
    probe_key{215, 1, 1345},  // present in the map
    probe_key{999, 4, 2000},  // missing entry
  };

  thrust::device_vector<bool> d_contains(d_queries.size());
  map.contains(d_queries.begin(), d_queries.end(), d_contains.begin());

  thrust::device_vector<value_type> d_found(d_queries.size());
  map.find(d_queries.begin(), d_queries.end(), d_found.begin());

  // Copy results back to host for printing
  thrust::host_vector<probe_key> h_queries = d_queries;
  thrust::host_vector<bool> h_contains     = d_contains;
  thrust::host_vector<value_type> h_found  = d_found;

  for (std::size_t i = 0; i < h_queries.size(); ++i) {
    auto const& query  = h_queries[i];
    auto const present = h_contains[i];
    std::cout << "Lookup tuple (sensor " << cuda::std::get<0>(query) << ", channel "
              << cuda::std::get<1>(query) << ", timestamp " << cuda::std::get<2>(query) << ") -> "
              << (present ? "found" : "missing");

    if (present) { std::cout << ", stored value = " << h_found[i]; }

    std::cout << "\n";
  }

  return 0;
}
