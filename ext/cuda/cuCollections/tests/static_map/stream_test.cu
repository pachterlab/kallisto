/*
 * Copyright (c) 2021-2025, NVIDIA CORPORATION.
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
#include <cuda/std/tuple>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>

TEMPLATE_TEST_CASE_SIG("static_map: unique sequence of keys on given stream",
                       "",
                       ((typename Key, typename Value), Key, Value),
                       (int32_t, int32_t),
                       (int32_t, int64_t),
                       (int64_t, int32_t),
                       (int64_t, int64_t))
{
  cudaStream_t stream;
  CUCO_CUDA_TRY(cudaStreamCreate(&stream));

  {  // Scope ensures map is destroyed before stream
    constexpr std::size_t num_keys{500'000};
    auto map = cuco::static_map{num_keys * 2,
                                cuco::empty_key<Key>{-1},
                                cuco::empty_value<Value>{-1},
                                {},
                                cuco::linear_probing<1, cuco::default_hash_function<Key>>{},
                                {},
                                {},
                                {},
                                stream};

    thrust::device_vector<Key> d_keys(num_keys);
    thrust::device_vector<Value> d_values(num_keys);

    thrust::sequence(thrust::device, d_keys.begin(), d_keys.end());
    thrust::sequence(thrust::device, d_values.begin(), d_values.end());

    auto pairs_begin = thrust::make_transform_iterator(
      thrust::make_counting_iterator<int>(0),
      cuda::proclaim_return_type<cuco::pair<Key, Value>>(
        [] __device__(auto i) { return cuco::pair<Key, Value>(i, i); }));

    // bulk function test cases
    SECTION("All inserted keys-value pairs should be correctly recovered during find")
    {
      thrust::device_vector<Value> d_results(num_keys);

      map.insert(pairs_begin, pairs_begin + num_keys, stream);
      map.find(d_keys.begin(), d_keys.end(), d_results.begin(), stream);
      auto zip = thrust::make_zip_iterator(cuda::std::tuple{d_results.begin(), d_values.begin()});

      REQUIRE(cuco::test::all_of(zip,
                                 zip + num_keys,
                                 cuda::proclaim_return_type<bool>([] __device__(auto const& p) {
                                   return cuda::std::get<0>(p) == cuda::std::get<1>(p);
                                 }),
                                 stream));
    }

    SECTION("All inserted keys-value pairs should be contained")
    {
      thrust::device_vector<bool> d_contained(num_keys);

      map.insert(pairs_begin, pairs_begin + num_keys, stream);
      map.contains(d_keys.begin(), d_keys.end(), d_contained.begin(), stream);

      REQUIRE(
        cuco::test::all_of(d_contained.begin(), d_contained.end(), cuda::std::identity{}, stream));
    }
  }  // map is destroyed here

  CUCO_CUDA_TRY(cudaStreamDestroy(stream));
}
