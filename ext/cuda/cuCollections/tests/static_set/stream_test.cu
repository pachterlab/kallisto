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

#include <cuco/static_set.cuh>

#include <cuda/functional>
#include <cuda/std/iterator>
#include <cuda/std/tuple>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>

TEMPLATE_TEST_CASE_SIG("static_set: operations on different stream than constructor",
                       "",
                       ((typename Key), Key),
                       (int32_t),
                       (int64_t))
{
  cudaStream_t constructor_stream;
  cudaStream_t operation_stream;
  CUCO_CUDA_TRY(cudaStreamCreate(&constructor_stream));
  CUCO_CUDA_TRY(cudaStreamCreate(&operation_stream));

  {  // Scope ensures set is destroyed before streams
    constexpr std::size_t num_keys{500'000};
    auto set = cuco::static_set{num_keys * 2,
                                cuco::empty_key<Key>{-1},
                                {},
                                cuco::linear_probing<1, cuco::default_hash_function<Key>>{},
                                {},
                                {},
                                {},
                                constructor_stream};

    thrust::device_vector<Key> d_keys(num_keys);
    thrust::sequence(thrust::device, d_keys.begin(), d_keys.end());

    SECTION("Insert and contains on different stream than constructor")
    {
      set.insert(d_keys.begin(), d_keys.end(), operation_stream);

      thrust::device_vector<bool> d_contained(num_keys);
      set.contains(d_keys.begin(), d_keys.end(), d_contained.begin(), operation_stream);

      REQUIRE(cuco::test::all_of(
        d_contained.begin(), d_contained.end(), cuda::std::identity{}, operation_stream));
    }

    SECTION("Insert and find on different stream than constructor")
    {
      set.insert(d_keys.begin(), d_keys.end(), operation_stream);

      thrust::device_vector<Key> d_results(num_keys);
      set.find(d_keys.begin(), d_keys.end(), d_results.begin(), operation_stream);

      auto zip =
        thrust::make_zip_iterator(cuda::std::make_tuple(d_results.begin(), d_keys.begin()));
      REQUIRE(cuco::test::all_of(zip,
                                 zip + num_keys,
                                 cuda::proclaim_return_type<bool>([] __device__(auto const& p) {
                                   return cuda::std::get<0>(p) == cuda::std::get<1>(p);
                                 }),
                                 operation_stream));
    }

    SECTION("Insert and retrieve on different stream than constructor")
    {
      set.insert(d_keys.begin(), d_keys.end(), operation_stream);

      thrust::device_vector<Key> d_probe_results(num_keys);
      thrust::device_vector<Key> d_match_results(num_keys);
      auto [d_probe_end, d_match_end] = set.retrieve(d_keys.begin(),
                                                     d_keys.end(),
                                                     d_probe_results.begin(),
                                                     d_match_results.begin(),
                                                     operation_stream);

      auto const num_retrieved = cuda::std::distance(d_probe_results.begin(), d_probe_end);
      REQUIRE(num_retrieved == num_keys);
    }

    SECTION("Insert and size on different stream than constructor")
    {
      set.insert(d_keys.begin(), d_keys.end(), operation_stream);

      auto const size = set.size(operation_stream);
      REQUIRE(size == num_keys);
    }

    SECTION("Insert on constructor stream and query on different stream")
    {
      set.insert(d_keys.begin(), d_keys.end(), constructor_stream);

      thrust::device_vector<bool> d_contained(num_keys);
      set.contains(d_keys.begin(), d_keys.end(), d_contained.begin(), operation_stream);

      REQUIRE(cuco::test::all_of(
        d_contained.begin(), d_contained.end(), cuda::std::identity{}, operation_stream));
    }
  }  // set is destroyed here

  CUCO_CUDA_TRY(cudaStreamDestroy(operation_stream));
  CUCO_CUDA_TRY(cudaStreamDestroy(constructor_stream));
}
