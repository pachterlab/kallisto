/*
 * Copyright (c) 2022-2024, NVIDIA CORPORATION.
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

#include <cuco/static_set.cuh>

#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("static_set size test", "")
{
  constexpr std::size_t num_keys{400};

  cuco::static_set<int> set{cuco::extent<std::size_t>{400}, cuco::empty_key{-1}};

  thrust::device_vector<int> d_keys(num_keys);

  thrust::sequence(thrust::device, d_keys.begin(), d_keys.end());

  auto const num_successes = set.insert(d_keys.begin(), d_keys.end());

  REQUIRE(set.size() == num_keys);
  REQUIRE(num_successes == num_keys);

  set.clear();

  REQUIRE(set.size() == 0);
}
