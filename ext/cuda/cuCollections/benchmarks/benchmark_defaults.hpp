/*
 * Copyright (c) 2023-2024, NVIDIA CORPORATION.
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

#pragma once

#include <cuco/hash_functions.cuh>

#include <nvbench/nvbench.cuh>

#include <cstdint>
#include <vector>

namespace cuco::benchmark::defaults {

using KEY_TYPE_RANGE   = nvbench::type_list<nvbench::int32_t, nvbench::int64_t>;
using VALUE_TYPE_RANGE = nvbench::type_list<nvbench::int32_t, nvbench::int64_t>;
using HASH_RANGE       = nvbench::type_list<cuco::identity_hash<char>,
                                            cuco::xxhash_32<char>,
                                            cuco::xxhash_64<char>,
                                            cuco::murmurhash3_32<char>>;  //,
// cuco::murmurhash3_x86_128<char>,
// cuco::murmurhash3_x64_128<char>>; // TODO handle tuple-like hash value

auto constexpr N             = 100'000'000;
auto constexpr OCCUPANCY     = 0.5;
auto constexpr MULTIPLICITY  = 1;
auto constexpr MATCHING_RATE = 1.0;
auto constexpr SKEW          = 0.5;
auto constexpr BATCH_SIZE    = 1'000'000;
auto constexpr INITIAL_SIZE  = 50'000'000;

auto const N_RANGE = nvbench::range(10'000'000, 100'000'000, 20'000'000);
auto const N_RANGE_CACHE =
  std::vector<nvbench::int64_t>{8'000, 80'000, 800'000, 8'000'000, 80'000'000};
auto const OCCUPANCY_RANGE     = nvbench::range(0.1, 0.9, 0.1);
auto const MULTIPLICITY_RANGE  = std::vector<nvbench::int64_t>{1, 2, 4, 8, 16};
auto const MATCHING_RATE_RANGE = nvbench::range(0.1, 1., 0.1);
auto const SKEW_RANGE          = nvbench::range(0.1, 1., 0.1);

}  // namespace cuco::benchmark::defaults
