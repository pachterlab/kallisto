/*
 * Copyright (c) 2024-2025, NVIDIA CORPORATION.
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

#include <cuda/std/array>

#include <vector>

namespace cuco::benchmark::defaults {

using BF_KEY  = nvbench::int64_t;
using BF_HASH = cuco::xxhash_64<char>;
using BF_WORD = nvbench::uint32_t;

static constexpr auto BF_N               = 1'000'000'000;
static constexpr auto BF_SIZE_MB         = 2'000;
static constexpr auto BF_WORDS_PER_BLOCK = 8;

auto const BF_SIZE_MB_RANGE_CACHE =
  std::vector<nvbench::int64_t>{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048};
auto const BF_PATTERN_BITS_RANGE = std::vector<nvbench::int64_t>{1, 2, 4, 6, 8, 16};

}  // namespace cuco::benchmark::defaults
