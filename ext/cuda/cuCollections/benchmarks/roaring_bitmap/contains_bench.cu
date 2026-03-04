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

#include <benchmark_defaults.hpp>
#include <benchmark_utils.hpp>

#include <cuco/roaring_bitmap.cuh>
#include <cuco/utility/key_generator.cuh>

#include <nvbench/nvbench.cuh>

#include <cuda/std/cstddef>
#include <cuda/std/cstdint>
#include <thrust/device_vector.h>
#include <thrust/universal_vector.h>

#include <filesystem>
#include <fstream>
#include <string>

using namespace cuco::benchmark;  // defaults
using namespace cuco::utility;    // key_generator, distribution

template <typename T>
void roaring_bitmap_contains(nvbench::state& state, nvbench::type_list<T>)
{
  auto const num_items   = state.get_int64("NumInputs");
  auto const bitmap_file = state.get_string_or_default("BitmapFile", {});

  std::ifstream file(bitmap_file, std::ios::binary);
  if (!file.is_open()) { state.skip("Bitmap file not found"); }

  // Get file size
  auto const file_size = std::filesystem::file_size(bitmap_file);

  thrust::universal_host_pinned_vector<cuda::std::byte> buffer(file_size);

  file.read(reinterpret_cast<char*>(thrust::raw_pointer_cast(buffer.data())), file_size);
  file.close();

  cuco::experimental::roaring_bitmap<T> roaring_bitmap(thrust::raw_pointer_cast(buffer.data()));

  thrust::device_vector<T> items(num_items);

  key_generator gen{};
  gen.generate(distribution::unique{}, items.begin(), items.end());

  thrust::device_vector<bool> contained(items.size(), false);

  state.add_element_count(items.size());
  state.add_global_memory_reads<T>(items.size(), "InputSize");

  auto& summ = state.add_summary("BitmapSizeMB");
  summ.set_string("hint", "BitmapSize");
  summ.set_string("short_name", "BitmapSizeMB");
  summ.set_string("description", "Bitmap size in MB");
  summ.set_float64("value", static_cast<double>(file_size) / (1024 * 1024));

  state.exec([&](nvbench::launch& launch) {
    roaring_bitmap.contains_async(
      items.begin(), items.end(), contained.begin(), {launch.get_stream()});
  });
}

NVBENCH_BENCH_TYPES(roaring_bitmap_contains,
                    NVBENCH_TYPE_AXES(nvbench::type_list<nvbench::uint32_t>))
  .set_name("roaring_bitmap_contains")
// Default benchmark is only available if the Roaring bitmap testdata has been downloaded
#ifdef CUCO_ROARING_DATA_DIR
  .add_string_axis("BitmapFile", {std::string(CUCO_ROARING_DATA_DIR) + "/bitmapwithruns.bin"})
#endif
  .add_int64_power_of_two_axis("NumInputs", {32});

NVBENCH_BENCH_TYPES(roaring_bitmap_contains,
                    NVBENCH_TYPE_AXES(nvbench::type_list<nvbench::uint64_t>))
  .set_name("roaring_bitmap_contains")
// Default benchmark is only available if the Roaring bitmap testdata has been downloaded
#ifdef CUCO_ROARING_DATA_DIR
  .add_string_axis("BitmapFile", {std::string(CUCO_ROARING_DATA_DIR) + "/portable_bitmap64.bin"})
#endif
  .add_int64_power_of_two_axis("NumInputs", {31});
