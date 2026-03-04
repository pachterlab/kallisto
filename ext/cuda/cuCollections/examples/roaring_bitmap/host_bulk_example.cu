/*
 * Copyright (c) 2025 NVIDIA CORPORATION.
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

#include <cuco/roaring_bitmap.cuh>
#include <cuco/utility/traits.hpp>

#include <cuda/std/cstddef>
#include <cuda/std/cstdint>
#include <cuda/std/type_traits>
#include <thrust/device_vector.h>
#include <thrust/logical.h>
#include <thrust/universal_vector.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/**
 * @file host_bulk_example.cu
 * @brief Demonstrates usage of the roaring_bitmap "bulk" lookup host APIs.
 *
 * In this example we load two 32-bit bitmaps and one 64-bit bitmap (portable format) from the
 * [RoaringBitmapFormatSpec](https://github.com/RoaringBitmap/RoaringFormatSpec) repository and
 * check if the bulk lookup API returns the correct results. Namely, we test the following files:
 * - [bitmapwithoutruns.bin
 * (32-bit)](https://github.com/RoaringBitmap/RoaringFormatSpec/blob/5177ad9/testdata/bitmapwithoutruns.bin)
 * - [bitmapwithruns.bin
 * (32-bit)](https://github.com/RoaringBitmap/RoaringFormatSpec/blob/5177ad9/testdata/bitmapwithruns.bin)
 * - [portable_bitmap64.bin
 * (64-bit)](https://github.com/RoaringBitmap/RoaringFormatSpec/blob/5177ad9/testdata64/portable_bitmap64.bin)
 *
 * @note This example requires the cmake option -DCUCO_DOWNLOAD_ROARING_TESTDATA=ON to be set.
 *
 */
template <typename KeyType>
bool check(std::string const& bitmap_file_path)
{
  auto generate_keys = []() -> thrust::device_vector<KeyType> {
    if constexpr (cuda::std::is_same_v<KeyType, cuda::std::uint32_t>) {
      // Create query keys for the bitmapwith{out}runs.bin files:
      // https://github.com/RoaringBitmap/RoaringFormatSpec/blob/5177ad9/testdata/README.md#test-data
      std::vector<cuda::std::uint32_t> keys;
      for (cuda::std::uint32_t k = 0; k < 100000; k += 1000) {
        keys.push_back(k);
      }
      for (int k = 100000; k < 200000; ++k) {
        keys.push_back(3 * k);
      }
      for (int k = 700000; k < 800000; ++k) {
        keys.push_back(k);
      }
      return thrust::device_vector<cuda::std::uint32_t>(keys.begin(), keys.end());
    } else if constexpr (cuda::std::is_same_v<KeyType, cuda::std::uint64_t>) {
      // Create query keys for the portable_bitmap64.bin file:
      // https://github.com/RoaringBitmap/RoaringFormatSpec/blob/5177ad9/testdata64/README.md#portable_bitmap64bin
      std::vector<cuda::std::uint64_t> keys;
      for (cuda::std::uint64_t k = 0x00000ull; k < 0x09000ull; ++k) {
        keys.push_back(k);
      }
      for (cuda::std::uint64_t k = 0x0A000ull; k < 0x10000ull; ++k) {
        keys.push_back(k);
      }
      keys.push_back(0x20000ull);
      keys.push_back(0x20005ull);
      for (cuda::std::uint64_t i = 0; i < 0x10000ull; i += 2ull) {
        keys.push_back(0x80000ull + i);
      }
      return thrust::device_vector<cuda::std::uint64_t>(keys.begin(), keys.end());
    } else {
      static_assert(cuco::dependent_false<KeyType>, "KeyType must be uint32_t or uint64_t");
      return {};
    }
  };

  // Open file
  std::ifstream file(bitmap_file_path, std::ios::binary);
  if (!file.is_open()) {
    std::cerr << "Failed to open " << bitmap_file_path << std::endl;
    return false;
  }

  // Get file size
  auto file_size = std::filesystem::file_size(bitmap_file_path);

  // Allocate host memory for the bitmap file
  thrust::universal_host_pinned_vector<cuda::std::byte> buffer(file_size);

  // Read file into memory
  file.read(reinterpret_cast<char*>(thrust::raw_pointer_cast(buffer.data())), file_size);
  file.close();

  // Create roaring bitmap from the file
  cuco::experimental::roaring_bitmap<KeyType> roaring_bitmap(
    thrust::raw_pointer_cast(buffer.data()));

  // Generate query keys (all should be contained in the bitmap)
  auto keys = generate_keys();

  // Create a vector to store the results
  thrust::device_vector<bool> contained(keys.size(), false);

  // Bulk-lookup query keys against the bitmap
  roaring_bitmap.contains(keys.begin(), keys.end(), contained.begin());

  // Check if all the keys are contained in the bitmap
  bool all_contained = thrust::all_of(contained.begin(), contained.end(), ::cuda::std::identity{});
  return all_contained;
}

int main()
{
#ifdef CUCO_ROARING_DATA_DIR
  std::string const data_dir = CUCO_ROARING_DATA_DIR;
  bool success               = check<cuda::std::uint32_t>(data_dir + "/bitmapwithoutruns.bin");
  success &= check<cuda::std::uint32_t>(data_dir + "/bitmapwithruns.bin");
  success &= check<cuda::std::uint64_t>(data_dir + "/portable_bitmap64.bin");

  std::cout << "success: " << std::boolalpha << success << std::endl;

  return success ? 0 : 1;
#else
  std::cerr << "This example requires CUCO_ROARING_DATA_DIR to be defined (build with cmake option "
               "-DCUCO_DOWNLOAD_ROARING_TESTDATA=ON)"
            << std::endl;
  return 1;
#endif
}
