/*
 * Copyright (c) 2023-2025, NVIDIA CORPORATION.
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

#include <cuco/detail/storage/storage.cuh>

#include <cuda/std/cstdint>

namespace cuco {
// Forward declaration to avoid circular dependency
template <typename T, int32_t BucketSize, typename Extent, typename Allocator>
class bucket_storage;
}  // namespace cuco

namespace cuco {

/**
 * @brief Public storage class.
 *
 * @note This is a public interface used to control storage bucket size. A bucket consists of one
 * or multiple contiguous slots. The bucket size defines the workload granularity for each CUDA
 * thread, i.e., how many slots a thread would concurrently operate on when performing modify or
 * lookup operations. cuCollections uses the array of bucket storage to supersede the raw flat slot
 * storage due to its superior granularity control: When bucket size equals one, array of buckets
 * performs the same as the flat storage. If the underlying operation is more memory bandwidth
 * bound, e.g., high occupancy multimap operations, a larger bucket size can reduce the length of
 * probing sequences thus improve runtime performance.
 *
 * @tparam BucketSize Number of elements per bucket storage
 */
template <int BucketSize>
class storage {
 public:
  /// Number of slots per bucket storage
  static constexpr cuda::std::int32_t bucket_size = BucketSize;

  /// Type of implementation details
  template <class T, class Extent, class Allocator>
  using impl = bucket_storage<T, BucketSize, Extent, Allocator>;
};
}  // namespace cuco
