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

#pragma once

#include <cuda/atomic>

namespace cuco {
namespace test {
namespace detail {

template <typename Iterator, typename Predicate>
__global__ void count_if(Iterator begin,
                         Iterator end,
                         cuda::atomic<int, cuda::thread_scope_device>* count,
                         Predicate p)
{
  auto tid = blockDim.x * blockIdx.x + threadIdx.x;
  auto it  = begin + tid;

  while (it < end) {
    count->fetch_add(static_cast<int>(p(*it)));
    it += gridDim.x * blockDim.x;
  }
}

template <typename Iterator1, typename Iterator2, typename Predicate>
__global__ void count_if(Iterator1 begin1,
                         Iterator1 end1,
                         Iterator2 begin2,
                         cuda::atomic<int, cuda::thread_scope_device>* count,
                         Predicate p)
{
  auto const n = end1 - begin1;
  auto tid     = blockDim.x * blockIdx.x + threadIdx.x;

  while (tid < n) {
    auto cmp = begin1 + tid;
    auto ref = begin2 + tid;
    count->fetch_add(static_cast<int>(p(*cmp, *ref)));
    tid += gridDim.x * blockDim.x;
  }
}

}  // namespace detail
}  // namespace test
}  // namespace cuco
