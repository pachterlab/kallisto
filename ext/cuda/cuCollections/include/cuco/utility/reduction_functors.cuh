/*
 * Copyright (c) 2024, NVIDIA CORPORATION.
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

namespace cuco::reduce {

/**
 * @brief Device functor performing sum reduction, used with `insert-or-apply`
 */
struct plus {
  /**
   * @brief Performs atomic fetch_add on payload and the new value to be inserted
   *
   * @tparam T The payload type
   * @tparam Scope The cuda::thread_scope used for atomic_ref
   *
   * @param payload_ref The atomic_ref pointing to payload part of the slot
   * @param val The new value to be applied as reduction to the current value
   * in the payload.
   */
  template <typename T, cuda::thread_scope Scope>
  __device__ void operator()(cuda::atomic_ref<T, Scope> payload_ref, const T& val)
  {
    payload_ref.fetch_add(val, cuda::memory_order_relaxed);
  }
};

/**
 * @brief Device functor performing max reduction, used with `insert-or-apply`
 */
struct max {
  /**
   * @brief Performs atomic fetch_max on payload and the new value to be inserted
   *
   * @tparam T The payload type
   * @tparam Scope The cuda::thread_scope used for atomic_ref
   *
   * @param payload_ref The atomic_ref pointing to payload part of the slot
   * @param val The new value to be applied as reduction to the current value
   * in the payload.
   */
  template <typename T, cuda::thread_scope Scope>
  __device__ void operator()(cuda::atomic_ref<T, Scope> payload_ref, const T& val)
  {
    payload_ref.fetch_max(val, cuda::memory_order_relaxed);
  }
};

/**
 * @brief Device functor performing min reduction, used with `insert-or-apply`
 */
struct min {
  /**
   * @brief Performs atomic fetch_min on payload and the new value to be inserted
   *
   * @tparam T The payload type
   * @tparam Scope The cuda::thread_scope used for atomic_ref
   *
   * @param payload_ref The atomic_ref pointing to payload part of the slot
   * @param val The new value to be applied as reduction to the current value
   * in the payload.
   */
  template <typename T, cuda::thread_scope Scope>
  __device__ void operator()(cuda::atomic_ref<T, Scope> payload_ref, const T& val)
  {
    payload_ref.fetch_min(val, cuda::memory_order_relaxed);
  }
};

}  // namespace cuco::reduce