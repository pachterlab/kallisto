/*
 * Copyright (c) 2023, NVIDIA CORPORATION.
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
 */

#pragma once

#include <cuda/std/type_traits>

namespace cuco {
namespace detail {

/**
 * @brief Ceiling of an integer division
 *
 * @tparam T Type of dividend
 * @tparam U Type of divisor
 *
 * @throw If `T` is not an integral type
 * @throw If `U` is not an integral type
 *
 * @param dividend Numerator
 * @param divisor Denominator
 *
 * @return Ceiling of the integer division
 */
template <typename T, typename U>
__host__ __device__ constexpr T int_div_ceil(T dividend, U divisor) noexcept
{
  static_assert(cuda::std::is_integral_v<T>);
  static_assert(cuda::std::is_integral_v<U>);
  return (dividend + divisor - 1) / divisor;
}

}  // namespace detail
}  // namespace cuco
