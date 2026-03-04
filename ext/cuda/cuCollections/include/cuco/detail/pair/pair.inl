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

#include <cuda/std/type_traits>
#include <cuda/std/utility>

namespace cuco {

template <typename First, typename Second>
__host__ __device__ constexpr pair<First, Second>::pair(First const& f, Second const& s)
  : first{f}, second{s}
{
}

template <typename First, typename Second>
template <typename F, typename S>
__host__ __device__ constexpr pair<First, Second>::pair(pair<F, S> const& p)
  : first{p.first}, second{p.second}
{
}

template <typename F, typename S>
__host__ __device__ constexpr pair<cuda::std::decay_t<F>, cuda::std::decay_t<S>> make_pair(
  F&& f, S&& s) noexcept
{
  return pair<cuda::std::decay_t<F>, cuda::std::decay_t<S>>(cuda::std::forward<F>(f),
                                                            cuda::std::forward<S>(s));
}

template <class T1, class T2, class U1, class U2>
__host__ __device__ constexpr bool operator==(cuco::pair<T1, T2> const& lhs,
                                              cuco::pair<U1, U2> const& rhs) noexcept
{
  return lhs.first == rhs.first and lhs.second == rhs.second;
}

}  // namespace cuco

namespace cuda::std {
#include <cuco/detail/pair/tuple_helpers.inl>
}  // namespace cuda::std
