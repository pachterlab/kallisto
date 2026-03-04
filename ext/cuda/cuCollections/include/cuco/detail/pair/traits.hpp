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
 */

#pragma once

#include <cuda/std/tuple>
#include <cuda/std/type_traits>
#include <thrust/device_reference.h>

#include <tuple>

namespace cuco::detail {

template <typename T, typename = void>
struct is_std_pair_like : cuda::std::false_type {};

template <typename T>
struct is_std_pair_like<T,
                        cuda::std::void_t<decltype(std::get<0>(cuda::std::declval<T>())),
                                          decltype(std::get<1>(cuda::std::declval<T>()))>>
  : cuda::std::
      conditional_t<std::tuple_size<T>::value == 2, cuda::std::true_type, cuda::std::false_type> {};

template <typename T, typename = void>
struct is_cuda_std_pair_like_impl : cuda::std::false_type {};

template <typename T>
struct is_cuda_std_pair_like_impl<
  T,
  cuda::std::void_t<decltype(cuda::std::get<0>(cuda::std::declval<T>())),
                    decltype(cuda::std::get<1>(cuda::std::declval<T>())),
                    decltype(cuda::std::tuple_size<T>::value)>>
  : cuda::std::conditional_t<cuda::std::tuple_size<T>::value == 2,
                             cuda::std::true_type,
                             cuda::std::false_type> {};

template <typename T>
struct is_cuda_std_pair_like
  : is_cuda_std_pair_like_impl<cuda::std::remove_reference_t<decltype(thrust::raw_reference_cast(
      cuda::std::declval<T>()))>> {};

}  // namespace cuco::detail
