/*
 * Copyright (c) 2020-2025, NVIDIA CORPORATION.
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

#include <cuco/detail/pair/helpers.cuh>
#include <cuco/detail/pair/traits.hpp>

#include <cuda/std/tuple>
#include <cuda/std/type_traits>
#include <thrust/device_reference.h>

#include <tuple>

namespace cuco {

/**
 * @brief Custom pair type
 *
 * @note This is necessary because `cuda::std::pair` is under aligned.
 *
 * @tparam First Type of the first value in the pair
 * @tparam Second Type of the second value in the pair
 */
template <typename First, typename Second>
struct alignas(detail::pair_alignment<First, Second>()) pair {
  using first_type  = First;   ///< Type of the first value in the pair
  using second_type = Second;  ///< Type of the second value in the pair

  pair()            = default;
  ~pair()           = default;
  pair(pair const&) = default;  ///< Copy constructor
  pair(pair&&)      = default;  ///< Move constructor

  /**
   * @brief Replaces the contents of the pair with another pair.
   *
   * @return Reference of the current pair object
   */
  pair& operator=(pair const&) = default;

  /**
   * @brief Replaces the contents of the pair with another pair.
   *
   * @return Reference of the current pair object
   */
  pair& operator=(pair&&) = default;

  /**
   * @brief Constructs a pair from objects `f` and `s`.
   *
   * @param f The object to copy into `first`
   * @param s The object to copy into `second`
   */
  __host__ __device__ constexpr pair(First const& f, Second const& s);

  /**
   * @brief Constructs a pair by copying from the given pair `p`.
   *
   * @tparam F Type of the first value of `p`
   * @tparam S Type of the second value of `p`
   *
   * @param p The pair to copy from
   */
  template <typename F, typename S>
  __host__ __device__ constexpr pair(pair<F, S> const& p);

  /**
   * @brief Constructs a pair from the given std::pair-like `p`.
   *
   * @tparam T Type of the pair to copy from
   *
   * @param p The input pair to copy from
   */
  template <typename T, cuda::std::enable_if_t<detail::is_std_pair_like<T>::value>* = nullptr>
  __host__ __device__ constexpr pair(T const& p)
    : pair{std::get<0>(thrust::raw_reference_cast(p)), std::get<1>(thrust::raw_reference_cast(p))}
  {
  }

  /**
   * @brief Constructs a pair from the given cuda::std::pair-like `p`.
   *
   * @tparam T Type of the pair to copy from
   *
   * @param p The input pair to copy from
   */
  template <typename T, cuda::std::enable_if_t<detail::is_cuda_std_pair_like<T>::value>* = nullptr>
  __host__ __device__ constexpr pair(T const& p)
    : pair{cuda::std::get<0>(thrust::raw_reference_cast(p)),
           cuda::std::get<1>(thrust::raw_reference_cast(p))}
  {
  }

  First first;    ///< The first value in the pair
  Second second;  ///< The second value in the pair
};

/**
 * @brief Creates a pair with the given first and second elements
 *
 * @tparam F Type of first element
 * @tparam S Type of second element
 *
 * @param f First element
 * @param s Second element
 *
 * @return A pair with first element `f` and second element `s`.
 */
template <typename F, typename S>
__host__ __device__ constexpr pair<cuda::std::decay_t<F>, cuda::std::decay_t<S>> make_pair(
  F&& f, S&& s) noexcept;

/**
 * @brief Tests if both elements of lhs and rhs are equal
 *
 * @tparam T1 Type of the first element of the left-hand side pair
 * @tparam T2 Type of the second element of the left-hand side pair
 * @tparam U1 Type of the first element of the right-hand side pair
 * @tparam U2 Type of the second element of the right-hand side pair
 *
 * @param lhs Left-hand side pair
 * @param rhs Right-hand side pair
 *
 * @return True if two pairs are equal. False otherwise
 */
template <class T1, class T2, class U1, class U2>
__host__ __device__ constexpr bool operator==(cuco::pair<T1, T2> const& lhs,
                                              cuco::pair<U1, U2> const& rhs) noexcept;

template <typename T>
struct is_cuco_pair : cuda::std::false_type {};

template <typename F, typename S>
struct is_cuco_pair<cuco::pair<F, S>> : cuda::std::true_type {};

/// Trait determing whether a given type is tuple-like or not
template <typename T>
struct is_tuple_like : cuda::std::disjunction<is_cuco_pair<T>,
                                              detail::is_std_pair_like<T>,
                                              detail::is_cuda_std_pair_like<T>> {};

}  // namespace cuco

#include <cuco/detail/pair/pair.inl>
