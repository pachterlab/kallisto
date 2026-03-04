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

#include <cuco/detail/error.hpp>
#include <cuco/detail/prime.hpp>  // TODO move to detail/extent/
#include <cuco/detail/utility/math.cuh>
#include <cuco/detail/utils.hpp>
#include <cuco/probing_scheme.cuh>
#include <cuco/storage.cuh>
#include <cuco/utility/fast_int.cuh>

#include <cuda/std/type_traits>

#include <algorithm>
#include <cstdint>
#include <limits>

namespace cuco {
template <typename SizeType, std::size_t N>
struct valid_extent {
  using value_type = SizeType;  ///< Extent value type

  __host__ __device__ constexpr value_type value() const noexcept { return N; }
  __host__ __device__ explicit constexpr operator value_type() const noexcept { return value(); }

 private:
  __host__ __device__ explicit constexpr valid_extent() noexcept {}
  __host__ __device__ explicit constexpr valid_extent(SizeType) noexcept {}

  // Friend declarations for all make_valid_extent overloads
  template <int32_t CGSize_, int32_t BucketSize_, typename SizeType_, std::size_t N_>
  friend auto constexpr make_valid_extent(extent<SizeType_, N_> ext);

  template <typename ProbingScheme, typename Storage, typename SizeType_, std::size_t N_>
  friend auto constexpr make_valid_extent(extent<SizeType_, N_> ext);

  template <template <typename> class ProbingScheme,
            typename Storage,
            typename SizeType_,
            std::size_t N_>
  friend auto constexpr make_valid_extent(extent<SizeType_, N_> ext);

  template <template <typename, typename> class ProbingScheme,
            typename Storage,
            typename SizeType_,
            std::size_t N_>
  friend auto constexpr make_valid_extent(extent<SizeType_, N_> ext);

  // Operator overloads
  template <typename Rhs>
  friend __host__ __device__ constexpr value_type operator-(valid_extent const& lhs,
                                                            Rhs rhs) noexcept
  {
    return lhs.value() - rhs;
  }

  template <typename Rhs>
  friend __host__ __device__ constexpr value_type operator/(valid_extent const& lhs,
                                                            Rhs rhs) noexcept
  {
    return lhs.value() / rhs;
  }

  template <typename Lhs>
  friend __host__ __device__ constexpr value_type operator%(Lhs lhs,
                                                            valid_extent const& rhs) noexcept
  {
    return lhs % rhs.value();
  }
};

template <typename SizeType>
struct valid_extent<SizeType, dynamic_extent> : cuco::utility::fast_int<SizeType> {
  using value_type =
    typename cuco::utility::fast_int<SizeType>::fast_int::value_type;  ///< Extent value type

 private:
  using cuco::utility::fast_int<SizeType>::fast_int;

  // Friend declarations for all make_valid_extent overloads
  template <int32_t CGSize_, int32_t BucketSize_, typename SizeType_, std::size_t N_>
  friend auto constexpr make_valid_extent(extent<SizeType_, N_> ext);

  template <typename ProbingScheme, typename Storage, typename SizeType_, std::size_t N_>
  friend auto constexpr make_valid_extent(extent<SizeType_, N_> ext);

  template <template <typename> class ProbingScheme,
            typename Storage,
            typename SizeType_,
            std::size_t N_>
  friend auto constexpr make_valid_extent(extent<SizeType_, N_> ext);

  template <template <typename, typename> class ProbingScheme,
            typename Storage,
            typename SizeType_,
            std::size_t N_>
  friend auto constexpr make_valid_extent(extent<SizeType_, N_> ext);
};

// Primary implementation for fixed CGSize and BucketSize
template <int32_t CGSize, int32_t BucketSize, typename SizeType, std::size_t N>
[[nodiscard]] auto constexpr make_valid_extent(extent<SizeType, N> ext)
{
  auto constexpr stride    = CGSize * BucketSize;
  auto constexpr max_prime = cuco::detail::primes.back();
  auto constexpr max_value =
    (static_cast<uint64_t>(cuda::std::numeric_limits<SizeType>::max()) < max_prime)
      ? cuda::std::numeric_limits<SizeType>::max()
      : static_cast<SizeType>(max_prime);
  auto const size = cuco::detail::int_div_ceil(
    cuda::std::max(static_cast<SizeType>(ext), static_cast<SizeType>(1)), stride);
  if (size > max_value) { CUCO_FAIL("Invalid input extent"); }

  if constexpr (N == dynamic_extent) {
    return valid_extent<SizeType, dynamic_extent>{static_cast<SizeType>(
      *cuco::detail::lower_bound(
        cuco::detail::primes.begin(), cuco::detail::primes.end(), static_cast<uint64_t>(size)) *
      stride)};
  } else {
    return valid_extent<SizeType,
                        static_cast<std::size_t>(
                          *cuco::detail::lower_bound(cuco::detail::primes.begin(),
                                                     cuco::detail::primes.end(),
                                                     static_cast<uint64_t>(size)) *
                          stride)>{};
  }
}

// Overload for SizeType without extent
template <int32_t CGSize, int32_t BucketSize, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType size)
{
  return make_valid_extent<CGSize, BucketSize, SizeType, dynamic_extent>(extent<SizeType>{size});
}

// Implementation for ProbingScheme and Storage types
template <typename ProbingScheme, typename Storage, typename SizeType, std::size_t N>
[[nodiscard]] auto constexpr make_valid_extent(extent<SizeType, N> ext)
{
  if constexpr (cuco::is_double_hashing<ProbingScheme>::value) {
    return make_valid_extent<ProbingScheme::cg_size, Storage::bucket_size, SizeType, N>(ext);
  } else {
    auto constexpr stride = ProbingScheme::cg_size * Storage::bucket_size;
    auto const size =
      cuco::detail::int_div_ceil(
        cuda::std::max(static_cast<SizeType>(ext), static_cast<SizeType>(1)), stride) +
      static_cast<SizeType>(ext == 0);

    if constexpr (N == dynamic_extent) {
      return valid_extent<SizeType, dynamic_extent>{size * stride};
    } else {
      return valid_extent<SizeType, size * stride>{};
    }
  }
}

// Overload for ProbingScheme and Storage with SizeType
template <typename ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(extent<SizeType> ext, double desired_load_factor)
{
  CUCO_EXPECTS(desired_load_factor > 0., "Desired occupancy must be larger than zero");
  CUCO_EXPECTS(desired_load_factor <= 1., "Desired occupancy must be no larger than one");

  auto const temp = cuda::std::ceil(static_cast<double>(SizeType{ext}) / desired_load_factor);
  if (temp > static_cast<double>(cuda::std::numeric_limits<SizeType>::max())) {
    CUCO_FAIL(
      "Invalid load factor: requested extent divided by load factor exceeds maximum representable "
      "value");
  }
  return make_valid_extent<ProbingScheme, Storage>(
    cuco::extent<SizeType>{static_cast<SizeType>(temp)});
}

template <typename ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType size, double desired_load_factor)
{
  return make_valid_extent<ProbingScheme, Storage>(cuco::extent<SizeType>{size},
                                                   desired_load_factor);
}

template <typename ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType size)
{
  return make_valid_extent<ProbingScheme, Storage, SizeType, dynamic_extent>(
    cuco::extent<SizeType>{size});
}

// Template template parameter overloads for single-type ProbingScheme
template <template <typename> class ProbingScheme,
          typename Storage,
          typename SizeType,
          std::size_t N>
[[nodiscard]] auto constexpr make_valid_extent(extent<SizeType, N> ext)
{
  using ProbeType = ProbingScheme<int>;
  return make_valid_extent<ProbeType, Storage, SizeType, N>(ext);
}

template <template <typename> class ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType size)
{
  using ProbeType = ProbingScheme<int>;
  return make_valid_extent<ProbeType, Storage, SizeType>(size);
}

// Template template parameter overloads for two-type ProbingScheme
template <template <typename, typename> class ProbingScheme,
          typename Storage,
          typename SizeType,
          std::size_t N>
[[nodiscard]] auto constexpr make_valid_extent(extent<SizeType, N> ext)
{
  using ProbeType = ProbingScheme<int, int>;
  return make_valid_extent<ProbeType, Storage, SizeType, N>(ext);
}

template <template <typename, typename> class ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType size)
{
  using ProbeType = ProbingScheme<int, int>;
  return make_valid_extent<ProbeType, Storage, SizeType>(size);
}

}  // namespace cuco
