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
 * limitations under the License.
 */

#pragma once

#include <cuco/detail/__config>

#include <cuda/std/bit>
#include <cuda/std/cstdint>
#include <cuda/std/limits>
#include <cuda/std/type_traits>

namespace cuco::utility {

/**
 * @brief Integer type with optimized division and modulo operators.
 *
 * @tparam T Underlying integer type
 */
template <typename T>
struct fast_int {
  static_assert(cuda::std::is_same_v<T, cuda::std::int32_t> or
                  cuda::std::is_same_v<T, cuda::std::uint32_t>
#if defined(CUCO_HAS_INT128)
                  or cuda::std::is_same_v<T, cuda::std::int64_t> or
                  cuda::std::is_same_v<T, cuda::std::uint64_t>
#endif
                ,
                "Unsupported integer type");

  using value_type = T;  ///< Underlying integer type

  /**
   * @brief Constructs a fast_int from an integer value.
   *
   * @param value Integer value
   */
  __host__ __device__ explicit constexpr fast_int(T value) noexcept : value_{value}
  {
    evaluate_magic_numbers();
  }

  /**
   * @brief Get the underlying integer value.
   *
   * @return Underlying value
   */
  __host__ __device__ constexpr value_type value() const noexcept { return value_; }

  /**
   * @brief Explicit conversion operator to the underlying value type.
   *
   * @return Underlying value
   */
  __host__ __device__ explicit constexpr operator value_type() const noexcept { return value_; }

 private:
  using intermediate_type =
    cuda::std::conditional_t<sizeof(value_type) == 4,
                             cuda::std::uint64_t,
                             unsigned __int128>;  ///< Intermediate type for multiplication
  using unsigned_value_type = cuda::std::make_unsigned_t<value_type>;  ///< Unsigned value type
  using signed_value_type   = cuda::std::make_signed_t<value_type>;    ///< Signed value type

  static constexpr value_type value_bits =
    CHAR_BIT * sizeof(value_type);  ///< Number of bits required to represent the value

  /**
   * @brief Computes the high bits of the multiplication of two unsigned integers.
   *
   * @param lhs Left-hand side of the multiplication
   * @param rhs Right-hand side of the multiplication
   *
   * @return High bits of the multiplication
   */
  __host__ __device__ constexpr value_type mulhi(unsigned_value_type lhs,
                                                 unsigned_value_type rhs) const noexcept
  {
#if defined(__CUDA_ARCH__)
    if constexpr (sizeof(value_type) == 4) {
      return __umulhi(lhs, rhs);
    } else {
      return __umul64hi(lhs, rhs);
    }
#else
    return (intermediate_type(lhs) * intermediate_type(rhs)) >> value_bits;
#endif
  }

  /**
   * @brief Computes the log2 of an unsigned integer.
   *
   * @param v Unsigned integer
   *
   * @return Log2 of the unsigned integer
   */
  __host__ __device__ constexpr value_type log2(value_type v) const noexcept
  {
    return cuda::std::bit_width(unsigned_value_type(v)) - 1;
  }

  /**
   * @brief Computes the magic numbers for the fast division.
   */
  __host__ __device__ constexpr void evaluate_magic_numbers() noexcept
  {
    // TODO assert(value_ > 0);
    auto const val_log2 = this->log2(value_);

    // if value_ is a power of 2, we can use a simple shift
    if (cuda::std::has_single_bit(unsigned_value_type(value_))) {
      magic_ = 0;
      shift_ = val_log2;
    } else {
      auto upper      = intermediate_type(1) << value_bits;
      auto lower      = intermediate_type(1);
      auto const lval = intermediate_type(value_);

      // compute the magic number and shift; see "Hacker's Delight" by Henry S. Warren, Jr., 10-2
      for (shift_ = 0; shift_ < val_log2; ++shift_, upper <<= 1, lower <<= 1) {
        if ((upper % lval) <= lower) { break; }
      }
      magic_ = upper / lval;
    }
  }

  value_type value_;  ///< Underlying integer value
  value_type magic_;  ///< Magic number for fast division
  value_type shift_;  ///< Shift for fast division

  template <typename Lhs>
  friend __host__ __device__ constexpr value_type operator/(Lhs lhs, fast_int const& rhs) noexcept
  {
    static_assert(cuda::std::is_same_v<Lhs, value_type>,
                  "Left-hand side operand must be of type value_type.");
    if (rhs.value_ == 1) { return lhs; }                // edge case for value_ == 1
    if (rhs.magic_ == 0) { return lhs >> rhs.shift_; }  // edge case for value_ == pow2
    auto const mul = (lhs == cuda::std::numeric_limits<T>::max()) ? lhs : lhs + 1;
    return rhs.mulhi(rhs.magic_, mul) >> rhs.shift_;
  }

  template <typename Rhs>
  friend __host__ __device__ constexpr auto operator-(fast_int const& lhs, Rhs rhs) noexcept
  {
    return lhs.value() - rhs;
  }

  template <typename Rhs>
  friend __host__ __device__ constexpr auto operator/(fast_int const& lhs, Rhs rhs) noexcept
  {
    return lhs.value() / rhs;
  }

  template <typename Lhs>
  friend __host__ __device__ constexpr value_type operator%(Lhs lhs, fast_int const& rhs) noexcept
  {
    return lhs - (lhs / rhs) * rhs.value_;
  }
};
}  // namespace cuco::utility
