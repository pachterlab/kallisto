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

#include <cuco/detail/hyperloglog/tuning.cuh>

#include <cuda/functional>
#include <cuda/std/cmath>
#include <cuda/std/limits>

#include <cstddef>

namespace cuco::hyperloglog_ns::detail {

/**
 * @brief Estimate correction algorithm based on HyperLogLog++.
 *
 * @note Variable names correspond to the definitions given in the HLL++ paper:
 * https://static.googleusercontent.com/media/research.google.com/de//pubs/archive/40671.pdf
 * @note Previcion must be >= 4.
 *
 */
class finalizer {
  // Note: Most of the types in this implementation are explicit instead of relying on `auto` to
  // avoid confusion with the reference implementation.

 public:
  /**
   * @brief Contructs an HLL finalizer object.
   *
   * @param precision HLL precision parameter
   */
  __host__ __device__ constexpr finalizer(int precision) : precision_{precision}, m_{1 << precision}
  {
  }

  /**
   * @brief Compute the bias-corrected cardinality estimate.
   *
   * @param z Geometric mean of registers
   * @param v Number of 0 registers
   *
   * @return Bias-corrected cardinality estimate
   */
  __host__ __device__ std::size_t operator()(double z, int v) const noexcept
  {
    double e = this->alpha_mm() / z;

    if (v > 0) {
      // Use linear counting for small cardinality estimates.
      double const h = this->m_ * log(static_cast<double>(this->m_) / v);
      // The threshold `2.5 * m` is from the original HLL algorithm.
      if (e <= 2.5 * this->m_) { return cuda::std::round(h); }

      if (this->precision_ < 19) {
        e = (h <= threshold(this->precision_)) ? h : this->bias_corrected_estimate(e);
      }
    } else {
      // HLL++ is defined only when p < 19, otherwise we need to fallback to HLL.
      if (this->precision_ < 19) { e = this->bias_corrected_estimate(e); }
    }

    return cuda::std::round(e);
  }

 private:
  __host__ __device__ constexpr double alpha_mm() const noexcept
  {
    switch (this->m_) {
      case 16: return 0.673 * this->m_ * this->m_;
      case 32: return 0.697 * this->m_ * this->m_;
      case 64: return 0.709 * this->m_ * this->m_;
      default: return (0.7213 / (1.0 + 1.079 / this->m_)) * this->m_ * this->m_;
    }
  }

  __host__ __device__ constexpr double bias_corrected_estimate(double e) const noexcept
  {
    return (e < 5.0 * this->m_) ? e - this->bias(e) : e;
  }

  __host__ __device__ constexpr double bias(double e) const noexcept
  {
    auto const anchor_index = this->interpolation_anchor_index(e);
    int const n             = raw_estimate_data_size(this->precision_);

    auto low  = cuda::std::max(anchor_index - k + 1, 0);
    auto high = cuda::std::min(low + k, n);
    // Keep moving bounds as long as the (exclusive) high bound is closer to the estimate than
    // the lower (inclusive) bound.
    while (high < n and this->distance(e, high) < this->distance(e, low)) {
      low += 1;
      high += 1;
    }

    auto biases     = bias_data(this->precision_);
    double bias_sum = 0.0;
    for (int i = low; i < high; ++i) {
      bias_sum += biases[i];
    }

    return bias_sum / (high - low);
  }

  __host__ __device__ constexpr double distance(double e, int i) const noexcept
  {
    auto const diff = e - raw_estimate_data(this->precision_)[i];
    return diff * diff;
  }

  __host__ __device__ constexpr int interpolation_anchor_index(double e) const noexcept
  {
    auto estimates      = raw_estimate_data(this->precision_);
    int const n         = raw_estimate_data_size(this->precision_);
    int left            = 0;
    int right           = static_cast<int>(n) - 1;
    int mid             = -1;
    int candidate_index = 0;  // Index of the closest element found

    while (left <= right) {
      mid = left + (right - left) / 2;

      if (estimates[mid] < e) {
        left = mid + 1;
      } else if (estimates[mid] > e) {
        right = mid - 1;
      } else {
        // Exact match found, no need to look further
        return mid;
      }
    }

    // At this point, 'left' is the insertion point. We need to compare the elements at 'left' and
    // 'left - 1' to find the closest one, taking care of boundary conditions.

    // Distance from 'e' to the element at 'left', if within bounds
    double const dist_lhs = left < static_cast<int>(n) ? cuda::std::abs(estimates[left] - e)
                                                       : cuda::std::numeric_limits<double>::max();
    // Distance from 'e' to the element at 'left - 1', if within bounds
    double const dist_rhs = left - 1 >= 0 ? cuda::std::abs(estimates[left - 1] - e)
                                          : cuda::std::numeric_limits<double>::max();

    candidate_index = (dist_lhs < dist_rhs) ? left : left - 1;

    return candidate_index;
  }

  static constexpr auto k = 6;  ///< Number of interpolation points to consider
  int precision_;
  int m_;
};
}  // namespace cuco::hyperloglog_ns::detail
