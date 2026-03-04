/*
 * Copyright (c) 2022-2025, NVIDIA CORPORATION.
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

#include <cstddef>
#include <cstdint>

namespace cuco {
static constexpr std::size_t dynamic_extent = static_cast<std::size_t>(-1);

/**
 * @brief Static extent class.
 *
 * @tparam SizeType Size type
 * @tparam N Extent
 */
template <typename SizeType, std::size_t N = dynamic_extent>
struct extent {
  using value_type = SizeType;  ///< Extent value type

  constexpr extent() = default;

  /// Constructs from `SizeType`
  __host__ __device__ constexpr extent(SizeType) noexcept {}

  /**
   * @brief Conversion to value_type.
   *
   * @return Extent size
   */
  __host__ __device__ constexpr operator value_type() const noexcept { return N; }
};

/**
 * @brief Dynamic extent class.
 *
 * @tparam SizeType Size type
 */
template <typename SizeType>
struct extent<SizeType, dynamic_extent> {
  using value_type = SizeType;  ///< Extent value type

  /**
   * @brief Constructs extent from a given `size`.
   *
   * @param size The extent size
   */
  __host__ __device__ constexpr extent(SizeType size) noexcept : value_{size} {}

  /**
   * @brief Conversion to value_type.
   *
   * @return Extent size
   */
  __host__ __device__ constexpr operator value_type() const noexcept { return value_; }

 private:
  value_type value_;  ///< Extent value
};

/**
 * @brief Computes valid bucket extent/capacity based on given parameters.
 *
 * @note The actual capacity of a hash table should be exclusively determined by the return
 * value of this utility since the output depends on the requested low-bound size, the probing
 * scheme, and the storage. This utility is used internally during container constructions while for
 * container ref constructions, it would be users' responsibility to use this function to determine
 * the capacity ctor argument for the container.
 *
 * @tparam ProbingScheme Type of probing scheme
 * @tparam Storage Type of storage
 * @tparam SizeType Size type
 * @tparam N Extent size (can be dynamic_extent)
 *
 * @param ext The input extent
 *
 * @throw If the input extent is invalid
 *
 * @return Resulting valid extent
 */
template <typename ProbingScheme, typename Storage, typename SizeType, std::size_t N>
[[nodiscard]] auto constexpr make_valid_extent(cuco::extent<SizeType, N> ext);

template <template <typename> class ProbingScheme,
          typename Storage,
          typename SizeType,
          std::size_t N>
[[nodiscard]] auto constexpr make_valid_extent(cuco::extent<SizeType, N> ext);

template <template <typename, typename> class ProbingScheme,
          typename Storage,
          typename SizeType,
          std::size_t N>
[[nodiscard]] auto constexpr make_valid_extent(cuco::extent<SizeType, N> ext);

/**
 * @brief Computes valid bucket extent/capacity based on given parameters.
 *
 * @note The actual capacity of a hash table should be exclusively determined by the return
 * value of this utility since the output depends on the requested low-bound size, the probing
 * scheme, and the storage. This utility is used internally during container constructions while for
 * container ref constructions, it would be users' responsibility to use this function to determine
 * the capacity ctor argument for the container.
 *
 * @tparam ProbingScheme Type of probing scheme
 * @tparam Storage Type of storage
 * @tparam SizeType Size type
 *
 * @param ext The input size value
 *
 * @throw If the input size is invalid
 *
 * @return Resulting valid bucket extent
 */
template <typename ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType ext);

template <template <typename> class ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType ext);

template <template <typename, typename> class ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType ext);

/**
 * @brief Computes valid bucket extent based on given parameters and desired load factor.
 *
 * @tparam ProbingScheme Type of probing scheme
 * @tparam Storage Type of storage
 * @tparam SizeType Size type
 *
 * @param ext The input extent
 * @param desired_load_factor The desired load factor (e.g., 0.5 for 50%)
 *
 * @throw If the desired load factor is invalid
 *
 * @return Resulting valid extent
 */
template <typename ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(cuco::extent<SizeType> ext,
                                               double desired_load_factor);

/**
 * @brief Computes valid bucket extent based on given size and desired load factor.
 *
 * @tparam ProbingScheme Type of probing scheme
 * @tparam Storage Type of storage
 * @tparam SizeType Size type
 *
 * @param ext The input extent
 * @param desired_load_factor The desired load factor (e.g., 0.5 for 50%)
 *
 * @throw If the desired load factor is invalid
 *
 * @return Resulting valid extent
 */
template <typename ProbingScheme, typename Storage, typename SizeType>
[[nodiscard]] auto constexpr make_valid_extent(SizeType ext, double desired_load_factor);
}  // namespace cuco

#include <cuco/detail/extent/extent.inl>
