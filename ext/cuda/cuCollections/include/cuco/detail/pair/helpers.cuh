/*
 * Copyright (c) 2021-2025, NVIDIA CORPORATION.
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

#include <cuda/functional>
#include <cuda/std/bit>
#include <cuda/std/type_traits>

#include <cstdint>

namespace cuco::detail {

/**
 * @brief Gives value to use as alignment for a pair type that is at least the
 * size of the sum of the size of the first type and second type, or 16,
 * whichever is smaller.
 */
template <typename First, typename Second>
__host__ __device__ constexpr std::size_t pair_alignment()
{
  constexpr std::size_t alignment = cuda::std::bit_ceil(sizeof(First) + sizeof(Second));
  return cuda::std::min(std::size_t{16}, alignment);
}

/**
 * @brief Denotes the equivalent packed type based on the size of the object.
 *
 * @tparam N The size of the object
 */
template <std::size_t N>
struct packed {
  using type = void;  ///< `void` type by default
};

/**
 * @brief Denotes the packed type when the size of the object is 8.
 */
template <>
struct packed<sizeof(uint64_t)> {
  using type = uint64_t;  ///< Packed type as `uint64_t` if the size of the object is 8
};

/**
 * @brief Denotes the packed type when the size of the object is 4.
 */
template <>
struct packed<sizeof(uint32_t)> {
  using type = uint32_t;  ///< Packed type as `uint32_t` if the size of the object is 4
};

template <typename Pair>
using packed_t = typename packed<sizeof(Pair)>::type;

/**
 * @brief Indicates if a pair type can be packed.
 *
 * When the size of the key,value pair being inserted into the hash table is
 * equal in size to a type where atomicCAS is natively supported, it is more
 * efficient to "pack" the pair and insert it with a single atomicCAS.
 *
 * Pair types whose key and value have the same object representation may be
 * packed. Also, the `Pair` must not contain any padding bits otherwise
 * accessing the packed value would be undefined.
 *
 * @tparam Pair The pair type that will be packed
 *
 * @return true If the pair type can be packed
 * @return false  If the pair type cannot be packed
 */
template <typename Pair>
__host__ __device__ constexpr bool is_packable()
{
  return not cuda::std::is_void<packed_t<Pair>>::value and
         cuda::std::has_unique_object_representations_v<Pair>;
}

/**
 * @brief Allows viewing a pair in a packed representation.
 *
 * Used as an optimization for inserting when a pair can be inserted with a
 * single atomicCAS
 */
template <typename Pair>
union pair_converter {
  using packed_type = packed_t<Pair>;  ///< The packed pair type
  packed_type packed;                  ///< The pair in the packed representation
  Pair pair;                           ///< The pair in the pair representation

  /**
   * @brief Constructs a pair converter by copying from `p`
   *
   * @tparam T Type that is convertible to `Pair`
   *
   * @param p The pair to copy from
   */
  template <typename T>
  __device__ pair_converter(T&& p) : pair{p}
  {
  }

  /**
   * @brief Constructs a pair converter by copying from `p`
   *
   * @param p The packed data to copy from
   */
  __device__ pair_converter(packed_type p) : packed{p} {}
};

}  // namespace cuco::detail
