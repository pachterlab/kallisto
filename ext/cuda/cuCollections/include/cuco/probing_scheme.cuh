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

#include <cuco/detail/probing_scheme/probing_scheme_base.cuh>
#include <cuco/pair.cuh>

#include <cuda/std/tuple>
#include <cuda/std/type_traits>

#include <cooperative_groups.h>

namespace cuco {
/**
 * @brief Public linear probing scheme class.
 *
 * @note Linear probing is efficient when few collisions are present, e.g., low occupancy or low
 * multiplicity.
 *
 * @note `Hash` should be callable object type.
 *
 * @tparam CGSize Size of CUDA Cooperative Groups
 * @tparam Hash Unary callable type
 */
template <int32_t CGSize, typename Hash>
class linear_probing : private detail::probing_scheme_base<CGSize> {
  using probing_scheme_base_type =
    detail::probing_scheme_base<CGSize>;  ///< The base probe scheme type

 public:
  using probing_scheme_base_type::cg_size;
  using hasher = Hash;  ///< Hash function type

  /**
   *@brief Constructs linear probing scheme with the hasher callable.
   *
   * @param hash Hasher
   */
  __host__ __device__ constexpr linear_probing(Hash const& hash = {});

  /**
   *@brief Makes a copy of the current probing method with the given hasher
   *
   * @tparam NewHash New hasher type
   *
   * @param hash Hasher
   *
   * @return Copy of the current probing method
   */
  template <typename NewHash>
  [[nodiscard]] __host__ __device__ constexpr auto rebind_hash_function(
    NewHash const& hash) const noexcept;

  /**
   * @brief Returns a probing iterator
   *
   * @tparam BucketSize Size of the bucket
   * @tparam ProbeKey Type of probing key
   * @tparam Extent Type of extent
   *
   * @param probe_key The probing key
   * @param upper_bound Upper bound of the iteration
   * @return An iterator whose value_type is convertible to slot index type
   */
  template <int32_t BucketSize, typename ProbeKey, typename Extent>
  __host__ __device__ constexpr auto make_iterator(ProbeKey probe_key,
                                                   Extent upper_bound) const noexcept;

  /**
   * @brief Returns a CG-based probing iterator
   *
   * @tparam BucketSize Size of the bucket
   * @tparam ProbeKey Type of probing key
   * @tparam Extent Type of extent
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g the Cooperative Group to generate probing iterator
   * @param probe_key The probing key
   * @param upper_bound Upper bound of the iteration
   * @return An iterator whose value_type is convertible to slot index type
   */
  template <int32_t BucketSize, typename ProbeKey, typename Extent, typename ParentCG>
  __host__ __device__ constexpr auto make_iterator(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> g,
    ProbeKey probe_key,
    Extent upper_bound) const noexcept;

  /**
   * @brief Gets the function used to hash keys
   *
   * @return The function used to hash keys
   */
  __host__ __device__ constexpr hasher hash_function() const noexcept;

 private:
  Hash hash_;
};

/**
 * @brief Public double hashing scheme class.
 *
 * @note Default probing scheme for cuco data structures. It shows superior performance over linear
 * probing especially when dealing with high multiplicty and/or high occupancy use cases.
 *
 * @note `Hash1` and `Hash2` should be callable object type.
 *
 * @note `Hash2` needs to be able to construct from an integer value to avoid secondary clustering.
 *
 * @tparam CGSize Size of CUDA Cooperative Groups
 * @tparam Hash1 Unary callable type
 * @tparam Hash2 Unary callable type
 */
template <int32_t CGSize, typename Hash1, typename Hash2 = Hash1>
class double_hashing : private detail::probing_scheme_base<CGSize> {
  using probing_scheme_base_type =
    detail::probing_scheme_base<CGSize>;  ///< The base probe scheme type

 public:
  using probing_scheme_base_type::cg_size;
  using hasher = cuda::std::tuple<Hash1, Hash2>;  ///< Hash function type

  /**
   *@brief Constructs double hashing probing scheme with the two hasher callables.
   *
   * @param hash1 First hasher
   * @param hash2 Second hasher
   */
  __host__ __device__ constexpr double_hashing(Hash1 const& hash1 = {}, Hash2 const& hash2 = {1});

  /**
   *@brief Constructs double hashing probing scheme with the hasher tuple
   *
   * @param hash Hasher tuple
   */
  __host__ __device__ constexpr double_hashing(cuda::std::tuple<Hash1, Hash2> const& hash);

  /**
   *@brief Makes a copy of the current probing method with the given hasher
   *
   * @tparam NewHash Tuple-like new hasher type
   *
   * @throw If `cuco::is_tuple_like_v<NewHash>` is `false`
   *
   * @param hash Hasher
   *
   * @return Copy of the current probing method
   */
  template <typename NewHash,
            typename Enable = cuda::std::enable_if_t<cuco::is_tuple_like<NewHash>::value>>
  [[nodiscard]] __host__ __device__ constexpr auto rebind_hash_function(NewHash const& hash) const;

  /**
   * @brief Returns a probing iterator
   *
   * @tparam BucketSize Size of the bucket
   * @tparam ProbeKey Type of probing key
   * @tparam Extent Type of extent
   *
   * @param probe_key The probing key
   * @param upper_bound Upper bound of the iteration
   * @return An iterator whose value_type is convertible to slot index type
   */
  template <int32_t BucketSize, typename ProbeKey, typename Extent>
  __host__ __device__ constexpr auto make_iterator(ProbeKey probe_key,
                                                   Extent upper_bound) const noexcept;

  /**
   * @brief Returns a CG-based probing iterator
   *
   * @tparam BucketSize Size of the bucket
   * @tparam ProbeKey Type of probing key
   * @tparam Extent Type of extent
   * @tparam ParentCG Type of parent Cooperative Group
   *
   * @param g the Cooperative Group to generate probing iterator
   * @param probe_key The probing key
   * @param upper_bound Upper bound of the iteration
   * @return An iterator whose value_type is convertible to slot index type
   */
  template <int32_t BucketSize, typename ProbeKey, typename Extent, typename ParentCG>
  __host__ __device__ constexpr auto make_iterator(
    cooperative_groups::thread_block_tile<cg_size, ParentCG> g,
    ProbeKey probe_key,
    Extent upper_bound) const noexcept;

  /**
   * @brief Gets the functions used to hash keys
   *
   * @return The functions used to hash keys
   */
  __host__ __device__ constexpr hasher hash_function() const noexcept;

 private:
  Hash1 hash1_;
  Hash2 hash2_;
};

/**
 * @brief Trait indicating whether the given probing scheme is of `double_hashing` type or not
 *
 * @tparam T Input probing scheme type
 */
template <typename T>
struct is_double_hashing : cuda::std::false_type {};

/**
 * @brief Trait indicating whether the given probing scheme is of `double_hashing` type or not
 *
 * @tparam CGSize Size of CUDA Cooperative Groups
 * @tparam Hash1 Unary callable type
 * @tparam Hash2 Unary callable type
 */
template <int32_t CGSize, typename Hash1, typename Hash2>
struct is_double_hashing<cuco::double_hashing<CGSize, Hash1, Hash2>> : cuda::std::true_type {};

}  // namespace cuco

#include <cuco/detail/probing_scheme/probing_scheme_impl.inl>
