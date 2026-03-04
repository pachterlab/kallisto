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

#include <cuco/detail/utils.cuh>
#include <cuco/pair.cuh>

#include <algorithm>
#include <cstdint>

namespace cuco {
namespace detail {

/**
 * @brief Probing iterator class.
 *
 * @tparam Extent Type of Extent
 */
template <typename Extent>
class probing_iterator {
 public:
  using extent_type = Extent;                            ///< Extent type
  using size_type   = typename extent_type::value_type;  ///< Size type

  /**
   * @brief Constructs an probing iterator
   *
   * @param start Iteration starting point
   * @param step_size Double hashing step size
   * @param upper_bound Upper bound of the iteration
   */
  __host__ __device__ constexpr probing_iterator(size_type start,
                                                 size_type step_size,
                                                 extent_type upper_bound) noexcept
    : curr_index_{start}, step_size_{step_size}, upper_bound_{upper_bound}
  {
    // TODO: revise this API when introducing quadratic probing into cuco
  }

  /**
   * @brief Dereference operator
   *
   * @return Current slot index
   */
  __host__ __device__ constexpr auto operator*() const noexcept { return curr_index_; }

  /**
   * @brief Prefix increment operator
   *
   * @return Current iterator
   */
  __host__ __device__ constexpr auto operator++() noexcept
  {
    // TODO: step_size_ can be a build time constant (e.g. linear probing)
    //  Worth passing another extent type?
    curr_index_ = (curr_index_ + step_size_) % upper_bound_;
    return *this;
  }

  /**
   * @brief Postfix increment operator
   *
   * @return Old iterator before increment
   */
  __host__ __device__ constexpr auto operator++(int32_t) noexcept
  {
    auto temp = *this;
    ++(*this);
    return temp;
  }

 private:
  size_type curr_index_;
  size_type step_size_;
  extent_type upper_bound_;
};
}  // namespace detail

template <int32_t CGSize, typename Hash>
__host__ __device__ constexpr linear_probing<CGSize, Hash>::linear_probing(Hash const& hash)
  : hash_{hash}
{
}

template <int32_t CGSize, typename Hash>
template <typename NewHash>
__host__ __device__ constexpr auto linear_probing<CGSize, Hash>::rebind_hash_function(
  NewHash const& hash) const noexcept
{
  return linear_probing<cg_size, NewHash>{hash};
}

template <int32_t CGSize, typename Hash>
template <int32_t BucketSize, typename ProbeKey, typename Extent>
__host__ __device__ constexpr auto linear_probing<CGSize, Hash>::make_iterator(
  ProbeKey probe_key, Extent upper_bound) const noexcept
{
  using size_type      = typename Extent::value_type;
  size_type const init = cuco::detail::sanitize_hash<size_type>(hash_(probe_key)) %
                         (upper_bound / BucketSize) * BucketSize;
  return detail::probing_iterator<Extent>{init, static_cast<size_type>(BucketSize), upper_bound};
}

template <int32_t CGSize, typename Hash>
template <int32_t BucketSize, typename ProbeKey, typename Extent, typename ParentCG>
__host__ __device__ constexpr auto linear_probing<CGSize, Hash>::make_iterator(
  cooperative_groups::thread_block_tile<cg_size, ParentCG> g,
  ProbeKey probe_key,
  Extent upper_bound) const noexcept
{
  using size_type            = typename Extent::value_type;
  size_type constexpr stride = cg_size * BucketSize;
  size_type const init =
    cuco::detail::sanitize_hash<size_type>(hash_(probe_key)) % (upper_bound / stride) * stride +
    g.thread_rank() * BucketSize;
  return detail::probing_iterator<Extent>{init, stride, upper_bound};
}

template <int32_t CGSize, typename Hash>
__host__ __device__ constexpr linear_probing<CGSize, Hash>::hasher
linear_probing<CGSize, Hash>::hash_function() const noexcept
{
  return hash_;
}

template <int32_t CGSize, typename Hash1, typename Hash2>
__host__ __device__ constexpr double_hashing<CGSize, Hash1, Hash2>::double_hashing(
  Hash1 const& hash1, Hash2 const& hash2)
  : hash1_{hash1}, hash2_{hash2}
{
}

template <int32_t CGSize, typename Hash1, typename Hash2>
__host__ __device__ constexpr double_hashing<CGSize, Hash1, Hash2>::double_hashing(
  cuda::std::tuple<Hash1, Hash2> const& hash)
  : hash1_{hash.first}, hash2_{hash.second}
{
}

template <int32_t CGSize, typename Hash1, typename Hash2>
template <typename NewHash, typename Enable>
__host__ __device__ constexpr auto double_hashing<CGSize, Hash1, Hash2>::rebind_hash_function(
  NewHash const& hash) const
{
  static_assert(cuco::is_tuple_like<NewHash>::value,
                "The given hasher must be a tuple-like object");

  auto const [hash1, hash2] = cuda::std::tuple{hash};
  using hash1_type          = cuda::std::decay_t<decltype(hash1)>;
  using hash2_type          = cuda::std::decay_t<decltype(hash2)>;
  return double_hashing<cg_size, hash1_type, hash2_type>{hash1, hash2};
}

template <int32_t CGSize, typename Hash1, typename Hash2>
template <int32_t BucketSize, typename ProbeKey, typename Extent>
__host__ __device__ constexpr auto double_hashing<CGSize, Hash1, Hash2>::make_iterator(
  ProbeKey probe_key, Extent upper_bound) const noexcept
{
  using size_type = typename Extent::value_type;
  return detail::probing_iterator<Extent>{
    static_cast<size_type>(cuco::detail::sanitize_hash<size_type>(hash1_(probe_key)) %
                           (upper_bound / BucketSize) * BucketSize),
    static_cast<size_type>(
      (cuco::detail::sanitize_hash<size_type>(hash2_(probe_key)) % (upper_bound / BucketSize - 1) +
       1) *
      BucketSize),  // step size in range [1, prime - 1]
    upper_bound};
}

template <int32_t CGSize, typename Hash1, typename Hash2>
template <int32_t BucketSize, typename ProbeKey, typename Extent, typename ParentCG>
__host__ __device__ constexpr auto double_hashing<CGSize, Hash1, Hash2>::make_iterator(
  cooperative_groups::thread_block_tile<cg_size, ParentCG> g,
  ProbeKey probe_key,
  Extent upper_bound) const noexcept
{
  int32_t const stride = cg_size * BucketSize;
  using size_type      = typename Extent::value_type;

  return detail::probing_iterator<Extent>{
    static_cast<size_type>(cuco::detail::sanitize_hash<size_type>(hash1_(probe_key)) %
                             (upper_bound / stride) * stride +
                           g.thread_rank() * BucketSize),
    static_cast<size_type>(
      (cuco::detail::sanitize_hash<size_type>(hash2_(probe_key)) % (upper_bound / stride - 1) + 1) *
      stride),
    upper_bound};  // TODO use fast_int operator
}

template <int32_t CGSize, typename Hash1, typename Hash2>
__host__ __device__ constexpr double_hashing<CGSize, Hash1, Hash2>::hasher
double_hashing<CGSize, Hash1, Hash2>::hash_function() const noexcept
{
  return {hash1_, hash2_};
}

}  // namespace cuco
