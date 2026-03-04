/*
 * Copyright (c) 2024-2025, NVIDIA CORPORATION.
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

#include <cuco/detail/hyperloglog/hyperloglog_impl.cuh>
#include <cuco/hash_functions.cuh>
#include <cuco/types.cuh>
#include <cuco/utility/cuda_thread_scope.cuh>

#include <cuda/std/cstddef>
#include <cuda/stream_ref>

#include <cooperative_groups.h>

namespace cuco {
/**
 * @brief A GPU-accelerated utility for approximating the number of distinct items in a multiset.
 *
 * @note This implementation is based on the HyperLogLog++ algorithm:
 * https://static.googleusercontent.com/media/research.google.com/de//pubs/archive/40671.pdf.
 *
 * @tparam T Type of items to count
 * @tparam Scope The scope in which operations will be performed by individual threads
 * @tparam Hash Hash function used to hash items
 */
template <class T,
          cuda::thread_scope Scope = cuda::thread_scope_device,
          class Hash               = cuco::xxhash_64<T>>
class hyperloglog_ref {
  using impl_type = detail::hyperloglog_impl<T, Scope, Hash>;

 public:
  static constexpr auto thread_scope = impl_type::thread_scope;  ///< CUDA thread scope

  using value_type    = typename impl_type::value_type;     ///< Type of items to count
  using hasher        = typename impl_type::hasher;         ///< Type of hash function
  using register_type = typename impl_type::register_type;  ///< HLL register type

  template <cuda::thread_scope NewScope>
  using with_scope = hyperloglog_ref<T, NewScope, Hash>;  ///< Ref type with different
                                                          ///< thread scope

  /**
   * @brief Constructs a non-owning `hyperloglog_ref` object.
   *
   * @throw If sketch size < 0.0625KB or 64B or standard deviation > 0.2765. Throws if called from
   * host; UB if called from device.
   * @throw If sketch storage has insufficient alignment. Throws if called from host; UB if called
   * from device.
   *
   * @param sketch_span Reference to sketch storage
   * @param hash The hash function used to hash items
   */
  __host__ __device__ constexpr hyperloglog_ref(cuda::std::span<cuda::std::byte> sketch_span,
                                                Hash const& hash = {});

  /**
   * @brief Resets the estimator, i.e., clears the current count estimate.
   *
   * @tparam CG CUDA Cooperative Group type
   *
   * @param group CUDA Cooperative group this operation is executed in
   */
  template <class CG>
  __device__ constexpr void clear(CG group) noexcept;

  /**
   * @brief Asynchronously resets the estimator, i.e., clears the current count estimate.
   *
   * @param stream CUDA stream this operation is executed in
   */
  __host__ constexpr void clear_async(cuda::stream_ref stream = cuda::stream_ref{
                                        cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief Resets the estimator, i.e., clears the current count estimate.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `clear_async`.
   *
   * @param stream CUDA stream this operation is executed in
   */
  __host__ constexpr void clear(cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Adds an item to the estimator.
   *
   * @param item The item to be counted
   */
  __device__ constexpr void add(T const& item) noexcept;

  /**
   * @brief Asynchronously adds to be counted items to the estimator.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * T></tt> is `true`
   *
   * @param first Beginning of the sequence of items
   * @param last End of the sequence of items
   * @param stream CUDA stream this operation is executed in
   */
  template <class InputIt>
  __host__ constexpr void add_async(
    InputIt first, InputIt last, cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Adds to be counted items to the estimator.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `add_async`.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * T></tt> is `true`
   *
   * @param first Beginning of the sequence of items
   * @param last End of the sequence of items
   * @param stream CUDA stream this operation is executed in
   */
  template <class InputIt>
  __host__ constexpr void add(InputIt first,
                              InputIt last,
                              cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Merges the result of `other` estimator reference into `*this` estimator reference.
   *
   * @throw If this->sketch_bytes() != other.sketch_bytes() then behavior is undefined
   *
   * @tparam CG CUDA Cooperative Group type
   * @tparam OtherScope Thread scope of `other` estimator
   *
   * @param group CUDA Cooperative group this operation is executed in
   * @param other Other estimator reference to be merged into `*this`
   */
  template <class CG, cuda::thread_scope OtherScope>
  __device__ constexpr void merge(CG group, hyperloglog_ref<T, OtherScope, Hash> const& other);

  /**
   * @brief Asynchronously merges the result of `other` estimator reference into `*this` estimator.
   *
   * @throw If this->sketch_bytes() != other.sketch_bytes()
   *
   * @tparam OtherScope Thread scope of `other` estimator
   *
   * @param other Other estimator reference to be merged into `*this`
   * @param stream CUDA stream this operation is executed in
   */
  template <cuda::thread_scope OtherScope>
  __host__ constexpr void merge_async(hyperloglog_ref<T, OtherScope, Hash> const& other,
                                      cuda::stream_ref stream = cuda::stream_ref{
                                        cudaStream_t{nullptr}});

  /**
   * @brief Merges the result of `other` estimator reference into `*this` estimator.
   *
   * @throw If this->sketch_bytes() != other.sketch_bytes()
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `merge_async`.
   *
   * @tparam OtherScope Thread scope of `other` estimator
   *
   * @param other Other estimator reference to be merged into `*this`
   * @param stream CUDA stream this operation is executed in
   */
  template <cuda::thread_scope OtherScope>
  __host__ constexpr void merge(hyperloglog_ref<T, OtherScope, Hash> const& other,
                                cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Compute the estimated distinct items count.
   *
   * @param group CUDA thread block group this operation is executed in
   *
   * @return Approximate distinct items count
   */
  [[nodiscard]] __device__ std::size_t estimate(
    cooperative_groups::thread_block const& group) const noexcept;

  /**
   * @brief Compute the estimated distinct items count.
   *
   * @note This function synchronizes the given stream.
   *
   * @param stream CUDA stream this operation is executed in
   *
   * @return Approximate distinct items count
   */
  [[nodiscard]] __host__ constexpr std::size_t estimate(cuda::stream_ref stream = cuda::stream_ref{
                                                          cudaStream_t{nullptr}}) const;

  /**
   * @brief Gets the hash function.
   *
   * @return The hash function
   */
  [[nodiscard]] __host__ __device__ constexpr auto hash_function() const noexcept;

  /**
   * @brief Gets the span of the sketch.
   *
   * @return The cuda::std::span of the sketch
   */
  [[nodiscard]] __host__ __device__ constexpr cuda::std::span<cuda::std::byte> sketch()
    const noexcept;

  /**
   * @brief Gets the number of bytes required for the sketch storage.
   *
   * @return The number of bytes required for the sketch
   */
  [[nodiscard]] __host__ __device__ constexpr std::size_t sketch_bytes() const noexcept;

  /**
   * @brief Gets the number of bytes required for the sketch storage.
   *
   * @param sketch_size_kb Upper bound sketch size in KB
   *
   * @return The number of bytes required for the sketch
   */
  [[nodiscard]] __host__ __device__ static constexpr std::size_t sketch_bytes(
    cuco::sketch_size_kb sketch_size_kb) noexcept;

  /**
   * @brief Gets the number of bytes required for the sketch storage.
   *
   * @param standard_deviation Upper bound standard deviation for approximation error
   *
   * @return The number of bytes required for the sketch
   */
  [[nodiscard]] __host__ __device__ static constexpr std::size_t sketch_bytes(
    cuco::standard_deviation standard_deviation) noexcept;

  /**
   * @brief Gets the alignment required for the sketch storage.
   *
   * @return The required alignment
   */
  [[nodiscard]] __host__ __device__ static constexpr std::size_t sketch_alignment() noexcept;

 private:
  impl_type impl_;  ///< Implementation object

  template <class T_, cuda::thread_scope Scope_, class Hash_>
  friend class hyperloglog_ref;
};
}  // namespace cuco

#include <cuco/detail/hyperloglog/hyperloglog_ref.inl>
