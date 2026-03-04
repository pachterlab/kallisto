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

#include <cuco/detail/bloom_filter/bloom_filter_impl.cuh>
#include <cuco/utility/cuda_thread_scope.cuh>

#include <cuda/atomic>
#include <cuda/stream_ref>

namespace cuco {

/**
 * @brief Non-owning "ref" type of `bloom_filter`.
 *
 * @note Ref types are trivially-copyable and are intended to be passed by value.
 *
 * @tparam Key Key type
 * @tparam Extent Size type that is used to determine the number of blocks in the filter
 * @tparam Scope The scope in which operations will be performed by individual threads
 * @tparam Policy Type that defines how to generate and store key fingerprints (see
 * `cuco/bloom_filter_policies.cuh`)
 */
template <class Key, class Extent, cuda::thread_scope Scope, class Policy>
class bloom_filter_ref {
  using impl_type = detail::bloom_filter_impl<Key, Extent, Scope, Policy>;  ///< Implementation type

 public:
  static constexpr auto thread_scope = impl_type::thread_scope;  ///< CUDA thread scope
  static constexpr auto words_per_block =
    impl_type::words_per_block;  ///< Number of words/segments in each filter block

  using key_type    = typename impl_type::key_type;      ///< Key Type
  using extent_type = typename impl_type::extent_type;   ///< Extent type
  using size_type   = typename extent_type::value_type;  ///< Underlying type of the extent type
  using word_type =
    typename impl_type::word_type;  ///< Underlying word/segment type of a filter block
  using filter_block_type =
    typename impl_type::filter_block_type;  ///< Opaque type of a filter block

  /**
   * @brief Constructs the ref object from existing storage.
   *
   * @note The storage span starting at `data` must have an extent of at least `num_blocks`
   * elements of type `filter_block_type`.
   *
   * @param data Pointer to the storage span of the filter
   * @param num_blocks Number of sub-filters or blocks
   * @param scope The scope in which operations will be performed
   * @param policy Fingerprint generation policy (see `cuco/bloom_filter_policies.cuh`)
   */
  __host__ __device__ explicit constexpr bloom_filter_ref(filter_block_type* data,
                                                          Extent num_blocks,
                                                          cuda_thread_scope<Scope> scope,
                                                          Policy const& policy);

  /**
   * @brief Constructs the ref object from existing storage.
   *
   * @note This overload is deprecated and will be removed in the near future.
   *
   * @param data Pointer to the storage span of the filter
   * @param num_blocks Number of sub-filters or blocks
   * @param scope The scope in which operations will be performed
   * @param policy Fingerprint generation policy (see `cuco/bloom_filter_policies.cuh`)
   */
  __host__ __device__ explicit constexpr bloom_filter_ref(word_type* data,
                                                          Extent num_blocks,
                                                          cuda_thread_scope<Scope> scope,
                                                          Policy const& policy);

  /**
   * @brief Device function that cooperatively erases all information from the filter.
   *
   * @tparam CG Cooperative Group type
   *
   * @param group The Cooperative Group this operation is executed with
   */
  template <class CG>
  __device__ constexpr void clear(CG group);

  /**
   * @brief Erases all information from the filter.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `clear_async`.
   *
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  __host__ constexpr void clear(cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously erases all information from the filter.
   *
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  __host__ constexpr void clear_async(cuda::stream_ref stream = cuda::stream_ref{
                                        cudaStream_t{nullptr}});

  /**
   * @brief Device function that adds a key to the filter.
   *
   * @tparam ProbeKey Input type that is implicitly convertible to `key_type`
   *
   * @param key The key to be added
   */
  template <class ProbeKey>
  __device__ void add(ProbeKey const& key);

  /**
   * @brief Device function that cooperatively adds a key to the filter.
   *
   * @note Best performance is achieved if the size of the CG is equal to `words_per_block`.
   *
   * @tparam CG Cooperative Group type
   * @tparam ProbeKey Input key type
   *
   * @param group The Cooperative Group this operation is executed with
   * @param key The key to be added
   */
  template <class CG, class ProbeKey>
  __device__ void add(CG group, ProbeKey const& key);

  /**
   * @brief Device function that adds all keys in the range `[first, last)` to the filter.
   *
   * @note Best performance is achieved if the size of the CG is larger than or equal to
   * `words_per_block`.
   *
   * @tparam CG Cooperative Group type
   * @tparam InputIt Device-accessible random access input key iterator
   *
   * @param group The Cooperative Group this operation is executed with
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   */
  template <class CG, class InputIt>
  __device__ void add(CG group, InputIt first, InputIt last);

  /**
   * @brief Adds all keys in the range `[first, last)` to the filter.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `add_async`.
   *
   * @tparam InputIt Device-accessible random access input key iterator
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt>
  __host__ constexpr void add(InputIt first,
                              InputIt last,
                              cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously adds all keys in the range `[first, last)` to the filter.
   *
   * @tparam InputIt Device-accessible random access input key iterator
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt>
  __host__ constexpr void add_async(
    InputIt first, InputIt last, cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Adds keys in the range `[first, last)` if `pred` of the corresponding `stencil` returns
   * `true`.
   *
   * @note The key `*(first + i)` is added if `pred( *(stencil + i) )` returns `true`.
   * @note This function synchronizes the given stream and returns the number of successful
   * insertions. For asynchronous execution use `add_if_async`.
   *
   * @tparam InputIt Device-accessible random access input key iterator
   * @tparam StencilIt Device-accessible random-access iterator whose `value_type` is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt, class StencilIt, class Predicate>
  __host__ constexpr void add_if(InputIt first,
                                 InputIt last,
                                 StencilIt stencil,
                                 Predicate pred,
                                 cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously adds keys in the range `[first, last)` if `pred` of the corresponding
   * `stencil` returns `true`.
   *
   * @note The key `*(first + i)` is added if `pred( *(stencil + i) )` returns `true`.
   *
   * @tparam InputIt Device-accessible random access input key iterator
   * @tparam StencilIt Device-accessible random-access iterator whose `value_type` is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt, class StencilIt, class Predicate>
  __host__ constexpr void add_if_async(InputIt first,
                                       InputIt last,
                                       StencilIt stencil,
                                       Predicate pred,
                                       cuda::stream_ref stream = cuda::stream_ref{
                                         cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief Device function that tests if a key's fingerprint is present in the filter.
   *
   * @tparam ProbeKey Probe key type
   *
   * @param key The key to be tested
   *
   * @return `true` iff the key's fingerprint was present in the filter
   */
  template <class ProbeKey>
  [[nodiscard]] __device__ bool contains(ProbeKey const& key) const;

  /**
   * @brief Device function that tests if a key's fingerprint is present in the filter.
   *
   * @note Best performance is achieved if the size of the CG is equal to `(words_per_block *
   * sizeof(word_type)) / 32`.
   *
   * @tparam CG Cooperative Group type
   * @tparam ProbeKey Input type that is implicitly convertible to `key_type`
   *
   * @param group The Cooperative Group this operation is executed with
   * @param key The key to be tested
   *
   * @return `true` iff the key's fingerprint was present in the filter
   */
  template <class CG, class ProbeKey>
  [[nodiscard]] __device__ bool contains(CG group, ProbeKey const& key) const;

  // TODO
  // template <class CG, class InputIt, class OutputIt>
  // __device__ void contains(CG group, InputIt first, InputIt last, OutputIt output_begin)
  // const;

  /**
   * @brief Tests all keys in the range `[first, last)` if their fingerprints are present in the
   * filter.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `contains_async`.
   *
   * @tparam InputIt Device-accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * bloom_filter<K>::key_type></tt> is `true`
   * @tparam OutputIt Device-accessible output iterator assignable from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt, class OutputIt>
  __host__ constexpr void contains(InputIt first,
                                   InputIt last,
                                   OutputIt output_begin,
                                   cuda::stream_ref stream = cuda::stream_ref{
                                     cudaStream_t{nullptr}}) const;

  /**
   * @brief Asynchronously tests all keys in the range `[first, last)` if their fingerprints are
   * present in the filter.
   *
   * @tparam InputIt Device-accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * bloom_filter<K>::key_type></tt> is `true`
   * @tparam OutputIt Device-accessible output iterator assignable from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt, class OutputIt>
  __host__ constexpr void contains_async(InputIt first,
                                         InputIt last,
                                         OutputIt output_begin,
                                         cuda::stream_ref stream = cuda::stream_ref{
                                           cudaStream_t{nullptr}}) const noexcept;

  /**
   * @brief Tests all keys in the range `[first, last)` if their fingerprints are present in the
   * filter if `pred` of the corresponding `stencil` returns `true`.
   *
   * @note The key `*(first + i)` is queried if `pred( *(stencil + i) )` returns `true`.
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `contains_if_async`.
   *
   * @tparam InputIt Device-accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * bloom_filter<K>::key_type></tt> is `true`
   * @tparam StencilIt Device-accessible random-access iterator whose `value_type` is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam OutputIt Device-accessible output iterator assignable from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt, class StencilIt, class Predicate, class OutputIt>
  __host__ constexpr void contains_if(InputIt first,
                                      InputIt last,
                                      StencilIt stencil,
                                      Predicate pred,
                                      OutputIt output_begin,
                                      cuda::stream_ref stream = cuda::stream_ref{
                                        cudaStream_t{nullptr}}) const;

  /**
   * @brief Asynchronously tests all keys in the range `[first, last)` if their fingerprints are
   * present in the filter if `pred` of the corresponding `stencil` returns `true`.
   *
   * @note The key `*(first + i)` is queried if `pred( *(stencil + i) )` returns `true`.
   *
   * @tparam InputIt Device-accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * bloom_filter<K>::key_type></tt> is `true`
   * @tparam StencilIt Device-accessible random-access iterator whose `value_type` is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam OutputIt Device-accessible output iterator assignable from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream CUDA stream used for device memory operations and kernel launches
   */
  template <class InputIt, class StencilIt, class Predicate, class OutputIt>
  __host__ constexpr void contains_if_async(InputIt first,
                                            InputIt last,
                                            StencilIt stencil,
                                            Predicate pred,
                                            OutputIt output_begin,
                                            cuda::stream_ref stream = cuda::stream_ref{
                                              cudaStream_t{nullptr}}) const noexcept;

  /**
   * @brief Gets a pointer to the underlying filter storage.
   *
   * @return Pointer to the underlying filter storage
   */
  [[nodiscard]] __host__ __device__ constexpr word_type* data() noexcept;

  /**
   * @brief Gets a pointer to the underlying filter storage.
   *
   * @return Pointer to the underlying filter storage
   */
  [[nodiscard]] __host__ __device__ constexpr word_type const* data() const noexcept;

  /**
   * @brief Gets the number of sub-filter blocks.
   *
   * @return Number of sub-filter blocks
   */
  [[nodiscard]] __host__ __device__ constexpr extent_type block_extent() const noexcept;

 private:
  impl_type impl_;  ///< Object containing the Blocked Bloom Filter implementation
};
}  // namespace cuco

#include <cuco/detail/bloom_filter/bloom_filter_ref.inl>