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

#include <cuco/detail/open_addressing/open_addressing_impl.cuh>
#include <cuco/extent.cuh>
#include <cuco/hash_functions.cuh>
#include <cuco/probing_scheme.cuh>
#include <cuco/static_multiset_ref.cuh>
#include <cuco/storage.cuh>
#include <cuco/types.cuh>
#include <cuco/utility/allocator.hpp>
#include <cuco/utility/cuda_thread_scope.cuh>
#include <cuco/utility/traits.hpp>

#include <cuda/atomic>
#include <cuda/std/functional>
#include <cuda/stream_ref>

#include <cstddef>
#include <memory>

namespace cuco {
/**
 * @brief A GPU-accelerated, unordered, associative container of possibly non-unique objects
 *
 * The `static_multiset` supports two types of operations:
 * - Host-side "bulk" operations
 * - Device-side "singular" operations
 *
 * The host-side bulk operations include `insert`, `contains`, etc. These APIs should be used when
 * there are a large number of keys to modify or lookup. For example, given a range of keys
 * specified by device-accessible iterators, the bulk `insert` function will insert all keys into
 * the set.
 *
 * The singular device-side operations allow individual threads (or cooperative groups) to perform
 * independent modify or lookup operations from device code. These operations are accessed through
 * non-owning, trivially copyable reference types (or "ref"). User can combine any arbitrary
 * operators (see options in `include/cuco/operator.hpp`) when creating the ref. Concurrent modify
 * and lookup will be supported if both kinds of operators are specified during the ref
 * construction.
 *
 * @note Allows constant time concurrent modify or lookup operations from threads in device code.
 * @note cuCollections data structures always place the slot keys on the right-hand side when
 * invoking the key comparison predicate, i.e., `pred(query_key, slot_key)`. Order-sensitive
 * `KeyEqual` should be used with caution.
 * @note `ProbingScheme::cg_size` indicates how many threads are used to handle one independent
 * device operation. `cg_size == 1` uses the scalar (or non-CG) code paths.
 *
 * @throw If the size of the given key type is larger than 8 bytes
 * @throw If the given key type doesn't have unique object representations, i.e.,
 * `cuco::bitwise_comparable_v<Key> == false`
 * @throw If the probing scheme type is not inherited from `cuco::detail::probing_scheme_base`
 *
 * @tparam Key Type used for keys. Requires `cuco::is_bitwise_comparable_v<Key>`
 * @tparam Extent Data structure size type
 * @tparam Scope The scope in which operations will be performed by individual threads.
 * @tparam KeyEqual Binary callable type used to compare two keys for equality
 * @tparam ProbingScheme Probing scheme (see `include/cuco/probing_scheme.cuh` for choices)
 * @tparam Allocator Type of allocator used for device storage
 * @tparam Storage Slot bucket storage type
 */
template <class Key,
          class Extent             = cuco::extent<std::size_t>,
          cuda::thread_scope Scope = cuda::thread_scope_device,
          class KeyEqual           = cuda::std::equal_to<Key>,
          class ProbingScheme      = cuco::double_hashing<4,  // CG size
                                                          cuco::default_hash_function<Key>>,
          class Allocator          = cuco::cuda_allocator<Key>,
          class Storage            = cuco::storage<2>>
class static_multiset {
  using impl_type = detail::
    open_addressing_impl<Key, Key, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>;

 public:
  static constexpr auto cg_size      = impl_type::cg_size;       ///< CG size used for probing
  static constexpr auto bucket_size  = impl_type::bucket_size;   ///< Bucket size used for probing
  static constexpr auto thread_scope = impl_type::thread_scope;  ///< CUDA thread scope

  using key_type       = typename impl_type::key_type;        ///< Key type
  using value_type     = typename impl_type::value_type;      ///< Key type
  using extent_type    = typename impl_type::extent_type;     ///< Extent type
  using size_type      = typename impl_type::size_type;       ///< Size type
  using key_equal      = typename impl_type::key_equal;       ///< Key equality comparator type
  using allocator_type = typename impl_type::allocator_type;  ///< Allocator type
  /// Non-owning bucket storage ref type
  using storage_ref_type    = typename impl_type::storage_ref_type;
  using probing_scheme_type = typename impl_type::probing_scheme_type;  ///< Probing scheme type
  using hasher              = typename probing_scheme_type::hasher;     ///< Hash function type
  template <typename... Operators>
  using ref_type = cuco::static_multiset_ref<key_type,
                                             thread_scope,
                                             key_equal,
                                             probing_scheme_type,
                                             storage_ref_type,
                                             Operators...>;  ///< Non-owning container ref type

  static_multiset(static_multiset const&)            = delete;
  static_multiset& operator=(static_multiset const&) = delete;

  static_multiset(static_multiset&&) = default;  ///< Move constructor

  /**
   * @brief Replaces the contents of the container with another container.
   *
   * @return Reference of the current multiset object
   */
  static_multiset& operator=(static_multiset&&) = default;
  ~static_multiset()                            = default;

  /**
   * @brief Constructs a statically-sized multiset with the specified initial capacity, sentinel
   * values and CUDA stream
   *
   * The actual multiset capacity depends on the given `capacity`, the probing scheme, CG size, and
   * the bucket size and it is computed via the `make_valid_extent` factory. Insert operations will
   * not automatically grow the set. Attempting to insert more unique keys than the capacity of the
   * multiset results in undefined behavior.
   *
   * @note Any `*_sentinel`s are reserved and behavior is undefined when attempting to insert
   * this sentinel value.
   * @note This constructor doesn't synchronize the given stream.
   *
   * @param capacity The requested lower-bound multiset size
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage Kind of storage to use
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the set
   */
  constexpr static_multiset(Extent capacity,
                            empty_key<Key> empty_key_sentinel,
                            KeyEqual const& pred                = {},
                            ProbingScheme const& probing_scheme = {},
                            cuda_thread_scope<Scope> scope      = {},
                            Storage storage                     = {},
                            Allocator const& alloc              = {},
                            cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Constructs a statically-sized multiset with the number of elements to insert `n`, the
   * desired load factor, etc
   *
   * @note This constructor helps users create a set based on the number of elements to insert and
   * the desired load factor without manually computing the desired capacity. The actual set
   * capacity will be a size no smaller than `ceil(n / desired_load_factor)`. It's determined by
   * multiple factors including the given `n`, the desired load factor, the probing scheme, the CG
   * size, and the bucket size and is computed via the `make_valid_extent` factory.
   * @note Insert operations will not automatically grow the container.
   * @note Attempting to insert more unique keys than the capacity of the container results in
   * undefined behavior.
   * @note Any `*_sentinel`s are reserved and behavior is undefined when attempting to insert
   * this sentinel value.
   * @note This constructor doesn't synchronize the given stream.
   * @note This overload will convert compile-time extents to runtime constants which might lead to
   * performance regressions.
   *
   * @throw If the desired occupancy is no bigger than zero
   * @throw If the desired occupancy is no smaller than one
   *
   * @param n The number of elements to insert
   * @param desired_load_factor The desired load factor of the container, e.g., 0.5 implies a 50%
   * load factor
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage Kind of storage to use
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the set
   */
  constexpr static_multiset(Extent n,
                            double desired_load_factor,
                            empty_key<Key> empty_key_sentinel,
                            KeyEqual const& pred                = {},
                            ProbingScheme const& probing_scheme = {},
                            cuda_thread_scope<Scope> scope      = {},
                            Storage storage                     = {},
                            Allocator const& alloc              = {},
                            cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Constructs a statically-sized set with the specified initial capacity, sentinel values
   * and CUDA stream.
   *
   * The actual set capacity depends on the given `capacity`, the probing scheme, CG size, and the
   * bucket size and it is computed via the `make_valid_extent` factory. Insert operations will not
   * automatically grow the set. Attempting to insert more unique keys than the capacity of the
   * multiset results in undefined behavior.
   *
   * @note Any `*_sentinel`s are reserved and behavior is undefined when attempting to insert
   * this sentinel value.
   * @note If a non-default CUDA stream is provided, the caller is responsible for synchronizing the
   * stream before the object is first used.
   *
   * @param capacity The requested lower-bound set size
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param erased_key_sentinel The reserved key to denote erased slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage Kind of storage to use
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the set
   */
  constexpr static_multiset(Extent capacity,
                            empty_key<Key> empty_key_sentinel,
                            erased_key<Key> erased_key_sentinel,
                            KeyEqual const& pred                = {},
                            ProbingScheme const& probing_scheme = {},
                            cuda_thread_scope<Scope> scope      = {},
                            Storage storage                     = {},
                            Allocator const& alloc              = {},
                            cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Erases all elements from the container. After this call, `size()` returns zero.
   * Invalidates any references, pointers, or iterators referring to contained elements.
   *
   * @param stream CUDA stream this operation is executed in
   */
  void clear(cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously erases all elements from the container. After this call, `size()` returns
   * zero. Invalidates any references, pointers, or iterators referring to contained elements.
   *
   * @param stream CUDA stream this operation is executed in
   */
  void clear_async(cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief Inserts all keys in the range `[first, last)`
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `insert_async`.
   *
   * // TODO: to be revised due to heterogeneous lookup
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multiset<K>::value_type></tt> is `true`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt>
  void insert(InputIt first,
              InputIt last,
              cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously inserts all keys in the range `[first, last)`.
   *
   * // TODO: to be revised due to heterogeneous lookup
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multiset<K>::value_type></tt> is `true`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt>
  void insert_async(InputIt first,
                    InputIt last,
                    cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief Inserts keys in the range `[first, last)` if `pred` of the corresponding stencil returns
   * true.
   *
   * @note The key `*(first + i)` is inserted if `pred( *(stencil + i) )` returns true.
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `insert_if_async`.
   *
   * @tparam InputIt Device accessible random access iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   *
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param stream CUDA stream used for the operation
   *
   * @return Number of successful insertions
   */
  template <typename InputIt, typename StencilIt, typename Predicate>
  size_type insert_if(InputIt first,
                      InputIt last,
                      StencilIt stencil,
                      Predicate pred,
                      cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously inserts keys in the range `[first, last)` if `pred` of the corresponding
   * stencil returns true.
   *
   * @note The key `*(first + i)` is inserted if `pred( *(stencil + i) )` returns true.
   *
   * @tparam InputIt Device accessible random access iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   *
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param stream CUDA stream used for the operation
   */
  template <typename InputIt, typename StencilIt, typename Predicate>
  void insert_if_async(InputIt first,
                       InputIt last,
                       StencilIt stencil,
                       Predicate pred,
                       cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief Indicates whether the keys in the range `[first, last)` are contained in the multiset.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `contains_async`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename OutputIt>
  void contains(InputIt first,
                InputIt last,
                OutputIt output_begin,
                cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Asynchronously indicates whether the keys in the range `[first, last)` are contained in
   * the multiset.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename OutputIt>
  void contains_async(InputIt first,
                      InputIt last,
                      OutputIt output_begin,
                      cuda::stream_ref stream = cuda::stream_ref{
                        cudaStream_t{nullptr}}) const noexcept;

  /**
   * @brief Indicates whether the keys in the range `[first, last)` are contained in the multiset if
   * `pred` of the corresponding stencil returns `true`.
   *
   * @note If `pred( *(stencil + i) )` is true, stores `true` or `false` to `(output_begin + i)`
   * indicating if the key `*(first + i)` is present in the multiset. If `pred( *(stencil + i) )` is
   * `false`, stores `false` to `(output_begin + i)`.
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `contains_if_async`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose value type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
  void contains_if(InputIt first,
                   InputIt last,
                   StencilIt stencil,
                   Predicate pred,
                   OutputIt output_begin,
                   cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Asynchronously indicates whether the keys in the range `[first, last)` are contained in
   * the multiset if `pred` of the corresponding stencil returns `true`.
   *
   * @note If `pred( *(stencil + i) )` is true, stores `true` or `false` to `(output_begin + i)`
   * indicating if the key `*(first + i)` is present in the multiset. If `pred( *(stencil + i) )` is
   * `false`, stores `false` to `(output_begin + i)`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose value type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
  void contains_if_async(InputIt first,
                         InputIt last,
                         StencilIt stencil,
                         Predicate pred,
                         OutputIt output_begin,
                         cuda::stream_ref stream = cuda::stream_ref{
                           cudaStream_t{nullptr}}) const noexcept;

  /**
   * @brief For all keys in the range `[first, last)`, finds an element with its key equivalent to
   * the query key.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use `find_async`.
   * @note If the key `*(first + i)` has a matched `element` in the multiset, copies `element` to
   * `(output_begin + i)`. Else, copies the empty key sentinel.
   * @note For a given key `*(first + i)`, if there are multiple matching elements in the multiset,
   * it copies the payload of one match (unspecified which) to `(output_begin + i)`. If no match is
   * found, it copies the empty key sentinel instead.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from the set's `key_type`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of elements retrieved for each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename OutputIt>
  void find(InputIt first,
            InputIt last,
            OutputIt output_begin,
            cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief For all keys in the range `[first, last)`, asynchronously finds an element with its key
   * equivalent to the query key.
   *
   * @note If the key `*(first + i)` has a matched `element` in the multiset, copies `element` to
   * `(output_begin + i)`. Else, copies the empty key sentinel.
   * @note For a given key `*(first + i)`, if there are multiple matching elements in the multiset,
   * it copies the payload of one match (unspecified which) to `(output_begin + i)`. If no match is
   * found, it copies the empty key sentinel instead.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from the set's `key_type`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of elements retrieved for each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename OutputIt>
  void find_async(InputIt first,
                  InputIt last,
                  OutputIt output_begin,
                  cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief For all keys in the range `[first, last)`, finds a match with its key equivalent to the
   * query key.
   *
   * @note If `pred( *(stencil + i) )` is true, stores the payload of the
   * matched key or the `empty_value_sentinel` to `(output_begin + i)`. If `pred( *(stencil + i) )`
   * is false, always stores the `empty_value_sentinel` to `(output_begin + i)`.
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `find_if_async`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose `value_type` is convertible to
   * Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam OutputIt Device accessible output iterator
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param output_begin Beginning of the sequence of matches retrieved for each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
  void find_if(InputIt first,
               InputIt last,
               StencilIt stencil,
               Predicate pred,
               OutputIt output_begin,
               cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief For all keys in the range `[first, last)`, asynchronously finds
   * a match with its key equivalent to the query key.
   *
   * @note If `pred( *(stencil + i) )` is true, stores the payload of the
   * matched key or the `empty_value_sentinel` to `(output_begin + i)`. If `pred( *(stencil + i) )`
   * is false, always stores the `empty_value_sentinel` to `(output_begin + i)`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose `value_type` is convertible to
   * Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam OutputIt Device accessible output iterator
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param output_begin Beginning of the sequence of matches retrieved for each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename StencilIt, typename Predicate, typename OutputIt>
  void find_if_async(InputIt first,
                     InputIt last,
                     StencilIt stencil,
                     Predicate pred,
                     OutputIt output_begin,
                     cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Applies the given function object `callback_op` to the copy of every filled slot in the
   * container
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @tparam CallbackOp Type of unary callback function object
   *
   * @param callback_op Function to apply to the copy of the filled slot
   * @param stream CUDA stream used for this operation
   */
  template <typename CallbackOp>
  void for_each(CallbackOp&& callback_op,
                cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Asynchronously applies the given function object `callback_op` to the copy of every
   * filled slot in the container
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @tparam CallbackOp Type of unary callback function object
   *
   * @param callback_op Function to apply to the copy of the filled slot
   * @param stream CUDA stream used for this operation
   */
  template <typename CallbackOp>
  void for_each_async(CallbackOp&& callback_op,
                      cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief For each key in the range [first, last), applies the function object `callback_op` to
   * the copy of all corresponding matches found in the container.
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @tparam InputIt Device accessible random access input iterator
   * @tparam CallbackOp Type of unary callback function object
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param callback_op Function to apply to the copy of the matched slot
   * @param stream CUDA stream used for this operation
   */
  template <typename InputIt, typename CallbackOp>
  void for_each(InputIt first,
                InputIt last,
                CallbackOp&& callback_op,
                cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief For each key in the range [first, last), asynchronously applies the function object
   * `callback_op` to the copy of all corresponding matches found in the container.
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @tparam InputIt Device accessible random access input iterator
   * @tparam CallbackOp Type of unary callback function object
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param callback_op Function to apply to the copy of the matched slot
   * @param stream CUDA stream used for this operation
   */
  template <typename InputIt, typename CallbackOp>
  void for_each_async(InputIt first,
                      InputIt last,
                      CallbackOp&& callback_op,
                      cuda::stream_ref stream = cuda::stream_ref{
                        cudaStream_t{nullptr}}) const noexcept;

  /**
   * @brief Counts the occurrences of keys in `[first, last)` contained in the multiset
   *
   * @note This function synchronizes the given stream.
   *
   * @tparam Input Device accessible input iterator
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param stream CUDA stream used for count
   *
   * @return The sum of total occurrences of all keys in `[first, last)`
   */
  template <typename InputIt>
  size_type count(InputIt first,
                  InputIt last,
                  cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Counts the occurrences of keys in `[first, last)` contained in the multiset
   *
   * @note This function synchronizes the given stream.
   *
   * @tparam Input Device accessible input iterator
   * @tparam ProbeKeyEqual Binary callable
   * @tparam ProbeHash Unary hash callable
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param probe_key_equal Binary callable to compare two keys for equality
   * @param probe_hash Unary callable to hash a given key
   * @param stream CUDA stream used for count
   *
   * @return The sum of total occurrences of all keys in `[first, last)`
   */
  template <typename InputIt, typename ProbeKeyEqual, typename ProbeHash>
  size_type count(InputIt first,
                  InputIt last,
                  ProbeKeyEqual const& probe_key_equal,
                  ProbeHash const& probe_hash,
                  cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Counts the occurrences of keys in `[first, last)` contained in the multiset
   *
   * @note This function synchronizes the given stream.
   * @note If a given key has no matches, its occurrence is 1.
   *
   * @tparam Input Device accessible input iterator
   * @tparam ProbeKeyEqual Binary callable
   * @tparam ProbeHash Unary hash callable
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param probe_key_equal Binary callable to compare two keys for equality
   * @param probe_hash Unary callable to hash a given key
   * @param stream CUDA stream used for count
   *
   * @return The sum of total occurrences of all keys in `[first, last)` where keys have no matches
   * are considered to have a single occurrence.
   */
  template <typename InputIt, typename ProbeKeyEqual, typename ProbeHash>
  size_type count_outer(InputIt first,
                        InputIt last,
                        ProbeKeyEqual const& probe_key_equal,
                        ProbeHash const& probe_hash,
                        cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Counts the number of occurrences of each query key in the multiset
   *
   * For each key in the input range `[first, last)`, this function computes the number of matching
   * elements in the multiset and writes the result to the corresponding position in the output
   * range starting at `output_begin`.
   *
   * @note The input and output ranges must be device-accessible and of the same length.
   * @note The behavior is undefined if the input and output ranges overlap.
   *
   * @tparam InputIt       Device-accessible input iterator type for query keys
   * @tparam ProbeKeyEqual Binary callable that compares two keys for equality
   * @tparam ProbeHash     Unary callable that computes the hash of a key
   * @tparam OutputIt      Device-accessible output iterator type for storing per-key counts
   *
   * @param first          Iterator to the beginning of the sequence of query keys
   * @param last           Iterator to the end of the sequence of query keys
   * @param probe_key_equal Predicate to compare a query key with a multiset key for equality
   * @param probe_hash     Hash function to compute the hash value of a query key
   * @param output_begin   Iterator to the beginning of the output range where per-key counts will
   * be stored
   * @param stream         CUDA stream on which to execute the counting operation
   */
  template <typename InputIt, typename ProbeKeyEqual, typename ProbeHash, typename OutputIt>
  void count_each(InputIt first,
                  InputIt last,
                  ProbeKeyEqual const& probe_key_equal,
                  ProbeHash const& probe_hash,
                  OutputIt output_begin,
                  cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Counts the number of occurrences of each query key in the multiset with outer semantics.
   *
   * For each key in the input range `[first, last)`, this function computes the number of matching
   * elements in the multiset and writes the result to the corresponding position in the output
   * range starting at `output_begin`.
   *
   * If a query key has no matches in the multiset, the result for that key will be 1 instead of 0.
   * Otherwise, the actual number of matches is returned.
   *
   * This provides "outer join"-like semantics, ensuring that every query key contributes at least 1
   * count.
   *
   * @note The input and output ranges must be device-accessible and of the same length.
   * @note The behavior is undefined if the input and output ranges overlap.
   *
   * @tparam InputIt       Device-accessible input iterator type for query keys
   * @tparam ProbeKeyEqual Binary callable that compares two keys for equality
   * @tparam ProbeHash     Unary callable that computes the hash of a key
   * @tparam OutputIt      Device-accessible output iterator type for storing per-key counts
   *
   * @param first          Iterator to the beginning of the sequence of query keys
   * @param last           Iterator to the end of the sequence of query keys
   * @param probe_key_equal Predicate to compare a query key with a multiset key for equality
   * @param probe_hash     Hash function to compute the hash value of a query key
   * @param output_begin   Iterator to the beginning of the output range where per-key counts will
   * be stored
   * @param stream         CUDA stream on which to execute the counting operation
   */
  template <typename InputIt, typename ProbeKeyEqual, typename ProbeHash, typename OutputIt>
  void count_each_outer(InputIt first,
                        InputIt last,
                        ProbeKeyEqual const& probe_key_equal,
                        ProbeHash const& probe_hash,
                        OutputIt output_begin,
                        cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Retrieves all the slots corresponding to all keys in the range `[first, last)`.
   *
   * If key `k = *(first + i)` exists in the container, copies `k` to `output_probe` and associated
   * slot contents to `output_match`, respectively. The output order is unspecified.
   *
   * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
   * Use `count()` to determine the size of the output range.
   *
   * This function synchronizes the given CUDA stream.
   *
   * @tparam InputProbeIt Device accessible input iterator
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the `InputProbeIt`'s `value_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   *
   * @param first Beginning of the input sequence of keys
   * @param last End of the input sequence of keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param stream CUDA stream this operation is executed in
   *
   * @return Iterator pair indicating the the end of the output sequences
   */
  template <class InputProbeIt, class OutputProbeIt, class OutputMatchIt>
  std::pair<OutputProbeIt, OutputMatchIt> retrieve(InputProbeIt first,
                                                   InputProbeIt last,
                                                   OutputProbeIt output_probe,
                                                   OutputMatchIt output_match,
                                                   cuda::stream_ref stream = cuda::stream_ref{
                                                     cudaStream_t{nullptr}}) const;

  /**
   * @brief Retrieves all the slots corresponding to all keys in the range `[first, last)`.
   *
   * If key `k = *(first + i)` exists in the container, copies `k` to `output_probe` and associated
   * slot contents to `output_match`, respectively. The output order is unspecified.
   *
   * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
   * Use `count()` to determine the size of the output range.
   *
   * This function synchronizes the given CUDA stream.
   *
   * @tparam InputProbeIt Device accessible input iterator
   * @tparam ProbeEqual Binary callable equal type
   * @tparam ProbeHash Unary callable hasher type that can be constructed from
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the `InputProbeIt`'s `value_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   *
   * @param first Beginning of the input sequence of keys
   * @param last End of the input sequence of keys
   * @param probe_equal The binary function to compare set keys and probe keys for equality
   * @param probe_hash The unary function to hash probe keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param stream CUDA stream this operation is executed in
   *
   * @return Iterator pair indicating the the end of the output sequences
   */
  template <class InputProbeIt,
            class ProbeEqual,
            class ProbeHash,
            class OutputProbeIt,
            class OutputMatchIt>
  std::pair<OutputProbeIt, OutputMatchIt> retrieve(InputProbeIt first,
                                                   InputProbeIt last,
                                                   ProbeEqual const& probe_equal,
                                                   ProbeHash const& probe_hash,
                                                   OutputProbeIt output_probe,
                                                   OutputMatchIt output_match,
                                                   cuda::stream_ref stream = cuda::stream_ref{
                                                     cudaStream_t{nullptr}}) const;

  /**
   * @brief Retrieves all the slots corresponding to all keys in the range `[first, last)`.
   *
   * If key `k = *(first + i)` exists in the container, copies `k` to `output_probe` and associated
   * slot contents to `output_match`, respectively. The output order is unspecified.
   *
   * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
   * Use `count_outer()` to determine the size of the output range.
   *
   * If a key `k` has no matches in the container, then `{key, empty_slot_sentinel}` will be added
   * to the output sequence.
   *
   * This function synchronizes the given CUDA stream.
   *
   * @tparam InputProbeIt Device accessible input iterator
   * @tparam ProbeEqual Binary callable equal type
   * @tparam ProbeHash Unary callable hasher type that can be constructed from
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the `InputProbeIt`'s `value_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   *
   * @param first Beginning of the input sequence of keys
   * @param last End of the input sequence of keys
   * @param probe_equal The binary function to compare set keys and probe keys for equality
   * @param probe_hash The unary function to hash probe keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param stream CUDA stream this operation is executed in
   *
   * @return Iterator pair indicating the the end of the output sequences
   */
  template <class InputProbeIt,
            class ProbeEqual,
            class ProbeHash,
            class OutputProbeIt,
            class OutputMatchIt>
  std::pair<OutputProbeIt, OutputMatchIt> retrieve_outer(InputProbeIt first,
                                                         InputProbeIt last,
                                                         ProbeEqual const& probe_equal,
                                                         ProbeHash const& probe_hash,
                                                         OutputProbeIt output_probe,
                                                         OutputMatchIt output_match,
                                                         cuda::stream_ref stream = cuda::stream_ref{
                                                           cudaStream_t{nullptr}}) const;

  /**
   * @brief Retrieves all keys contained in the multiset
   *
   * @note This API synchronizes the given stream.
   * @note The order in which keys are returned is implementation defined and not guaranteed to be
   * consistent between subsequent calls to `retrieve_all`.
   * @note Behavior is undefined if the range beginning at `output_begin` is smaller than the return
   * value of `size()`.
   *
   * @tparam OutputIt Device accessible random access output iterator whose `value_type` is
   * convertible from the container's `key_type`.
   *
   * @param output_begin Beginning output iterator for keys
   * @param stream CUDA stream used for this operation
   *
   * @return Iterator indicating the end of the output
   */
  template <typename OutputIt>
  OutputIt retrieve_all(OutputIt output_begin,
                        cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Regenerates the container.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `rehash_async`.
   *
   * @param stream CUDA stream used for this operation
   */
  void rehash(cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Reserves at least the specified number of slots and regenerates the container
   *
   * @note Changes the number of slots to a value that is not less than `capacity`, then
   * rehashes the container, i.e. puts the elements into appropriate slots considering
   * that the total number of slots has changed.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `rehash_async`.
   *
   * @note Behavior is undefined if the desired `capacity` is insufficient to store all of the
   * contained elements.
   *
   * @note This function is not available if the container's `extent_type` is static.
   *
   * @param capacity New capacity of the container
   * @param stream CUDA stream used for this operation
   */
  void rehash(size_type capacity,
              cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously regenerates the container.
   *
   * @param stream CUDA stream used for this operation
   */
  void rehash_async(cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously reserves at least the specified number of slots and regenerates the
   * container
   *
   * @note Changes the number of slots to a value that is not less than `capacity`, then
   * rehashes the container, i.e. puts the elements into appropriate slots considering
   * that the total number of slots has changed.
   *
   * @note Behavior is undefined if the desired `capacity` is insufficient to store all of the
   * contained elements.
   *
   * @note This function is not available if the container's `extent_type` is static.
   *
   * @param capacity New capacity of the container
   * @param stream CUDA stream used for this operation
   */
  void rehash_async(size_type capacity,
                    cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Gets the number of elements in the container.
   *
   * @note This function synchronizes the given stream.
   *
   * @param stream CUDA stream used to get the number of inserted elements
   * @return The number of elements in the container
   */
  [[nodiscard]] size_type size(cuda::stream_ref stream = cuda::stream_ref{
                                 cudaStream_t{nullptr}}) const;

  /**
   * @brief Gets the maximum number of elements the multiset can hold.
   *
   * @return The maximum number of elements the multiset can hold
   */
  [[nodiscard]] constexpr auto capacity() const noexcept;

  /**
   * @brief Gets a pointer to the underlying slot storage.
   *
   * @return Pointer to the underlying slot storage
   */
  [[nodiscard]] __host__ value_type* data() const;

  /**
   * @brief Gets the sentinel value used to represent an empty key slot.
   *
   * @return The sentinel value used to represent an empty key slot
   */
  [[nodiscard]] constexpr key_type empty_key_sentinel() const noexcept;

  /**
   * @brief Gets the sentinel value used to represent an erased key slot.
   *
   * @return The sentinel value used to represent an erased key slot
   */
  [[nodiscard]] constexpr key_type erased_key_sentinel() const noexcept;

  /**
   * @brief Gets the function used to compare keys for equality
   *
   * @return The function used to compare keys for equality
   */
  [[nodiscard]] constexpr key_equal key_eq() const noexcept;

  /**
   * @brief Gets the function(s) used to hash keys
   *
   * @return The function(s) used to hash keys
   */
  [[nodiscard]] constexpr hasher hash_function() const noexcept;

  /**
   * @brief Get device ref with operators.
   *
   * @tparam Operators Set of `cuco::op` to be provided by the ref
   *
   * @param ops List of operators, e.g., `cuco::insert`
   *
   * @return Device ref of the current `static_multiset` object
   */
  template <typename... Operators>
  [[nodiscard]] auto ref(Operators... ops) const noexcept;

 private:
  std::unique_ptr<impl_type> impl_;
};
}  // namespace cuco

#include <cuco/detail/static_multiset/static_multiset.inl>
