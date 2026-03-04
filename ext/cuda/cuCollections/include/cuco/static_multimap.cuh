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
 * limitations under the License.
 */

#pragma once

#include <cuco/detail/__config>
#include <cuco/detail/open_addressing/open_addressing_impl.cuh>
#include <cuco/detail/prime.hpp>
#include <cuco/hash_functions.cuh>
#include <cuco/probe_sequences.cuh>
#include <cuco/static_multimap_ref.cuh>
#include <cuco/types.cuh>
#include <cuco/utility/allocator.hpp>
#include <cuco/utility/traits.hpp>

#include <cuda/std/atomic>
#include <cuda/std/functional>
#include <cuda/stream_ref>

#if defined(CUCO_HAS_CUDA_BARRIER)
#include <cuda/barrier>
#endif

#include <cooperative_groups.h>

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace cuco {
namespace experimental {
/**
 * @brief A GPU-accelerated, unordered, associative container of key-value pairs that supports
 * equivalent keys.
 *
 * The `static_multimap` supports two types of operations:
 * - Host-side "bulk" operations
 * - Device-side "singular" operations
 *
 * The host-side bulk operations include `insert`, `contains`, etc. These APIs should be used when
 * there are a large number of keys to modify or lookup. For example, given a range of keys
 * specified by device-accessible iterators, the bulk `insert` function will insert all keys into
 * the map.
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
 * @throw If the size of the given payload type is larger than 8 bytes
 * @throw If the size of the given slot type is larger than 16 bytes
 * @throw If the given key type doesn't have unique object representations, i.e.,
 * `cuco::bitwise_comparable_v<Key> == false`
 * @throw If the given mapped type doesn't have unique object representations, i.e.,
 * `cuco::bitwise_comparable_v<T> == false`
 * @throw If the probing scheme type is not inherited from `cuco::detail::probing_scheme_base`
 *
 * @tparam Key Type used for keys. Requires `cuco::is_bitwise_comparable_v<Key>`
 * @tparam T Type of the mapped values
 * @tparam Extent Data structure size type
 * @tparam Scope The scope in which operations will be performed by individual threads.
 * @tparam KeyEqual Binary callable type used to compare two keys for equality
 * @tparam ProbingScheme Probing scheme (see `include/cuco/probing_scheme.cuh` for choices)
 * @tparam Allocator Type of allocator used for device storage
 * @tparam Storage Slot bucket storage type
 */
template <class Key,
          class T,
          class Extent             = cuco::extent<std::size_t>,
          cuda::thread_scope Scope = cuda::thread_scope_device,
          class KeyEqual           = cuda::std::equal_to<Key>,
          class ProbingScheme      = cuco::double_hashing<8,  // CG size
                                                          cuco::default_hash_function<Key>>,
          class Allocator          = cuco::cuda_allocator<cuco::pair<Key, T>>,
          class Storage            = cuco::storage<2>>
class static_multimap {
  static_assert(sizeof(Key) <= 8, "Container does not support key types larger than 8 bytes.");

  static_assert(sizeof(T) <= 8, "Container does not support payload types larger than 8 bytes.");

  static_assert(cuco::is_bitwise_comparable_v<T>,
                "Mapped type must have unique object representations or have been explicitly "
                "declared as safe for bitwise comparison via specialization of "
                "cuco::is_bitwise_comparable_v<T>.");

  using impl_type = cuco::detail::open_addressing_impl<Key,
                                                       cuco::pair<Key, T>,
                                                       Extent,
                                                       Scope,
                                                       KeyEqual,
                                                       ProbingScheme,
                                                       Allocator,
                                                       Storage>;

 public:
  static constexpr auto cg_size      = impl_type::cg_size;       ///< CG size used for probing
  static constexpr auto bucket_size  = impl_type::bucket_size;   ///< Bucket size used for probing
  static constexpr auto thread_scope = impl_type::thread_scope;  ///< CUDA thread scope

  using key_type       = typename impl_type::key_type;        ///< Key type
  using value_type     = typename impl_type::value_type;      ///< Key-value pair type
  using extent_type    = typename impl_type::extent_type;     ///< Extent type
  using size_type      = typename impl_type::size_type;       ///< Size type
  using key_equal      = typename impl_type::key_equal;       ///< Key equality comparator type
  using allocator_type = typename impl_type::allocator_type;  ///< Allocator type
  /// Non-owning bucket storage ref type
  using storage_ref_type    = typename impl_type::storage_ref_type;
  using probing_scheme_type = typename impl_type::probing_scheme_type;  ///< Probing scheme type
  using hasher              = typename probing_scheme_type::hasher;     ///< Hash function type

  using mapped_type = T;  ///< Payload type
  template <typename... Operators>
  using ref_type = cuco::static_multimap_ref<key_type,
                                             mapped_type,
                                             thread_scope,
                                             key_equal,
                                             probing_scheme_type,
                                             storage_ref_type,
                                             Operators...>;  ///< Non-owning container ref type

  static_multimap(static_multimap const&)            = delete;
  static_multimap& operator=(static_multimap const&) = delete;

  static_multimap(static_multimap&&) = default;  ///< Move constructor

  /**
   * @brief Replaces the contents of the container with another container.
   *
   * @return Reference of the current map object
   */
  static_multimap& operator=(static_multimap&&) = default;
  ~static_multimap()                            = default;

  /**
   * @brief Constructs a statically-sized map with the specified initial capacity, sentinel values
   * and CUDA stream
   *
   * The actual map capacity depends on the given `capacity`, the probing scheme, CG size, and the
   * bucket size and it is computed via the `make_valid_extent` factory. Insert operations will not
   * automatically grow the map. Attempting to insert more unique keys than the capacity of the map
   * results in undefined behavior.
   *
   * @note Any `*_sentinel`s are reserved and behavior is undefined when attempting to insert
   * this sentinel value.
   * @note This constructor doesn't synchronize the given stream.
   *
   * @param capacity The requested lower-bound map size
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param empty_value_sentinel The reserved mapped value for empty slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage Kind of storage to use
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the map
   */
  constexpr static_multimap(Extent capacity,
                            empty_key<Key> empty_key_sentinel,
                            empty_value<T> empty_value_sentinel,
                            KeyEqual const& pred                = {},
                            ProbingScheme const& probing_scheme = {},
                            cuda_thread_scope<Scope> scope      = {},
                            Storage storage                     = {},
                            Allocator const& alloc              = {},
                            cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Constructs a statically-sized map with the number of elements to insert `n`, the desired
   * load factor, etc
   *
   * @note This constructor helps users create a map based on the number of elements to insert and
   * the desired load factor without manually computing the desired capacity. The actual map
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
   * @param empty_value_sentinel The reserved mapped value for empty slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage Kind of storage to use
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the map
   */
  constexpr static_multimap(Extent n,
                            double desired_load_factor,
                            empty_key<Key> empty_key_sentinel,
                            empty_value<T> empty_value_sentinel,
                            KeyEqual const& pred                = {},
                            ProbingScheme const& probing_scheme = {},
                            cuda_thread_scope<Scope> scope      = {},
                            Storage storage                     = {},
                            Allocator const& alloc              = {},
                            cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Constructs a statically-sized map with the specified initial capacity, sentinel values
   * and CUDA stream.
   *
   * The actual map capacity depends on the given `capacity`, the probing scheme, CG size, and the
   * bucket size and it is computed via the `make_valid_extent` factory. Insert operations will not
   * automatically grow the map. Attempting to insert more unique keys than the capacity of the map
   * results in undefined behavior.
   *
   * @note Any `*_sentinel`s are reserved and behavior is undefined when attempting to insert
   * this sentinel value.
   * @note If a non-default CUDA stream is provided, the caller is responsible for synchronizing the
   * stream before the object is first used.
   *
   * @param capacity The requested lower-bound map size
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param empty_value_sentinel The reserved mapped value for empty slots
   * @param erased_key_sentinel The reserved key to denote erased slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage Kind of storage to use
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the map
   */
  constexpr static_multimap(Extent capacity,
                            empty_key<Key> empty_key_sentinel,
                            empty_value<T> empty_value_sentinel,
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
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multimap<K, V>::value_type></tt> is `true`
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
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multimap<K, V>::value_type></tt> is `true`
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
   * @note This function synchronizes the given stream and returns the number of successful
   * insertions. For asynchronous execution use `insert_if_async`.
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
   * @brief Indicates whether the keys in the range `[first, last)` are contained in the map.
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
   * the map.
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
   * @brief Indicates whether the keys in the range `[first, last)` are contained in the map if
   * `pred` of the corresponding stencil returns true.
   *
   * @note If `pred( *(stencil + i) )` is true, stores `true` or `false` to `(output_begin + i)`
   * indicating if the key `*(first + i)` is present in the map. If `pred( *(stencil + i) )` is
   * false, stores `false` to `(output_begin + i)`.
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `contains_if_async`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose value_type is
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
   * the map if `pred` of the corresponding stencil returns true.
   *
   * @note If `pred( *(stencil + i) )` is true, stores `true` or `false` to `(output_begin + i)`
   * indicating if the key `*(first + i)` is present in the map. If `pred( *(stencil + i) )` is
   * false, stores `false` to `(output_begin + i)`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose value_type is
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
   * @brief For all keys in the range `[first, last)`, finds a payload with its key equivalent to
   * the query key.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use `find_async`.
   * @note If the key `*(first + i)` has a matched `element` in the map, copies the payload of
   * `element` to `(output_begin + i)`. Else, copies the empty value sentinel.
   * @note For a given key `*(first + i)`, if there are multiple matching elements in the multimap,
   * it copies the payload of one match (unspecified which) to `(output_begin + i)`. If no match is
   * found, it copies the empty value sentinel instead.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from the map's `mapped_type`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of payloads retrieved for each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename OutputIt>
  void find(InputIt first,
            InputIt last,
            OutputIt output_begin,
            cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief For all keys in the range `[first, last)`, asynchronously finds a payload with its key
   * equivalent to the query key.
   *
   * @note If the key `*(first + i)` has a matched `element` in the map, copies the payload of
   * `element` to `(output_begin + i)`. Else, copies the empty value sentinel.
   * @note For a given key `*(first + i)`, if there are multiple matching elements in the multimap,
   * it copies the payload of one match (unspecified which) to `(output_begin + i)`. If no match is
   * found, it copies the empty value sentinel instead.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from the map's `mapped_type`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of payloads retrieved for each key
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
   * @brief Counts the occurrences of keys in `[first, last)` contained in the multimap
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
   * @brief Counts the occurrences of keys in `[first, last)` contained in the
   * multimap with custom hash and key comparator
   *
   * @note This function synchronizes the given stream.
   *
   * @tparam Input Device accessible input iterator
   * @tparam ProbeEqual Binary callable
   * @tparam ProbeHash Unary hash callable
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param probe_equal Binary callable to compare two keys for equality
   * @param probe_hash Unary callable to hash a given key
   * @param stream CUDA stream used for count
   *
   * @return The sum of total occurrences of all keys in `[first, last)`
   */
  template <typename InputIt, typename ProbeEqual, typename ProbeHash>
  size_type count(InputIt first,
                  InputIt last,
                  ProbeEqual const& probe_equal,
                  ProbeHash const& probe_hash,
                  cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}}) const;

  /**
   * @brief Retrieves the matched key-value pair in the multimap corresponding to all probe keys in
   * the range `[first, last)`
   *
   * If key `k = *(first + i)` has a match `m` in the multimap, copies a `cuco::pair{k, m}` to
   * unspecified locations in `[output_begin, output_end)`. Else, does nothing.
   *
   * @note This function synchronizes the given stream.
   * @note Behavior is undefined if the size of the output range exceeds
   * `std::distance(output_begin, output_end)`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputProbeIt Device accessible output iterator whose `value_type` can be constructed
   * from `ProbeKey`
   * @tparam OutputMatchIt Device accessible output iterator whose `value_type` can be constructed
   * from multimap's `value_type`
   *
   * @param first Beginning of the sequence of probe keys
   * @param last End of the sequence of probe keys
   * @param output_probe Beginning of the sequence of the probe keys that have a match
   * @param output_match Beginning of the sequence of the matched key-value pairs
   * @param stream CUDA stream used for retrieve
   *
   * @return The iterator indicating the last valid pair in the output
   */
  template <typename InputIt, typename OutputProbeIt, typename OutputMatchIt>
  std::pair<OutputProbeIt, OutputMatchIt> retrieve(InputIt first,
                                                   InputIt last,
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
   * @brief Retrieves all of the keys and their associated values contained in the multimap
   *
   * @note This API synchronizes the given stream.
   * @note The order in which keys are returned is implementation defined and not guaranteed to be
   * consistent between subsequent calls to `retrieve_all`.
   * @note Behavior is undefined if the range beginning at `keys_out` or `values_out` is smaller
   * than the return value of `size()`.
   *
   * @tparam KeyOut Device accessible random access output iterator whose `value_type` is
   * convertible from `key_type`.
   * @tparam ValueOut Device accessible random access output iterator whose `value_type` is
   * convertible from `mapped_type`.
   *
   * @param keys_out Beginning output iterator for keys
   * @param values_out Beginning output iterator for associated values
   * @param stream CUDA stream used for this operation
   *
   * @return Pair of iterators indicating the last elements in the output
   */
  template <typename KeyOut, typename ValueOut>
  std::pair<KeyOut, ValueOut> retrieve_all(KeyOut keys_out,
                                           ValueOut values_out,
                                           cuda::stream_ref stream = cuda::stream_ref{
                                             cudaStream_t{nullptr}}) const;

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
   * @brief Gets the maximum number of elements the multimap can hold.
   *
   * @return The maximum number of elements the multimap can hold
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
   * @brief Gets the sentinel value used to represent an empty value slot.
   *
   * @return The sentinel value used to represent an empty value slot
   */
  [[nodiscard]] constexpr mapped_type empty_value_sentinel() const noexcept;

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
   * @return Device ref of the current `static_multimap` object
   */
  template <typename... Operators>
  [[nodiscard]] auto ref(Operators... ops) const noexcept;

 private:
  std::unique_ptr<impl_type> impl_;   ///< Static map implementation
  mapped_type empty_value_sentinel_;  ///< Sentinel value that indicates an empty payload
};
}  // namespace experimental

/**
 * @brief A GPU-accelerated, unordered, associative container of key-value pairs that supports
 * equivalent keys.
 *
 * Allows constant time concurrent inserts or concurrent find operations from threads in device
 * code. Concurrent insert/find is allowed only when
 * <tt>static_multimap<Key, Value>::supports_concurrent_insert_find()</tt> is true.
 *
 * Current limitations:
 * - Requires keys and values where `cuco::is_bitwise_comparable_v<T>` is true
 * - Comparisons against the "sentinel" values will always be done with bitwise comparisons
 * Therefore, the objects must have unique, bitwise object representations (e.g., no padding bits).
 * - Does not support erasing keys
 * - Capacity is fixed and will not grow automatically
 * - Requires the user to specify sentinel values for both key and mapped value
 * to indicate empty slots
 * - Concurrent insert/find is only supported when
 * <tt>static_multimap<Key, Value>::supports_concurrent_insert_find()</tt> is true`
 *
 * The `static_multimap` supports two types of operations:
 * - Host-side "bulk" operations
 * - Device-side "singular" operations
 *
 * The host-side bulk operations include `insert`, `contains`, `count`, `retrieve` and their
 * variants. These APIs should be used when there are a large number of keys to insert or lookup in
 * the map. For example, given a range of keys specified by device-accessible iterators, the bulk
 * `insert` function will insert all keys into the map.
 *
 * The singular device-side operations allow individual threads to perform
 * independent operations (e.g. `insert`, etc.) from device code. These
 * operations are accessed through non-owning, trivially copyable "view" types:
 * `device_view` and `device_mutable_view`. The `device_view` class is an
 * immutable view that allows only non-modifying operations such as `count` or
 * `contains`. The `device_mutable_view` class only allows `insert` operations.
 * The two types are separate to prevent erroneous concurrent insert/find
 * operations.
 *
 * By default, when querying for a Key `k` in operations like `count` or `retrieve`, if `k` is not
 * present in the map, it will not contribute to the output. Query APIs with the `_outer` suffix
 * will include non-matching keys in the output. See the relevant API documentation for more
 * information.
 *
 * Typical associative container query APIs like `retrieve` look up values by solely by key, e.g.,
 * `count` for a Key `k` will count all values whose associated key `k'` matches `k` as determined
 * by `key_equal(k, k')`. In some cases, one may want to consider both key _and_ value when
 * determining if a key-value pair should contribute to the output. `static_multimap` supports this
 * use case with APIs prefixed with `pair_`, e.g., `pair_count` is given a key-value pair
 * `{k,v}` and only counts key-value pairs, `{k', v'}`, in the map where `pair_equal({k,v}, {k',
 * v'})` is true. See the relevant API documentation for more information.
 *
 * Example:
 * \code{.cpp}
 * int empty_key_sentinel = -1;
 * int empty_value_sentinel = -1;
 *
 * // Constructs a multimap with 100,000 slots using -1 and -1 as the empty key/value
 * // sentinels. Note the capacity is chosen knowing we will insert 50,000 keys,
 * // for an load factor of 50%.
 * static_multimap<int, int> m{100'000, empty_key_sentinel, empty_value_sentinel};
 *
 * // Create a sequence of pairs {{0,0}, {1,1}, ... {i,i}}
 * thrust::device_vector<cuco::pair<int,int>> pairs(50,000);
 * thrust::transform(thrust::make_counting_iterator(0),
 *                   thrust::make_counting_iterator(pairs.size()),
 *                   pairs.begin(),
 *                   []__device__(auto i){ return cuco::pair{i,i}; };
 *
 * // Inserts all pairs into the map
 * m.insert(pairs.begin(), pairs.end());
 *
 * // Get a `device_view` and passes it to a kernel where threads may perform
 * // `contains/count/retrieve` lookups
 * kernel<<<...>>>(m.get_device_view());
 * \endcode
 *
 *
 * @tparam Key Type used for keys. Requires `cuco::is_bitwise_comparable_v<Key>`
 * @tparam Value Type of the mapped values. Requires `cuco::is_bitwise_comparable_v<Value>`
 * @tparam Scope The scope in which multimap operations will be performed by
 * individual threads
 * @tparam ProbeSequence Probe sequence chosen between `cuco::legacy::linear_probing`
 * and `cuco::legacy::double_hashing`. (see `probe_sequences.cuh`)
 * @tparam Allocator Type of allocator used for device storage
 */
template <typename Key,
          typename Value,
          cuda::thread_scope Scope = cuda::thread_scope_device,
          typename Allocator       = cuco::cuda_allocator<char>,
          class ProbeSequence = cuco::legacy::double_hashing<8, cuco::default_hash_function<Key>>>
class static_multimap {
  static_assert(
    cuco::is_bitwise_comparable_v<Key>,
    "Key type must have unique object representations or have been explicitly declared as safe for "
    "bitwise comparison via specialization of cuco::is_bitwise_comparable_v<Key>.");

  static_assert(
    cuco::is_bitwise_comparable_v<Value>,
    "Value type must have unique object representations or have been explicitly declared as safe "
    "for bitwise comparison via specialization of cuco::is_bitwise_comparable_v<Value>.");

  static_assert(std::is_base_of_v<cuco::legacy::detail::probe_sequence_base<ProbeSequence::cg_size>,
                                  ProbeSequence>,
                "ProbeSequence must be a specialization of either cuco::legacy::double_hashing or "
                "cuco::legacy::linear_probing.");

 public:
  using value_type         = cuco::pair<Key, Value>;            ///< Type of key/value pairs
  using key_type           = Key;                               ///< Key type
  using mapped_type        = Value;                             ///< Type of mapped values
  using size_type          = std::size_t;                       ///< Size type
  using atomic_key_type    = cuda::atomic<key_type, Scope>;     ///< Type of atomic keys
  using atomic_mapped_type = cuda::atomic<mapped_type, Scope>;  ///< Type of atomic mapped values
  using pair_atomic_type =
    cuco::pair<atomic_key_type,
               atomic_mapped_type>;  ///< Pair type of atomic key and atomic mapped value
  using allocator_type = typename std::allocator_traits<Allocator>::template rebind_alloc<
    pair_atomic_type>;  ///< Type of the allocator to (de)allocate slots
  using probe_sequence_type =
    cuco::legacy::detail::probe_sequence<ProbeSequence, Key, Value, Scope>;  ///< Probe scheme type

  static_multimap(static_multimap const&)            = delete;
  static_multimap& operator=(static_multimap const&) = delete;

  static_multimap(static_multimap&&) = default;  ///< Move constructor

  /**
   * @brief Replaces the contents of the map with another map.
   *
   * @return Reference of the current map object
   */
  static_multimap& operator=(static_multimap&&) = default;
  ~static_multimap()                            = default;

  /**
   * @brief Indicate if concurrent insert/find is supported for the key/value types.
   *
   * @return Boolean indicating if concurrent insert/find is supported.
   */
  __host__ __device__ __forceinline__ static constexpr bool
  supports_concurrent_insert_find() noexcept
  {
    return cuco::detail::is_packable<value_type>();
  }

  /**
   * @brief The size of the CUDA cooperative thread group.
   *
   * @return The CG size.
   */
  __host__ __device__ __forceinline__ static constexpr uint32_t cg_size() noexcept
  {
    return ProbeSequence::cg_size;
  }

  /**
   * @brief Construct a statically-sized map with the specified initial capacity,
   * sentinel values and CUDA stream.
   *
   * The capacity of the map is fixed. Insert operations will not automatically
   * grow the map. Attempting to insert more unique keys than the capacity of
   * the map results in undefined behavior.
   *
   * Performance begins to degrade significantly beyond a load factor of ~70%.
   * For best performance, choose a capacity that will keep the load factor
   * below 70%. E.g., if inserting `N` unique keys, choose a capacity of
   * `N * (1/0.7)`.
   *
   * The `empty_key_sentinel` and `empty_value_sentinel` values are reserved and
   * undefined behavior results from attempting to insert any key/value pair
   * that contains either.
   *
   * @param capacity The total number of slots in the map
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param empty_value_sentinel The reserved mapped value for empty slots
   * @param stream CUDA stream used to initialize the map
   * @param alloc Allocator used for allocating device storage
   */
  static_multimap(std::size_t capacity,
                  empty_key<Key> empty_key_sentinel,
                  empty_value<Value> empty_value_sentinel,
                  cudaStream_t stream    = 0,
                  Allocator const& alloc = Allocator{});

  /**
   * @brief Inserts all key/value pairs in the range `[first, last)`.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multimap<K, V>::value_type></tt> is `true`
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt>
  void insert(InputIt first, InputIt last, cudaStream_t stream = 0);

  /**
   * @brief Inserts key/value pairs in the range `[first, first + n)` if `pred`
   * of the corresponding stencil returns true.
   *
   * The key/value pair `*(first + i)` is inserted if `pred( *(stencil + i) )` returns true.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multimap<K, V>::value_type></tt> is `true`
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>.
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt, typename StencilIt, typename Predicate>
  void insert_if(
    InputIt first, InputIt last, StencilIt stencil, Predicate pred, cudaStream_t stream = 0);

  /**
   * @brief Indicates whether the keys in the range `[first, last)` are contained in the map.
   *
   * Stores `true` or `false` to `(output + i)` indicating if the key `*(first + i)` exists in the
   * map.
   *
   * ProbeSequence hashers should be callable with both
   * <tt>std::iterator_traits<InputIt>::value_type</tt> and Key type.
   * <tt>std::invoke_result<KeyEqual, std::iterator_traits<InputIt>::value_type, Key></tt> must be
   * well-formed.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   * @tparam KeyEqual Binary callable type used to compare two keys for equality
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the output sequence indicating whether each key is present
   * @param key_equal The binary function to compare two keys for equality
   * @param stream CUDA stream used for contains
   */
  template <typename InputIt, typename OutputIt, typename KeyEqual = cuda::std::equal_to<key_type>>
  void contains(InputIt first,
                InputIt last,
                OutputIt output_begin,
                KeyEqual key_equal  = KeyEqual{},
                cudaStream_t stream = 0) const;

  /**
   * @brief Indicates whether the pairs in the range `[first, last)` are contained in the map.
   *
   * Stores `true` or `false` to `(output + i)` indicating if the pair `*(first + i)` exists in
   * the map.
   *
   * ProbeSequence hashers should be callable with both
   * <tt>std::iterator_traits<InputIt>::value_type::first_type</tt>
   * and Key type. <tt>std::invoke_result<KeyEqual,
   * std::iterator_traits<InputIt>::value_type::first_type, Key></tt>
   * must be well-formed.
   *
   * @tparam InputIt Device accessible random access input iterator
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   * @tparam PairEqual Binary callable type used to compare input pair and slot content for equality
   *
   * @param first Beginning of the sequence of pairs
   * @param last End of the sequence of pairs
   * @param output_begin Beginning of the output sequence indicating whether each pair is present
   * @param pair_equal The binary function to compare input pair and slot content for equality
   * @param stream CUDA stream used for contains
   */
  template <typename InputIt, typename OutputIt, typename PairEqual>
  void pair_contains(InputIt first,
                     InputIt last,
                     OutputIt output_begin,
                     PairEqual pair_equal,
                     cudaStream_t stream = 0) const;

  /**
   * @brief Counts the occurrences of keys in `[first, last)` contained in the multimap.
   *
   * For each key, `k = *(first + i)`, counts all matching keys, `k'`, as determined by
   * `key_equal(k, k')` and returns the sum of all matches for all keys.
   *
   * @tparam Input Device accessible input iterator whose `value_type` is convertible to `key_type`
   * @tparam KeyEqual Binary callable
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param stream CUDA stream used for count
   * @param key_equal Binary function to compare two keys for equality
   * @return The sum of total occurrences of all keys in `[first, last)`
   */
  template <typename InputIt, typename KeyEqual = cuda::std::equal_to<key_type>>
  std::size_t count(InputIt first,
                    InputIt last,
                    cudaStream_t stream = 0,
                    KeyEqual key_equal  = KeyEqual{}) const;

  /**
   * @brief Counts the occurrences of keys in `[first, last)` contained in the multimap.
   *
   * For each key, `k = *(first + i)`, counts all matching keys, `k'`, as determined by
   * `key_equal(k, k')` and returns the sum of all matches for all keys. If `k` does not have any
   * matches, it contributes 1 to the final sum.
   *
   * @tparam Input Device accessible input iterator whose `value_type` is convertible to `key_type`
   * @tparam KeyEqual Binary callable
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param stream CUDA stream used for count_outer
   * @param key_equal Binary function to compare two keys for equality
   * @return The sum of total occurrences of all keys in `[first, last)` where keys without matches
   * are considered to have a single occurrence.
   */
  template <typename InputIt, typename KeyEqual = cuda::std::equal_to<key_type>>
  std::size_t count_outer(InputIt first,
                          InputIt last,
                          cudaStream_t stream = 0,
                          KeyEqual key_equal  = KeyEqual{}) const;

  /**
   * @brief Counts the occurrences of key/value pairs in `[first, last)` contained in the multimap.
   *
   * For key-value pair, `kv = *(first + i)`, counts all matching key-value pairs, `kv'`, as
   * determined by `pair_equal(kv, kv')` and returns the sum of all matches for all key-value pairs.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multimap<K, V>::value_type></tt> is `true`
   * @tparam PairEqual Binary callable
   * @param first Beginning of the sequence of pairs to count
   * @param last End of the sequence of pairs to count
   * @param pair_equal Binary function to compare two pairs for equality
   * @param stream CUDA stream used for pair_count
   * @return The sum of total occurrences of all pairs in `[first, last)`
   */
  template <typename InputIt, typename PairEqual>
  std::size_t pair_count(InputIt first,
                         InputIt last,
                         PairEqual pair_equal,
                         cudaStream_t stream = 0) const;

  /**
   * @brief Counts the occurrences of key/value pairs in `[first, last)` contained in the multimap.
   *
   * For key-value pair, `kv = *(first + i)`, counts all matching key-value pairs, `kv'`, as
   * determined by `pair_equal(kv, kv')` and returns the sum of all matches for all key-value pairs.
   * if `kv` does not have any matches, it contributes 1 to the final sum.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multimap<K, V>::value_type></tt> is `true`
   * @tparam PairEqual Binary callable
   * @param first Beginning of the sequence of pairs to count
   * @param last End of the sequence of pairs to count
   * @param pair_equal Binary function to compare two pairs for equality
   * @param stream CUDA stream used for pair_count_outer
   * @return The sum of total occurrences of all pairs in `[first, last)` where a key-value pair
   * without a match is considered to have a single occurrence
   */
  template <typename InputIt, typename PairEqual>
  std::size_t pair_count_outer(InputIt first,
                               InputIt last,
                               PairEqual pair_equal,
                               cudaStream_t stream = 0) const;

  /**
   * @brief Retrieves all the values corresponding to all keys in the range `[first, last)`.
   *
   * If key `k = *(first + i)` exists in the map, copies `k` and all associated values to
   * unspecified locations in `[output_begin, output_end)`. Else, does nothing.
   *
   * Behavior is undefined if the size of the output range exceeds `std::distance(output_begin,
   * output_end)`. Use `count()` to determine the size of the output range.
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `key_type`
   * @tparam OutputIt Device accessible output iterator whose `value_type` is
   * constructible from the map's `value_type`
   * @tparam KeyEqual Binary callable type
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of key/value pairs retrieved for each key
   * @param stream CUDA stream used for retrieve
   * @param key_equal The binary function to compare two keys for equality
   * @return The iterator indicating the last valid key/value pairs in the output
   */
  template <typename InputIt, typename OutputIt, typename KeyEqual = cuda::std::equal_to<key_type>>
  OutputIt retrieve(InputIt first,
                    InputIt last,
                    OutputIt output_begin,
                    cudaStream_t stream = 0,
                    KeyEqual key_equal  = KeyEqual{}) const;

  /**
   * @brief Retrieves all the matches corresponding to all keys in the range `[first, last)`.
   *
   * If key `k = *(first + i)` exists in the map, copies `k` and all associated values to
   * unspecified locations in `[output_begin, output_end)`. Else, copies `k` and
   * `empty_value_sentinel`.
   *
   * Behavior is undefined if the size of the output range exceeds `std::distance(output_begin,
   * output_end)`. Use `count_outer()` to determine the size of the output range.
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `key_type`
   * @tparam OutputIt Device accessible output iterator whose `value_type` is
   * constructible from the map's `value_type`
   * @tparam KeyEqual Binary callable type
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of key/value pairs retrieved for each key
   * @param stream CUDA stream used for retrieve_outer
   * @param key_equal The binary function to compare two keys for equality
   * @return The iterator indicating the last valid key/value pairs in the output
   */
  template <typename InputIt, typename OutputIt, typename KeyEqual = cuda::std::equal_to<key_type>>
  OutputIt retrieve_outer(InputIt first,
                          InputIt last,
                          OutputIt output_begin,
                          cudaStream_t stream = 0,
                          KeyEqual key_equal  = KeyEqual{}) const;

  /**
   * @brief Retrieves all pairs matching the input probe pair in the range `[first, last)`.
   *
   * The `pair_` prefix indicates that the input data type is convertible to the map's
   * `value_type`. If pair_equal(*(first + i), slot[j]) returns true, then *(first+i) is
   * stored to `probe_output_begin`, and slot[j] is stored to `contained_output_begin`.
   *
   * Behavior is undefined if the size of the output range exceeds
   * `std::distance(probe_output_begin, probe_output_end)` (or
   * `std::distance(contained_output_begin, contained_output_end)`). Use
   * `pair_count()` to determine the size of the output range.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multimap<K, V>::value_type></tt> is `true`
   * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
   * `InputIt`s `value_type`.
   * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
   * the map's `value_type`.
   * @tparam PairEqual Binary callable type
   * @param first Beginning of the sequence of pairs
   * @param last End of the sequence of pairs
   * @param probe_output_begin Beginning of the sequence of the matched probe pairs
   * @param contained_output_begin Beginning of the sequence of the matched contained pairs
   * @param pair_equal The binary function to compare two pairs for equality
   * @param stream CUDA stream used for pair_retrieve
   * @return Pair of iterators pointing to the last elements in the output
   */
  template <typename InputIt, typename OutputIt1, typename OutputIt2, typename PairEqual>
  std::pair<OutputIt1, OutputIt2> pair_retrieve(InputIt first,
                                                InputIt last,
                                                OutputIt1 probe_output_begin,
                                                OutputIt2 contained_output_begin,
                                                PairEqual pair_equal,
                                                cudaStream_t stream = 0) const;

  /**
   * @brief Retrieves all pairs matching the input probe pair in the range `[first, last)`.
   *
   * The `pair_` prefix indicates that the input data type is convertible to the map's `value_type`.
   * If pair_equal(*(first + i), slot[j]) returns true, then *(first+i) is stored to
   * `probe_output_begin`, and slot[j] is stored to `contained_output_begin`. If *(first+i) doesn't
   * have matches in the map, copies *(first + i) in `probe_output_begin` and a pair of
   * `empty_key_sentinel` and `empty_value_sentinel` in `contained_output_begin`.
   *
   * Behavior is undefined if the size of the output range exceeds
   * `std::distance(probe_output_begin, probe_output_end)` (or
   * `std::distance(contained_output_begin, contained_output_end)`). Use
   * `pair_count()` to determine the size of the output range.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_multimap<K, V>::value_type></tt> is `true`
   * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
   * `InputIt`s `value_type`.
   * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
   * the map's `value_type`.
   * @tparam PairEqual Binary callable type
   * @param first Beginning of the sequence of pairs
   * @param last End of the sequence of pairs
   * @param probe_output_begin Beginning of the sequence of the matched probe pairs
   * @param contained_output_begin Beginning of the sequence of the matched contained pairs
   * @param pair_equal The binary function to compare two pairs for equality
   * @param stream CUDA stream used for pair_retrieve_outer
   * @return Pair of iterators pointing to the last elements in the output
   */
  template <typename InputIt, typename OutputIt1, typename OutputIt2, typename PairEqual>
  std::pair<OutputIt1, OutputIt2> pair_retrieve_outer(InputIt first,
                                                      InputIt last,
                                                      OutputIt1 probe_output_begin,
                                                      OutputIt2 contained_output_begin,
                                                      PairEqual pair_equal,
                                                      cudaStream_t stream = 0) const;

 private:
  /**
   * @brief Indicates if vector-load is used.
   *
   * Users have no explicit control on whether vector-load is used.
   *
   * @return Boolean indicating if vector-load is used.
   */
  static __host__ __device__ constexpr bool uses_vector_load() noexcept
  {
    return cuco::detail::is_packable<value_type>();
  }

  /**
   * @brief Returns the number of pairs loaded with each vector-load
   */
  static __host__ __device__ constexpr uint32_t vector_width() noexcept
  {
    return ProbeSequence::vector_width();
  }

  /**
   * @brief Returns the warp size.
   */
  static __host__ __device__ constexpr uint32_t warp_size() noexcept { return 32u; }

  /**
   * @brief Custom deleter for unique pointer of slots.
   */
  struct slot_deleter {
    slot_deleter(allocator_type& a, size_t& c, cuda::stream_ref s)
      : allocator{a}, capacity{c}, stream{s}
    {
    }

    slot_deleter(slot_deleter const&) = default;

    void operator()(pair_atomic_type* ptr) { allocator.deallocate(ptr, capacity, stream); }

    allocator_type& allocator;
    size_t& capacity;
    cuda::stream_ref stream;
  };

  class device_view_impl_base;
  class device_mutable_view_impl;
  class device_view_impl;

  template <typename ViewImpl>
  class device_view_base {
   protected:
    // Import member type definitions from `static_multimap`
    using value_type          = value_type;
    using key_type            = Key;
    using mapped_type         = Value;
    using pair_atomic_type    = pair_atomic_type;
    using iterator            = pair_atomic_type*;
    using const_iterator      = pair_atomic_type const*;
    using probe_sequence_type = probe_sequence_type;

    __host__ __device__ device_view_base(pair_atomic_type* slots,
                                         std::size_t capacity,
                                         empty_key<Key> empty_key_sentinel,
                                         empty_value<Value> empty_value_sentinel) noexcept
      : impl_{slots, capacity, empty_key_sentinel.value, empty_value_sentinel.value}
    {
    }

   public:
    /**
     * @brief Gets slots array.
     *
     * @return Slots array
     */
    __device__ __forceinline__ pair_atomic_type* get_slots() noexcept { return impl_.get_slots(); }

    /**
     * @brief Gets slots array.
     *
     * @return Slots array
     */
    __device__ __forceinline__ pair_atomic_type const* get_slots() const noexcept
    {
      return impl_.get_slots();
    }

    /**
     * @brief Gets the maximum number of elements the hash map can hold.
     *
     * @return The maximum number of elements the hash map can hold
     */
    __host__ __device__ __forceinline__ std::size_t get_capacity() const noexcept
    {
      return impl_.get_capacity();
    }

    /**
     * @brief Gets the sentinel value used to represent an empty key slot.
     *
     * @return The sentinel value used to represent an empty key slot
     */
    __host__ __device__ __forceinline__ Key get_empty_key_sentinel() const noexcept
    {
      return impl_.get_empty_key_sentinel();
    }

    /**
     * @brief Gets the sentinel value used to represent an empty value slot.
     *
     * @return The sentinel value used to represent an empty value slot
     */
    __host__ __device__ __forceinline__ Value get_empty_value_sentinel() const noexcept
    {
      return impl_.get_empty_value_sentinel();
    }

   protected:
    ViewImpl impl_;
  };  // class device_view_base

 public:
  /**
   * @brief Mutable, non-owning view-type that may be used in device code to
   * perform singular inserts into the map.
   *
   * `device_mutable_view` is trivially-copyable and is intended to be passed by
   * value.
   *
   * Example:
   * \code{.cpp}
   * cuco::static_multimap<int,int> m{100'000, -1, -1};
   *
   * // Inserts a sequence of pairs {{0,0}, {1,1}, ... {i,i}}
   * thrust::for_each(thrust::make_counting_iterator(0),
   *                  thrust::make_counting_iterator(50'000),
   *                  [map = m.get_device_mutable_view()]
   *                  __device__ (auto i) mutable {
   *                     map.insert(cuco::pair{i,i});
   *                  });
   * \endcode
   */
  class device_mutable_view : public device_view_base<device_mutable_view_impl> {
   public:
    using view_base_type =
      device_view_base<device_mutable_view_impl>;              ///< Base view implementation type
    using value_type  = typename view_base_type::value_type;   ///< Type of key/value pairs
    using key_type    = typename view_base_type::key_type;     ///< Key type
    using mapped_type = typename view_base_type::mapped_type;  ///< Type of the mapped values
    using iterator =
      typename view_base_type::iterator;  ///< Type of the forward iterator to `value_type`
    using const_iterator =
      typename view_base_type::const_iterator;  ///< Type of the forward iterator to `const
                                                ///< value_type`

    /**
     * @brief Construct a mutable view of the first `capacity` slots of the
     * slots array pointed to by `slots`.
     *
     * @param slots Pointer to beginning of initialized slots array
     * @param capacity The number of slots viewed by this object
     * @param empty_key_sentinel The reserved value for keys to represent empty
     * slots
     * @param empty_value_sentinel The reserved value for mapped values to
     * represent empty slots
     */
    __host__ __device__ device_mutable_view(pair_atomic_type* slots,
                                            std::size_t capacity,
                                            empty_key<Key> empty_key_sentinel,
                                            empty_value<Value> empty_value_sentinel) noexcept
      : view_base_type{slots, capacity, empty_key_sentinel, empty_value_sentinel}
    {
    }

    /**
     * @brief Inserts the specified key/value pair into the map.
     *
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param g The Cooperative Group that performs the insert
     * @param insert_pair The pair to insert
     */
    template <typename ParentCG>
    __device__ __forceinline__ void insert(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
      value_type const& insert_pair) noexcept;

   private:
    using device_view_base<device_mutable_view_impl>::impl_;
  };  // class device mutable view

  /**
   * @brief Non-owning view-type that may be used in device code to
   * perform singular find and contains operations for the map.
   *
   * `device_view` is trivially-copyable and is intended to be passed by
   * value.
   *
   */
  class device_view : public device_view_base<device_view_impl> {
   public:
    using view_base_type = device_view_base<device_view_impl>;    ///< Base view implementation type
    using value_type     = typename view_base_type::value_type;   ///< Type of key/value pairs
    using key_type       = typename view_base_type::key_type;     ///< Key type
    using mapped_type    = typename view_base_type::mapped_type;  ///< Type of the mapped values
    using iterator =
      typename view_base_type::iterator;  ///< Type of the forward iterator to `value_type`
    using const_iterator =
      typename view_base_type::const_iterator;  ///< Type of the forward iterator to `const
                                                ///< value_type`

    /**
     * @brief Construct a view of the first `capacity` slots of the
     * slots array pointed to by `slots`.
     *
     * @param slots Pointer to beginning of initialized slots array
     * @param capacity The number of slots viewed by this object
     * @param empty_key_sentinel The reserved value for keys to represent empty
     * slots
     * @param empty_value_sentinel The reserved value for mapped values to
     * represent empty slots
     */
    __host__ __device__ device_view(pair_atomic_type* slots,
                                    std::size_t capacity,
                                    empty_key<Key> empty_key_sentinel,
                                    empty_value<Value> empty_value_sentinel) noexcept
      : view_base_type{slots, capacity, empty_key_sentinel, empty_value_sentinel}
    {
    }

    /**
     * @brief Makes a copy of given `device_view` using non-owned memory.
     *
     * This function is intended to be used to create shared memory copies of small static maps,
     * although global memory can be used as well.
     *
     * @tparam CG The type of the cooperative thread group
     * @param g The cooperative thread group used to copy the slots
     * @param source_device_view `device_view` to copy from
     * @param memory_to_use Array large enough to support `capacity` elements. Object does not
     * take the ownership of the memory
     * @return Copy of passed `device_view`
     */
    template <typename CG>
    __device__ __forceinline__ static device_view make_copy(
      CG g, pair_atomic_type* const memory_to_use, device_view source_device_view) noexcept;

    /**
     * @brief Flushes per-CG buffer into the output sequence.
     *
     * A given CUDA Cooperative Group, `g`, loads `num_outputs` key-value pairs from `output_buffer`
     * and writes them into global memory in a coalesced fashion. CG-wide `memcpy_sync` is used if
     * `thrust::is_contiguous_iterator_v<OutputIt>` returns true. All threads of `g` must be active
     * due to implicit CG-wide synchronization during flushing.
     *
     * @tparam CG Cooperative Group type
     * @tparam atomicT Type of atomic storage
     * @tparam OutputIt Device accessible output iterator whose `value_type` is
     * constructible from the map's `value_type`
     * @param g The Cooperative Group used to flush output buffer
     * @param num_outputs Number of valid output in the buffer
     * @param output_buffer Buffer of the key/value pair sequence
     * @param num_matches Size of the output sequence
     * @param output_begin Beginning of the output sequence of key/value pairs
     */
    template <typename CG, typename atomicT, typename OutputIt>
    __device__ __forceinline__ void flush_output_buffer(CG g,
                                                        uint32_t const num_outputs,
                                                        value_type* output_buffer,
                                                        atomicT* num_matches,
                                                        OutputIt output_begin) noexcept;

    /**
     * @brief Flushes per-CG buffer into the output sequences.
     *
     * A given CUDA Cooperative Group, `g`, loads `num_outputs` elements from `probe_output_buffer`
     * and `num_outputs` elements from `contained_output_buffer`, then writes them into global
     * memory started from `probe_output_begin` and `contained_output_begin` respectively. All
     * threads of `g` must be active due to implicit CG-wide synchronization during flushing.
     *
     * @tparam CG Cooperative Group type
     * @tparam atomicT Type of atomic storage
     * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
     * `InputIt`s `value_type`.
     * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
     * the map's `value_type`.
     * @param g The Cooperative Group used to flush output buffer
     * @param num_outputs Number of valid output in the buffer
     * @param probe_output_buffer Buffer of the matched probe pair sequence
     * @param contained_output_buffer Buffer of the matched contained pair sequence
     * @param num_matches Size of the output sequence
     * @param probe_output_begin Beginning of the output sequence of the matched probe pairs
     * @param contained_output_begin Beginning of the output sequence of the matched contained
     * pairs
     */
    template <typename CG, typename atomicT, typename OutputIt1, typename OutputIt2>
    __device__ __forceinline__ void flush_output_buffer(CG g,
                                                        uint32_t const num_outputs,
                                                        value_type* probe_output_buffer,
                                                        value_type* contained_output_buffer,
                                                        atomicT* num_matches,
                                                        OutputIt1 probe_output_begin,
                                                        OutputIt2 contained_output_begin) noexcept;

    /**
     * @brief Indicates whether the key `k` exists in the map.
     *
     * If the key `k` was inserted into the map, `contains` returns
     * true. Otherwise, it returns false. Uses the CUDA Cooperative Groups API to
     * to leverage multiple threads to perform a single `contains` operation. This provides a
     * significant boost in throughput compared to the non Cooperative Group
     * `contains` at moderate to high load factors.
     *
     * ProbeSequence hashers should be callable with both ProbeKey and Key type.
     * `std::invoke_result<KeyEqual, ProbeKey, Key>` must be well-formed.
     *
     * If `key_equal(probe_key, slot_key)` returns true, `hash(probe_key) == hash(slot_key)` must
     * also be true.
     *
     * @tparam ProbeKey Probe key type
     * @tparam KeyEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param g The Cooperative Group used to perform the contains operation
     * @param k The key to search for
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return A boolean indicating whether the key/value pair
     * containing `k` was inserted
     */
    template <typename ProbeKey,
              typename KeyEqual = cuda::std::equal_to<key_type>,
              typename ParentCG = void>
    __device__ __forceinline__ bool contains(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
      ProbeKey const& k,
      KeyEqual key_equal = KeyEqual{}) const noexcept;

    /**
     * @brief Indicates whether the pair `p` exists in the map.
     *
     * If the pair `p` was inserted into the map, `contains` returns
     * true. Otherwise, it returns false. Uses the CUDA Cooperative Groups API to
     * to leverage multiple threads to perform a single `contains` operation. This provides a
     * significant boost in throughput compared to the non Cooperative Group
     * `contains` at moderate to high load factors.
     *
     * ProbeSequence hashers should be callable with both ProbePair::first_type and Key type.
     * `std::invoke_result<KeyEqual, ProbePair::first_type, Key>` must be well-formed.
     *
     * If `pair_equal(p, slot_content)` returns true, `hash(p.first) == hash(slot_key)` must
     * also be true.
     *
     * @tparam ProbePair Probe pair type
     * @tparam PairEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param g The Cooperative Group used to perform the contains operation
     * @param p The pair to search for
     * @param pair_equal The binary callable used to compare input pair and slot content
     * for equality
     * @return A boolean indicating whether the input pair was inserted in the map
     */
    template <typename ProbePair, typename PairEqual, typename ParentCG>
    __device__ __forceinline__ bool pair_contains(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
      ProbePair const& p,
      PairEqual pair_equal) const noexcept;

    /**
     * @brief Counts the occurrence of a given key contained in multimap.
     *
     * For a given key, `k`, counts all matching keys, `k'`, as determined by `key_equal(k, k')` and
     * returns the sum of all matches for `k`.
     *
     * @tparam KeyEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param g The Cooperative Group used to perform the count operation
     * @param k The key to search for
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return Number of matches found by the current thread
     */
    template <typename KeyEqual = cuda::std::equal_to<key_type>, typename ParentCG = void>
    __device__ __forceinline__ std::size_t count(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
      Key const& k,
      KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Counts the occurrence of a given key contained in multimap. If no
     * matches can be found for a given key, the corresponding occurrence is 1.
     *
     * For a given key, `k`, counts all matching keys, `k'`, as determined by `key_equal(k, k')` and
     * returns the sum of all matches for `k`. If `k` does not have any matches, returns 1.
     *
     * @tparam KeyEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param g The Cooperative Group used to perform the count operation
     * @param k The key to search for
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return Number of matches found by the current thread
     */
    template <typename KeyEqual = cuda::std::equal_to<key_type>, typename ParentCG = void>
    __device__ __forceinline__ std::size_t count_outer(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
      Key const& k,
      KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Counts the occurrence of a given key/value pair contained in multimap.
     *
     * For a given pair, `p`, counts all matching pairs, `p'`, as determined by `pair_equal(p, p')`
     * and returns the sum of all matches for `p`.
     *
     * @tparam PairEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param g The Cooperative Group used to perform the pair_count operation
     * @param pair The pair to search for
     * @param pair_equal The binary callable used to compare two pairs
     * for equality
     * @return Number of matches found by the current thread
     */
    template <typename PairEqual, typename ParentCG>
    __device__ __forceinline__ std::size_t pair_count(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
      value_type const& pair,
      PairEqual pair_equal) noexcept;

    /**
     * @brief Counts the occurrence of a given key/value pair contained in multimap.
     * If no matches can be found for a given key, the corresponding occurrence is 1.
     *
     * For a given pair, `p`, counts all matching pairs, `p'`, as determined by `pair_equal(p, p')`
     * and returns the sum of all matches for `p`. If `p` does not have any matches, returns 1.
     *
     * @tparam PairEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param g The Cooperative Group used to perform the pair_count operation
     * @param pair The pair to search for
     * @param pair_equal The binary callable used to compare two pairs
     * for equality
     * @return Number of matches found by the current thread
     */
    template <typename PairEqual, typename ParentCG>
    __device__ __forceinline__ std::size_t pair_count_outer(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> g,
      value_type const& pair,
      PairEqual pair_equal) noexcept;

    /**
     * @brief Retrieves all the matches of a given key contained in multimap with per-flushing-CG
     * shared memory buffer.
     *
     * For key `k` existing in the map, copies `k` and all associated values to unspecified
     * locations in `[output_begin, output_end)`.
     *
     * @tparam buffer_size Size of the output buffer
     * @tparam FlushingCG Type of Cooperative Group used to flush output buffer
     * @tparam atomicT Type of atomic storage
     * @tparam OutputIt Device accessible output iterator whose `value_type` is
     * constructible from the map's `value_type`
     * @tparam KeyEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param flushing_cg The Cooperative Group used to flush output buffer
     * @param probing_cg The Cooperative Group used to retrieve
     * @param k The key to search for
     * @param flushing_cg_counter Pointer to flushing_cg counter
     * @param output_buffer Shared memory buffer of the key/value pair sequence
     * @param num_matches Size of the output sequence
     * @param output_begin Beginning of the output sequence of key/value pairs
     * @param key_equal The binary callable used to compare two keys
     * for equality
     */
    template <uint32_t buffer_size,
              typename FlushingCG,
              typename atomicT,
              typename OutputIt,
              typename KeyEqual = cuda::std::equal_to<key_type>,
              typename ParentCG = void>
    __device__ __forceinline__ void retrieve(
      FlushingCG flushing_cg,
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
      Key const& k,
      uint32_t* flushing_cg_counter,
      value_type* output_buffer,
      atomicT* num_matches,
      OutputIt output_begin,
      KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Retrieves all the matches of a given key contained in multimap with per-flushing-CG
     * shared memory buffer.
     *
     * For key `k` existing in the map, copies `k` and all associated values to unspecified
     * locations in `[output_begin, output_end)`. If `k` does not have any matches, copies `k` and
     * `empty_value_sentinel()` into the output.
     *
     * @tparam buffer_size Size of the output buffer
     * @tparam FlushingCG Type of Cooperative Group used to flush output buffer
     * @tparam atomicT Type of atomic storage
     * @tparam OutputIt Device accessible output iterator whose `value_type` is
     * constructible from the map's `value_type`
     * @tparam KeyEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param flushing_cg The Cooperative Group used to flush output buffer
     * @param probing_cg The Cooperative Group used to retrieve
     * @param k The key to search for
     * @param flushing_cg_counter Pointer to flushing_cg counter
     * @param output_buffer Shared memory buffer of the key/value pair sequence
     * @param num_matches Size of the output sequence
     * @param output_begin Beginning of the output sequence of key/value pairs
     * @param key_equal The binary callable used to compare two keys
     * for equality
     */
    template <uint32_t buffer_size,
              typename FlushingCG,
              typename atomicT,
              typename OutputIt,
              typename KeyEqual = cuda::std::equal_to<key_type>,
              typename ParentCG = void>
    __device__ __forceinline__ void retrieve_outer(
      FlushingCG flushing_cg,
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
      Key const& k,
      uint32_t* flushing_cg_counter,
      value_type* output_buffer,
      atomicT* num_matches,
      OutputIt output_begin,
      KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Retrieves all the matches of a given pair
     *
     * For pair `p` with `n = pair_count(cg, p, pair_equal)` matching pairs, if `pair_equal(p,
     * slot)` returns true, stores `probe_key_begin[j] = p.first`, `probe_val_begin[j] = p.second`,
     * `contained_key_begin[j] = slot.first`, and `contained_val_begin[j] = slot.second` for an
     * unspecified value of `j` where `0 <= j < n`.
     *
     * Concurrent reads or writes to any of the output ranges results in undefined behavior.
     *
     * Behavior is undefined if the extent of any of the output ranges is less than `n`.
     *
     * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
     * `pair`'s `Key` type.
     * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
     * `pair`'s `Value` type.
     * @tparam OutputIt3 Device accessible output iterator whose `value_type` is constructible from
     * the map's `key_type`.
     * @tparam OutputIt4 Device accessible output iterator whose `value_type` is constructible from
     * the map's `mapped_type`.
     * @tparam PairEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param probing_cg The Cooperative Group used to retrieve
     * @param pair The pair to search for
     * @param probe_key_begin Beginning of the output sequence of the matched probe keys
     * @param probe_val_begin Beginning of the output sequence of the matched probe values
     * @param contained_key_begin Beginning of the output sequence of the matched contained keys
     * @param contained_val_begin Beginning of the output sequence of the matched contained values
     * @param pair_equal The binary callable used to compare two pairs for equality
     */
    template <typename OutputIt1,
              typename OutputIt2,
              typename OutputIt3,
              typename OutputIt4,
              typename PairEqual,
              typename ParentCG>
    __device__ __forceinline__ void pair_retrieve(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
      value_type const& pair,
      OutputIt1 probe_key_begin,
      OutputIt2 probe_val_begin,
      OutputIt3 contained_key_begin,
      OutputIt4 contained_val_begin,
      PairEqual pair_equal) noexcept;

    /**
     * @brief Retrieves all the matches of a given pair contained in multimap with per-flushing-CG
     * shared memory buffer.
     *
     * For pair `p`, if pair_equal(p, slot[j]) returns true, copies `p` to unspecified locations
     * in `[probe_output_begin, probe_output_end)` and copies slot[j] to unspecified locations in
     * `[contained_output_begin, contained_output_end)`.
     *
     * @tparam buffer_size Size of the output buffer
     * @tparam FlushingCG Type of Cooperative Group used to flush output buffer
     * @tparam atomicT Type of atomic storage
     * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
     * `InputIt`s `value_type`.
     * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
     * the map's `value_type`.
     * @tparam PairEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param flushing_cg The Cooperative Group used to flush output buffer
     * @param probing_cg The Cooperative Group used to retrieve
     * @param pair The pair to search for
     * @param warp_counter Pointer to the warp counter
     * @param probe_output_buffer Buffer of the matched probe pair sequence
     * @param contained_output_buffer Buffer of the matched contained pair sequence
     * @param num_matches Size of the output sequence
     * @param probe_output_begin Beginning of the output sequence of the matched probe pairs
     * @param contained_output_begin Beginning of the output sequence of the matched contained
     * pairs
     * @param pair_equal The binary callable used to compare two pairs for equality
     */
    template <uint32_t buffer_size,
              typename FlushingCG,
              typename atomicT,
              typename OutputIt1,
              typename OutputIt2,
              typename PairEqual,
              typename ParentCG>
    __device__ __forceinline__ void pair_retrieve(
      FlushingCG flushing_cg,
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
      value_type const& pair,
      uint32_t* warp_counter,
      value_type* probe_output_buffer,
      value_type* contained_output_buffer,
      atomicT* num_matches,
      OutputIt1 probe_output_begin,
      OutputIt2 contained_output_begin,
      PairEqual pair_equal) noexcept;

    /**
     * @brief Retrieves all the matches of a given pair
     *
     * For pair `p` with `n = pair_count_outer(cg, p, pair_equal)` matching pairs, if `pair_equal(p,
     * slot)` returns true, stores `probe_key_begin[j] = p.first`, `probe_val_begin[j] = p.second`,
     * `contained_key_begin[j] = slot.first`, and `contained_val_begin[j] = slot.second` for an
     * unspecified value of `j` where `0 <= j < n`. If `p` does not have any matches, stores
     * `probe_key_begin[0] = p.first`, `probe_val_begin[0] = p.second`, `contained_key_begin[0] =
     * empty_key_sentinel`, and `contained_val_begin[0] = empty_value_sentinel`.
     *
     * Concurrent reads or writes to any of the output ranges results in undefined behavior.
     *
     * Behavior is undefined if the extent of any of the output ranges is less than `n`.
     *
     * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
     * `pair`'s `Key` type.
     * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
     * `pair`'s `Value` type.
     * @tparam OutputIt3 Device accessible output iterator whose `value_type` is constructible from
     * the map's `key_type`.
     * @tparam OutputIt4 Device accessible output iterator whose `value_type` is constructible from
     * the map's `mapped_type`.
     * @tparam PairEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param probing_cg The Cooperative Group used to retrieve
     * @param pair The pair to search for
     * @param probe_key_begin Beginning of the output sequence of the matched probe keys
     * @param probe_val_begin Beginning of the output sequence of the matched probe values
     * @param contained_key_begin Beginning of the output sequence of the matched contained keys
     * @param contained_val_begin Beginning of the output sequence of the matched contained values
     * @param pair_equal The binary callable used to compare two pairs for equality
     */
    template <typename OutputIt1,
              typename OutputIt2,
              typename OutputIt3,
              typename OutputIt4,
              typename PairEqual,
              typename ParentCG>
    __device__ __forceinline__ void pair_retrieve_outer(
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
      value_type const& pair,
      OutputIt1 probe_key_begin,
      OutputIt2 probe_val_begin,
      OutputIt3 contained_key_begin,
      OutputIt4 contained_val_begin,
      PairEqual pair_equal) noexcept;

    /**
     * @brief Retrieves all the matches of a given pair contained in multimap with per-flushing-CG
     * shared memory buffer.
     *
     * For pair `p`, if pair_equal(p, slot[j]) returns true, copies `p` to unspecified locations
     * in `[probe_output_begin, probe_output_end)` and copies slot[j] to unspecified locations in
     * `[contained_output_begin, contained_output_end)`. If `p` does not have any matches, copies
     * `p` and a pair of `empty_key_sentinel` and `empty_value_sentinel` into the output.
     *
     * @tparam buffer_size Size of the output buffer
     * @tparam FlushingCG Type of Cooperative Group used to flush output buffer
     * @tparam atomicT Type of atomic storage
     * @tparam OutputIt1 Device accessible output iterator whose `value_type` is constructible from
     * `InputIt`s `value_type`.
     * @tparam OutputIt2 Device accessible output iterator whose `value_type` is constructible from
     * the map's `value_type`.
     * @tparam PairEqual Binary callable type
     * @tparam ParentCG Type of parent Cooperative Group
     *
     * @param flushing_cg The Cooperative Group used to flush output buffer
     * @param probing_cg The Cooperative Group used to retrieve
     * @param pair The pair to search for
     * @param flushing_cg_counter Pointer to the flushing CG counter
     * @param probe_output_buffer Buffer of the matched probe pair sequence
     * @param contained_output_buffer Buffer of the matched contained pair sequence
     * @param num_matches Size of the output sequence
     * @param probe_output_begin Beginning of the output sequence of the matched probe pairs
     * @param contained_output_begin Beginning of the output sequence of the matched contained
     * pairs
     * @param pair_equal The binary callable used to compare two pairs for equality
     */
    template <uint32_t buffer_size,
              typename FlushingCG,
              typename atomicT,
              typename OutputIt1,
              typename OutputIt2,
              typename PairEqual,
              typename ParentCG>
    __device__ __forceinline__ void pair_retrieve_outer(
      FlushingCG flushing_cg,
      cooperative_groups::thread_block_tile<ProbeSequence::cg_size, ParentCG> probing_cg,
      value_type const& pair,
      uint32_t* flushing_cg_counter,
      value_type* probe_output_buffer,
      value_type* contained_output_buffer,
      atomicT* num_matches,
      OutputIt1 probe_output_begin,
      OutputIt2 contained_output_begin,
      PairEqual pair_equal) noexcept;

   private:
    using device_view_base<device_view_impl>::impl_;  ///< Implementation detail of `device_view`
  };  // class device_view

  /**
   * @brief Return the raw pointer of the hash map slots.
   *
   * @return Raw pointer of the hash map slots
   */
  value_type* raw_slots() noexcept
  {
    // Unsafe access to the slots stripping away their atomic-ness to allow non-atomic access.
    // TODO: to be replace by atomic_ref when it's ready
    return reinterpret_cast<value_type*>(slots_.get());
  }

  /**
   * @brief Return the raw pointer of the hash map slots.
   *
   * @return Raw pointer of the hash map slots
   */
  value_type const* raw_slots() const noexcept
  {
    // Unsafe access to the slots stripping away their atomic-ness to allow non-atomic access.
    // TODO: to be replace by atomic_ref when it's ready
    return reinterpret_cast<value_type const*>(slots_.get());
  }

  /**
   * @brief Gets the maximum number of elements the hash map can hold.
   *
   * @return The maximum number of elements the hash map can hold
   */
  std::size_t get_capacity() const noexcept { return capacity_; }

  /**
   * @brief Gets the number of elements in the hash map.
   *
   * @param stream CUDA stream used to get the number of inserted elements
   * @return The number of elements in the map
   */
  std::size_t get_size(cudaStream_t stream = 0) const noexcept;

  /**
   * @brief Gets the load factor of the hash map.
   *
   * @param stream CUDA stream used to get the load factor
   * @return The load factor of the hash map
   */
  float get_load_factor(cudaStream_t stream = 0) const noexcept;

  /**
   * @brief Gets the sentinel value used to represent an empty key slot.
   *
   * @return The sentinel value used to represent an empty key slot
   */
  Key get_empty_key_sentinel() const noexcept { return empty_key_sentinel_; }

  /**
   * @brief Gets the sentinel value used to represent an empty value slot.
   *
   * @return The sentinel value used to represent an empty value slot
   */
  Value get_empty_value_sentinel() const noexcept { return empty_value_sentinel_; }

  /**
   * @brief Constructs a device_view object based on the members of the `static_multimap`
   * object.
   *
   * @return A device_view object based on the members of the `static_multimap` object
   */
  device_view get_device_view() const noexcept
  {
    return device_view(slots_.get(),
                       capacity_,
                       empty_key<Key>{empty_key_sentinel_},
                       empty_value<Value>{empty_value_sentinel_});
  }

  /**
   * @brief Constructs a device_mutable_view object based on the members of the
   * `static_multimap` object
   *
   * @return A device_mutable_view object based on the members of the `static_multimap` object
   */
  device_mutable_view get_device_mutable_view() const noexcept
  {
    return device_mutable_view(slots_.get(),
                               capacity_,
                               empty_key<Key>{empty_key_sentinel_},
                               empty_value<Value>{empty_value_sentinel_});
  }

 private:
  std::size_t capacity_{};        ///< Total number of slots
  Key empty_key_sentinel_{};      ///< Key value that represents an empty slot
  Value empty_value_sentinel_{};  ///< Initial value of empty slot
  allocator_type allocator_{};    ///< Allocator used to allocate slots
  slot_deleter delete_slots_;     ///< Custom slots deleter
  std::unique_ptr<pair_atomic_type, slot_deleter> slots_{};  ///< Pointer to flat slots storage
};  // class static_multimap

}  // namespace cuco

#include <cuco/detail/static_multimap/device_view_impl.inl>
#include <cuco/detail/static_multimap/static_multimap.inl>
