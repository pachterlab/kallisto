/*
 * Copyright (c) 2020-2025, NVIDIA CORPORATION.
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
#include <cuco/detail/static_map_kernels.cuh>
#include <cuco/hash_functions.cuh>
#include <cuco/pair.cuh>
#include <cuco/static_map_ref.cuh>
#include <cuco/types.cuh>
#include <cuco/utility/allocator.hpp>
#include <cuco/utility/cuda_thread_scope.cuh>
#include <cuco/utility/traits.hpp>

#include <cuda/std/atomic>
#include <cuda/std/functional>
#include <cuda/std/utility>
#include <cuda/stream_ref>

#if defined(CUCO_HAS_CUDA_BARRIER)
#include <cuda/barrier>
#endif

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

namespace cuco {
/**
 * @brief A GPU-accelerated, unordered, associative container of key-value pairs with unique keys.
 *
 * The `static_map` supports two types of operations:
 * - Host-side "bulk" operations
 * - Device-side "singular" operations
 *
 * The host-side bulk operations include `insert`, `contains`, etc. These APIs should be used when
 * there are a large number of keys to modify or lookup. For example, given a range of keys
 * specified by device-accessible iterators, the bulk `insert` function inserts all keys into
 * the map.
 *
 * The singular device-side operations allow individual threads (or cooperative groups) to perform
 * independent modify or lookup operations from device code. These operations are accessed through
 * non-owning, trivially copyable reference types (or "ref"). Users can combine any arbitrary
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
          class ProbingScheme      = cuco::linear_probing<4,  // CG size
                                                          cuco::default_hash_function<Key>>,
          class Allocator          = cuco::cuda_allocator<cuco::pair<Key, T>>,
          class Storage            = cuco::storage<1>>
class static_map {
  static_assert(sizeof(Key) <= 8, "Container does not support key types larger than 8 bytes.");

  static_assert(sizeof(T) <= 8, "Container does not support payload types larger than 8 bytes.");

  static_assert(cuco::is_bitwise_comparable_v<T>,
                "Mapped type must have unique object representations or have been explicitly "
                "declared as safe for bitwise comparison via specialization of "
                "cuco::is_bitwise_comparable_v<T>.");

  using impl_type = detail::open_addressing_impl<Key,
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
  using ref_type = cuco::static_map_ref<key_type,
                                        mapped_type,
                                        thread_scope,
                                        key_equal,
                                        probing_scheme_type,
                                        storage_ref_type,
                                        Operators...>;  ///< Non-owning container ref type

  static_map(static_map const&)            = delete;
  static_map& operator=(static_map const&) = delete;

  static_map(static_map&&) = default;  ///< Move constructor

  /**
   * @brief Replaces the contents of the container with another container.
   *
   * @return Reference of the current map object
   */
  static_map& operator=(static_map&&) = default;
  ~static_map()                       = default;

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
  constexpr static_map(Extent capacity,
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
  constexpr static_map(Extent n,
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
  constexpr static_map(Extent capacity,
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
   * @brief Inserts all keys in the range `[first, last)` and returns the number of successful
   * insertions.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `insert_async`.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_map<K, V>::value_type></tt> is `true`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream CUDA stream used for insert
   *
   * @return Number of successful insertions
   */
  template <typename InputIt>
  size_type insert(InputIt first,
                   InputIt last,
                   cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously inserts all keys in the range `[first, last)`.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_map<K, V>::value_type></tt> is `true`
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
   * @brief Asynchronously inserts all elements in the range `[first, last)`.
   *
   * @note: For a given element `*(first + i)`, if the container doesn't already contain an element
   * with an equivalent key, inserts the element at a location pointed by `iter` and writes
   * `iter` to `found_begin + i` and writes `true` to `inserted_begin + i`. Otherwise, finds the
   * location of the equivalent element, `iter` and writes `iter` to `found_begin + i` and writes
   * `false` to `inserted_begin + i`.
   *
   * @tparam InputIt Device accessible random access input iterator
   * @tparam FoundIt Device accessible random access output iterator whose `value_type`
   * is constructible from `map::iterator` type
   * @tparam InsertedIt Device accessible random access output iterator whose `value_type`
   * is constructible from `bool`
   *
   * @param first Beginning of the sequence of elements
   * @param last End of the sequence of elements
   * @param found_begin Beginning of the sequence of elements found for each key
   * @param inserted_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt, typename FoundIt, typename InsertedIt>
  void insert_and_find_async(InputIt first,
                             InputIt last,
                             FoundIt found_begin,
                             InsertedIt inserted_begin,
                             cuda::stream_ref stream = cuda::stream_ref{
                               cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief Inserts all elements in the range `[first, last)`.
   *
   * @note: For a given element `*(first + i)`, if the container doesn't already contain an element
   * with an equivalent key, inserts the element at a location pointed by `iter` and writes
   * `iter` to `found_begin + i` and writes `true` to `inserted_begin + i`. Otherwise, finds the
   * location of the equivalent element, `iter` and writes `iter` to `found_begin + i` and writes
   * `false` to `inserted_begin + i`.
   *
   * @tparam InputIt Device accessible random access input iterator
   * @tparam FoundIt Device accessible random access output iterator whose `value_type`
   * is constructible from `map::iterator` type
   * @tparam InsertedIt Device accessible random access output iterator whose `value_type`
   * is constructible from `bool`
   *
   * @param first Beginning of the sequence of elements
   * @param last End of the sequence of elements
   * @param found_begin Beginning of the sequence of elements found for each key
   * @param inserted_begin Beginning of the sequence of booleans for the presence of each key
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt, typename FoundIt, typename InsertedIt>
  void insert_and_find(InputIt first,
                       InputIt last,
                       FoundIt found_begin,
                       InsertedIt inserted_begin,
                       cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief For any key-value pair `{k, v}` in the range `[first, last)`, if a key equivalent to `k`
   * already exists in the container, assigns `v` to the mapped_type corresponding to the key `k`.
   * If the key does not exist, inserts the pair as if by insert.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `insert_or_assign_async`.
   * @note If multiple pairs in `[first, last)` compare equal, it is unspecified which pair is
   * inserted or assigned.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_map<K, V>::value_type></tt> is `true`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt>
  void insert_or_assign(InputIt first,
                        InputIt last,
                        cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief For any key-value pair `{k, v}` in the range `[first, last)`, if a key equivalent to `k`
   * already exists in the container, assigns `v` to the mapped_type corresponding to the key `k`.
   * If the key does not exist, inserts the pair as if by insert.
   *
   * @note If multiple pairs in `[first, last)` compare equal, it is unspecified which pair is
   * inserted or assigned.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_map<K, V>::value_type></tt> is `true`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt>
  void insert_or_assign_async(InputIt first,
                              InputIt last,
                              cuda::stream_ref stream = cuda::stream_ref{
                                cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief For any key-value pair `{k, v}` in the range `[first, first + n)`, if a key equivalent
   * to `k` already exists in the container, then binary operation is applied using `op` callable
   * object on the existing value at slot and the element to insert. If the key does not exist,
   * inserts the pair as if by insert.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `insert_or_apply_async`.
   * @note Callable object to perform binary operation should be able to invoke as
   *  Op(cuda::atomic_ref<T, Scope>, T>)
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_map<K, V>::value_type></tt> is `true`
   * @tparam Op Callable type used to peform `apply` operation.
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param op Callable object to perform apply operation.
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt, typename Op>
  void insert_or_apply(InputIt first,
                       InputIt last,
                       Op op,
                       cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief For any key-value pair `{k, v}` in the range `[first, first + n)`, if a key equivalent
   * to `k` already exists in the container, then binary operation is applied using `op` callable
   * object on the existing value at slot and the element to insert. If the key does not exist,
   * inserts the pair as if by insert.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `insert_or_apply_async`.
   * @note Callable object to perform binary operation should be able to invoke as
   *  Op(cuda::atomic_ref<T, Scope>, T>)
   * @note There could be performance improvements if `init` value passed here equals to the
   *  `sentinel value` of the map.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_map<K, V>::value_type></tt> is `true`
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type used to peform `apply` operation.
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param init The identity value of the op
   * @param op Callable object to perform apply operation.
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt, typename Init, typename Op>
  void insert_or_apply(InputIt first,
                       InputIt last,
                       Init init,
                       Op op,
                       cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief For any key-value pair `{k, v}` in the range `[first, first + n)`, if a key equivalent
   * to `k` already exists in the container, then binary operation is applied using `op` callable
   * object on the existing value at slot and the element to insert. If the key does not exist,
   * inserts the pair as if by insert.
   *
   * @note Callable object to perform binary operation should be able to invoke as
   *  Op(cuda::atomic_ref<T, Scope>, T>)
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_map<K, V>::value_type></tt> is `true`
   * @tparam Op Callable type used to peform `apply` operation.
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param op Callable object to perform apply operation.
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt, typename Op>
  void insert_or_apply_async(InputIt first,
                             InputIt last,
                             Op op,
                             cuda::stream_ref stream = cuda::stream_ref{
                               cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief For any key-value pair `{k, v}` in the range `[first, first + n)`, if a key equivalent
   * to `k` already exists in the container, then binary operation is applied using `op` callable
   * object on the existing value at slot and the element to insert. If the key does not exist,
   * inserts the pair as if by insert.
   *
   * @note Callable object to perform binary operation should be able to invoke as
   *  Op(cuda::atomic_ref<T, Scope>, T>)
   * @note There could be performance improvements if `init` value passed here equals to the
   *  `sentinel value` of the map.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * static_map<K, V>::value_type></tt> is `true`
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type used to peform `apply` operation.
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param init The identity value of the op
   * @param op Callable object to perform apply operation.
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt,
            typename Init,
            typename Op,
            typename = std::enable_if_t<std::is_convertible_v<Init, T>>>
  void insert_or_apply_async(InputIt first,
                             InputIt last,
                             Init init,
                             Op op,
                             cuda::stream_ref stream = cuda::stream_ref{
                               cudaStream_t{nullptr}}) noexcept;

  /**
   * @brief Erases keys in the range `[first, last)`.
   *
   * @note For each key `k` in `[first, last)`, if contains(k) returns true, removes `k` and it's
   * associated value from the map. Else, no effect.
   *
   * @note This function synchronizes `stream`.
   *
   * @note Side-effects:
   *  - `contains(k) == false`
   *  - `find(k) == end()`
   *  - `insert({k,v}) == true`
   *  - `size()` is reduced by the total number of erased keys
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `key_type`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream Stream used for executing the kernels
   *
   * @throw std::runtime_error if a unique erased key sentinel value was not
   * provided at construction
   */
  template <typename InputIt>
  void erase(InputIt first,
             InputIt last,
             cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously erases keys in the range `[first, last)`.
   *
   * @note For each key `k` in `[first, last)`, if contains(k) returns true, removes `k` and it's
   * associated value from the map. Else, no effect.
   *
   * @note Side-effects:
   *  - `contains(k) == false`
   *  - `find(k) == end()`
   *  - `insert({k,v}) == true`
   *  - `size()` is reduced by the total number of erased keys
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `key_type`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stream Stream used for executing the kernels
   *
   * @throw std::runtime_error if a unique erased key sentinel value was not
   * provided at construction
   */
  template <typename InputIt>
  void erase_async(InputIt first,
                   InputIt last,
                   cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

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
   * false, stores false to `(output_begin + i)`.
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
   * false, stores false to `(output_begin + i)`.
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
   * @brief For all keys in the range `[first, last)`, asynchronously finds an element with key
   * equivalent to the query key.
   *
   * @note If the key `*(first + i)` has a matched `element` in the map, copies the payload of
   * `element` to `(output_begin + i)`. Else, copies the empty value sentinel.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam ProbeEqual Binary callable equal type
   * @tparam ProbeHash Unary callable hasher type
   * @tparam OutputIt Device accessible output iterator assignable from the map's `mapped_type`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param probe_equal The binary function to compare map keys and probe keys for equality
   * @param probe_hash The unary function to hash probe keys
   * @param output_begin Beginning of the sequence of elements retrieved for each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename ProbeEqual, typename ProbeHash, typename OutputIt>
  void find_async(InputIt first,
                  InputIt last,
                  ProbeEqual const& probe_equal,
                  ProbeHash const& probe_hash,
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
   * @tparam ProbeEqual Binary callable equal type
   * @tparam ProbeHash Unary callable hasher type
   * @tparam OutputIt Device accessible output iterator
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param probe_equal The binary function to compare map keys and probe keys for equality
   * @param probe_hash The unary function to hash probe keys
   * @param output_begin Beginning of the sequence of matches retrieved for each key
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt,
            typename StencilIt,
            typename Predicate,
            typename ProbeEqual,
            typename ProbeHash,
            typename OutputIt>
  void find_if_async(InputIt first,
                     InputIt last,
                     StencilIt stencil,
                     Predicate pred,
                     ProbeEqual const& probe_equal,
                     ProbeHash const& probe_hash,
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
   * @brief Counts the occurrences of keys in `[first, last)` contained in the map
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
   * @brief Retrieves the matched key-value pair in the map corresponding to all probe keys in the
   * range
   * `[first, last)`
   *
   * If key `k = *(first + i)` has a match `m` in the map, copies a `cuco::pair{k, m}` to
   * unspecified locations in `[output_begin, output_end)`. Else, does nothing.
   *
   * @note This function synchronizes the given stream.
   * @note Behavior is undefined if the size of the output range exceeds
   * `std::distance(output_begin, output_end)`.
   * @note Behavior is undefined if the given key has multiple matches in the set.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputProbeIt Device accessible output iterator whose `value_type` can be constructed
   * from `ProbeKey`
   * @tparam OutputMatchIt Device accessible output iterator whose `value_type` can be constructed
   * from map's `value_type`
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
   * @brief Retrieves all of the keys and their associated values contained in the map
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
   * @brief Gets the maximum number of elements the hash map can hold.
   *
   * @return The maximum number of elements the hash map can hold
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
   * @return Device ref of the current `static_map` object
   */
  template <typename... Operators>
  [[nodiscard]] auto ref(Operators... ops) const noexcept;

 private:
  /**
   * @brief Dispatches to shared memory map kernel if `num_elements_per_thread > 2`, else
   * fallbacks to global memory map kernel.
   *
   * @tparam HasInit Boolean to dispatch based on init parameter
   * @tparam CGSize Number of threads in each CG
   * @tparam Allocator Allocator type used to created shared_memory map
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the `value_type` of the data structure
   * @tparam Init Type of init value convertible to payload type
   * @tparam Op Callable type used to peform `apply` operation.
   * @tparam Ref Type of non-owning device ref allowing access to storage
   *
   * @param first Beginning of the sequence of input elements
   * @param last End of the sequence of input elements
   * @param init The init value of the `op`
   * @param op Callable object to perform apply operation.
   * @param ref Non-owning container device ref used to access the slot storage
   * @param stream CUDA stream used for insert_or_apply operation
   */
  template <bool HasInit,
            int32_t CGSize,
            typename AllocatorType,
            typename InputIt,
            typename InitType,
            typename OpType,
            typename RefType>
  void dispatch_insert_or_apply(
    InputIt first, InputIt last, InitType init, OpType op, RefType ref, cuda::stream_ref stream);

  std::unique_ptr<impl_type> impl_;   ///< Static map implementation
  mapped_type empty_value_sentinel_;  ///< Sentinel value that indicates an empty payload
};

namespace experimental {
template <class Key,
          class T,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
class dynamic_map;
}

template <typename Key, typename Value, cuda::thread_scope Scope, typename Allocator>
class dynamic_map;

namespace legacy {

/**
 * @brief A GPU-accelerated, unordered, associative container of key-value
 * pairs with unique keys.
 *
 * Allows constant time concurrent inserts or concurrent find operations from threads in device
 * code. Concurrent insert and find are supported only if the pair type is packable (see
 * `cuco::detail::is_packable` constexpr).
 *
 * Current limitations:
 * - Requires keys and values that where `cuco::is_bitwise_comparable_v<T>` is true
 *    - Comparisons against the "sentinel" values will always be done with bitwise comparisons.
 * - Capacity is fixed and will not grow automatically
 * - Requires the user to specify sentinel values for both key and mapped value to indicate empty
 * slots
 * - Conditionally support concurrent insert and find operations
 *
 * The `static_map` supports two types of operations:
 * - Host-side "bulk" operations
 * - Device-side "singular" operations
 *
 * The host-side bulk operations include `insert`, `erase`, `find`, and `contains`. These
 * APIs should be used when there are a large number of keys to insert, erase or lookup
 * in the map. For example, given a range of keys specified by device-accessible
 * iterators, the bulk `insert` function will insert all keys into the map. Note that in order
 * for a `static_map` instance to support `erase`, the user must provide an `erased_key_sentinel`
 * which is distinct from the `empty_key_sentinel` at construction. If `erase` is called on a
 * `static_map` which was not constructed in this way, a runtime error will be generated.
 *
 * The singular device-side operations allow individual threads to perform
 * independent insert or find/contains operations from device code. These
 * operations are accessed through non-owning, trivially copyable "view" types:
 * `device_view` and `device_mutable_view`. The `device_view` class is an
 * immutable view that allows only non-modifying operations such as `find` or
 * `contains`. The `device_mutable_view` class only allows `insert` and `erase` operations.
 * The two types are separate to prevent erroneous concurrent insert/erase/find
 * operations. Note that the device-side `erase` may only be called if the corresponding
 * `device_mutable_view` was constructed with a user-provided `erased_key_sentinel`. It is
 * up to the user to ensure this condition is met.
 *
 * Example:
 * \code{.cpp}
 * int empty_key_sentinel = -1;
 * int empty_value_sentinel = -1;
 * int erased_key_sentinel = -2;
 *
 * // Constructs a map with 100,000 slots using -1 and -1 as the empty key/value
 * // sentinels. The supplied erased key sentinel of -2 must be a different value from the empty
 * // key sentinel. If erase functionality is not needed, you may elect to not supply an erased
 * // key sentinel to the constructor. Note the capacity is chosen knowing we will insert 50,000
 * // keys, for an load factor of 50%.
 * static_map<int, int> m{100'000, empty_key_sentinel, empty_value_sentinel, erased_value_sentinel};
 *
 * // Create a sequence of pairs {{0,0}, {1,1}, ... {i,i}}
 * thrust::device_vector<cuda::std::pair<int,int>> pairs(50,000);
 * thrust::transform(thrust::make_counting_iterator(0),
 *                   thrust::make_counting_iterator(pairs.size()),
 *                   pairs.begin(),
 *                   []__device__(auto i){ return cuco::pair{i,i}; };
 *
 *
 * // Inserts all pairs into the map
 * m.insert(pairs.begin(), pairs.end());
 *
 * // Get a `device_view` and passes it to a kernel where threads may perform
 * // `find/contains` lookups
 * kernel<<<...>>>(m.get_device_view());
 * \endcode
 *
 *
 * @tparam Key Arithmetic type used for key
 * @tparam Value Type of the mapped values
 * @tparam Scope The scope in which insert/find operations will be performed by
 * individual threads.
 * @tparam Allocator Type of allocator used for device storage
 */
template <typename Key,
          typename Value,
          cuda::thread_scope Scope = cuda::thread_scope_device,
          typename Allocator       = cuco::cuda_allocator<char>>
class static_map {
  static_assert(
    cuco::is_bitwise_comparable_v<Key>,
    "Key type must have unique object representations or have been explicitly declared as safe for "
    "bitwise comparison via specialization of cuco::is_bitwise_comparable_v<Key>.");

  static_assert(cuco::is_bitwise_comparable_v<Value>,
                "Value type must have unique object representations or have been explicitly "
                "declared as safe for bitwise comparison via specialization of "
                "cuco::is_bitwise_comparable_v<Value>.");

  friend class dynamic_map<Key, Value, Scope, Allocator>;  ///< Dynamic map as friend class

 public:
  using value_type         = cuco::pair<Key, Value>;            ///< Type of key/value pairs
  using key_type           = Key;                               ///< Key type
  using mapped_type        = Value;                             ///< Type of mapped values
  using atomic_key_type    = cuda::atomic<key_type, Scope>;     ///< Type of atomic keys
  using atomic_mapped_type = cuda::atomic<mapped_type, Scope>;  ///< Type of atomic mapped values
  using pair_atomic_type =
    cuco::pair<atomic_key_type,
               atomic_mapped_type>;  ///< Pair type of atomic key and atomic mapped value
  using slot_type           = pair_atomic_type;                  ///< Type of hash map slots
  using atomic_ctr_type     = cuda::atomic<std::size_t, Scope>;  ///< Atomic counter type
  using allocator_type      = Allocator;                         ///< Allocator type
  using slot_allocator_type = typename std::allocator_traits<Allocator>::template rebind_alloc<
    pair_atomic_type>;  ///< Type of the allocator to (de)allocate slots
  using counter_allocator_type = typename std::allocator_traits<Allocator>::template rebind_alloc<
    atomic_ctr_type>;  ///< Type of the allocator to (de)allocate atomic counters

#if !defined(CUCO_HAS_INDEPENDENT_THREADS)
  static_assert(atomic_key_type::is_always_lock_free,
                "A key type larger than 8B is supported for only sm_70 and up.");
  static_assert(atomic_mapped_type::is_always_lock_free,
                "A value type larger than 8B is supported for only sm_70 and up.");
#endif

  static_map(static_map const&) = delete;
  static_map(static_map&&)      = delete;

  static_map& operator=(static_map const&) = delete;
  static_map& operator=(static_map&&)      = delete;

  /**
   * @brief Indicates if concurrent insert/find is supported for the key/value types.
   *
   * @return Boolean indicating if concurrent insert/find is supported.
   */
  __host__ __device__ static constexpr bool supports_concurrent_insert_find() noexcept
  {
    return cuco::detail::is_packable<value_type>();
  }

  /**
   * @brief Constructs a statically sized map with the specified number of slots
   * and sentinel values.
   *
   * The capacity of the map is fixed. Insert operations will not automatically
   * grow the map. Attempting to insert equal to or more unique keys than the capacity
   * of the map results in undefined behavior (there should be at least one empty slot).
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
   * @param alloc Allocator used for allocating device storage
   * @param stream Stream used for executing the kernels
   */
  static_map(std::size_t capacity,
             empty_key<Key> empty_key_sentinel,
             empty_value<Value> empty_value_sentinel,
             Allocator const& alloc = Allocator{},
             cudaStream_t stream    = 0);

  /**
   * @brief Constructs a fixed-size map with erase capability.
   * empty_key_sentinel and erased_key_sentinel must be different values.
   *
   * @throw std::runtime error if the empty key sentinel and erased key sentinel
   * are the same value
   *
   * @param capacity The total number of slots in the map
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param empty_value_sentinel The reserved mapped value for empty slots
   * @param erased_key_sentinel The reserved value to denote erased slots
   * @param alloc Allocator used for allocating device storage
   * @param stream Stream used for executing the kernels
   */
  static_map(std::size_t capacity,
             empty_key<Key> empty_key_sentinel,
             empty_value<Value> empty_value_sentinel,
             erased_key<Key> erased_key_sentinel,
             Allocator const& alloc = Allocator{},
             cudaStream_t stream    = 0);

  /**
   * @brief Destroys the map and frees its contents.
   *
   */
  ~static_map();

  /**
   * @brief Inserts all key/value pairs in the range `[first, last)`.
   *
   * This function synchronizes `stream`.
   *
   * If multiple keys in `[first, last)` compare equal, it is unspecified which
   * element is inserted.
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `value_type`
   * @tparam Hash Unary callable type
   * @tparam KeyEqual Binary callable type
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param hash The unary function to apply to hash each key
   * @param key_equal The binary function to compare two keys for equality
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt,
            typename Hash     = cuco::default_hash_function<key_type>,
            typename KeyEqual = cuda::std::equal_to<key_type>>
  void insert(InputIt first,
              InputIt last,
              Hash hash           = Hash{},
              KeyEqual key_equal  = KeyEqual{},
              cudaStream_t stream = 0);

  /**
   * @brief Inserts key/value pairs in the range `[first, last)` if `pred`
   * of the corresponding stencil returns true.
   *
   * The key/value pair `*(first + i)` is inserted if `pred( *(stencil + i) )` returns true.
   *
   * @tparam InputIt Device accessible random access iterator whose `value_type` is
   * convertible to the map's `value_type`
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam Hash Unary callable type
   * @tparam KeyEqual Binary callable type
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param hash The unary function to hash each key
   * @param key_equal The binary function to compare two keys for equality
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt,
            typename StencilIt,
            typename Predicate,
            typename Hash     = cuco::default_hash_function<key_type>,
            typename KeyEqual = cuda::std::equal_to<key_type>>
  void insert_if(InputIt first,
                 InputIt last,
                 StencilIt stencil,
                 Predicate pred,
                 Hash hash           = Hash{},
                 KeyEqual key_equal  = KeyEqual{},
                 cudaStream_t stream = 0);

  /**
   * @brief Erases keys in the range `[first, last)`.
   *
   * For each key `k` in `[first, last)`, if `contains(k) == true), removes `k` and it's
   * associated value from the map. Else, no effect.
   *
   *  Side-effects:
   *  - `contains(k) == false`
   *  - `find(k) == end()`
   *  - `insert({k,v}) == true`
   *  - `get_size()` is reduced by the total number of erased keys
   *
   * This function synchronizes `stream`.
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `value_type`
   * @tparam Hash Unary callable type
   * @tparam KeyEqual Binary callable type
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param hash The unary function to apply to hash each key
   * @param key_equal The binary function to compare two keys for equality
   * @param stream Stream used for executing the kernels
   *
   * @throw std::runtime_error if a unique erased key sentinel value was not
   * provided at construction
   */
  template <typename InputIt,
            typename Hash     = cuco::default_hash_function<key_type>,
            typename KeyEqual = cuda::std::equal_to<key_type>>
  void erase(InputIt first,
             InputIt last,
             Hash hash           = Hash{},
             KeyEqual key_equal  = KeyEqual{},
             cudaStream_t stream = 0);

  /**
   * @brief Finds the values corresponding to all keys in the range `[first, last)`.
   *
   * If the key `*(first + i)` exists in the map, copies its associated value to `(output_begin +
   * i)`. Else, copies the empty value sentinel.
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `key_type`
   * @tparam OutputIt Device accessible output iterator whose `value_type` is
   * convertible to the map's `mapped_type`
   * @tparam Hash Unary callable type
   * @tparam KeyEqual Binary callable type
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of values retrieved for each key
   * @param hash The unary function to apply to hash each key
   * @param key_equal The binary function to compare two keys for equality
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt,
            typename OutputIt,
            typename Hash     = cuco::default_hash_function<key_type>,
            typename KeyEqual = cuda::std::equal_to<key_type>>
  void find(InputIt first,
            InputIt last,
            OutputIt output_begin,
            Hash hash           = Hash{},
            KeyEqual key_equal  = KeyEqual{},
            cudaStream_t stream = 0);

  /**
   * @brief Retrieves all of the keys and their associated values.
   *
   * The order in which keys are returned is implementation defined and not guaranteed to be
   * consistent between subsequent calls to `retrieve_all`.
   *
   * Behavior is undefined if the range beginning at `keys_out` or `values_out` is less than
   * `get_size()`
   *
   * @tparam KeyOut Device accessible random access output iterator whose `value_type` is
   * convertible from `key_type`.
   * @tparam ValueOut Device accessible random access output iterator whose `value_type` is
   * convertible from `mapped_type`.
   * @param keys_out Beginning output iterator for keys
   * @param values_out Beginning output iterator for values
   * @param stream CUDA stream used for this operation
   * @return Pair of iterators indicating the last elements in the output
   */
  template <typename KeyOut, typename ValueOut>
  std::pair<KeyOut, ValueOut> retrieve_all(KeyOut keys_out,
                                           ValueOut values_out,
                                           cudaStream_t stream = 0) const;

  /**
   * @brief Indicates whether the keys in the range `[first, last)` are contained in the map.
   *
   * Writes a `bool` to `(output + i)` indicating if the key `*(first + i)` exists in the map.
   *
   * Hash should be callable with both <tt>std::iterator_traits<InputIt>::value_type</tt> and Key
   * type. <tt>std::invoke_result<KeyEqual, std::iterator_traits<InputIt>::value_type, Key></tt>
   * must be well-formed.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   * @tparam Hash Unary callable type
   * @tparam KeyEqual Binary callable type
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param hash The unary function to apply to hash each key
   * @param key_equal The binary function to compare two keys for equality
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt,
            typename OutputIt,
            typename Hash     = cuco::default_hash_function<key_type>,
            typename KeyEqual = cuda::std::equal_to<key_type>>
  void contains(InputIt first,
                InputIt last,
                OutputIt output_begin,
                Hash hash           = Hash{},
                KeyEqual key_equal  = KeyEqual{},
                cudaStream_t stream = 0) const;

 private:
  class device_view_base {
   protected:
    // Import member type definitions from `static_map`
    using value_type     = value_type;
    using key_type       = Key;
    using mapped_type    = Value;
    using iterator       = pair_atomic_type*;
    using const_iterator = pair_atomic_type const*;
    using slot_type      = slot_type;

    Key empty_key_sentinel_{};      ///< Key value that represents an empty slot
    Key erased_key_sentinel_{};     ///< Key value that represents an erased slot
    Value empty_value_sentinel_{};  ///< Initial Value of empty slot
    pair_atomic_type* slots_{};     ///< Pointer to flat slots storage
    std::size_t capacity_{};        ///< Total number of slots

    __host__ __device__ device_view_base(pair_atomic_type* slots,
                                         std::size_t capacity,
                                         empty_key<Key> empty_key_sentinel,
                                         empty_value<Value> empty_value_sentinel) noexcept
      : slots_{slots},
        capacity_{capacity},
        empty_key_sentinel_{empty_key_sentinel.value},
        erased_key_sentinel_{empty_key_sentinel.value},
        empty_value_sentinel_{empty_value_sentinel.value}
    {
    }

    __host__ __device__ device_view_base(pair_atomic_type* slots,
                                         std::size_t capacity,
                                         empty_key<Key> empty_key_sentinel,
                                         empty_value<Value> empty_value_sentinel,
                                         erased_key<Key> erased_key_sentinel) noexcept
      : slots_{slots},
        capacity_{capacity},
        empty_key_sentinel_{empty_key_sentinel.value},
        erased_key_sentinel_{erased_key_sentinel.value},
        empty_value_sentinel_{empty_value_sentinel.value}
    {
    }

    /**
     * @brief Returns the initial slot for a given key `k`
     *
     * @tparam ProbeKey Probe key type
     * @tparam Hash Unary callable type
     *
     * @param k The key to get the slot for
     * @param hash The unary callable used to hash the key
     * @return Pointer to the initial slot for `k`
     */
    template <typename ProbeKey, typename Hash>
    __device__ iterator initial_slot(ProbeKey const& k, Hash hash) noexcept
    {
      return &slots_[hash(k) % capacity_];
    }

    /**
     * @brief Returns the initial slot for a given key `k`
     *
     * @tparam ProbeKey Probe key type
     * @tparam Hash Unary callable type
     *
     * @param k The key to get the slot for
     * @param hash The unary callable used to hash the key
     * @return Pointer to the initial slot for `k`
     */
    template <typename ProbeKey, typename Hash>
    __device__ const_iterator initial_slot(ProbeKey const& k, Hash hash) const noexcept
    {
      return &slots_[hash(k) % capacity_];
    }

    /**
     * @brief Returns the initial slot for a given key `k`
     *
     * To be used for Cooperative Group based probing.
     *
     * @tparam CG Cooperative Group type
     * @tparam ProbeKey Probe key type
     * @tparam Hash Unary callable type
     *
     * @param g the Cooperative Group for which the initial slot is needed
     * @param k The key to get the slot for
     * @param hash The unary callable used to hash the key
     * @return Pointer to the initial slot for `k`
     */
    template <typename CG, typename ProbeKey, typename Hash>
    __device__ iterator initial_slot(CG g, ProbeKey const& k, Hash hash) noexcept
    {
      return &slots_[(hash(k) + g.thread_rank()) % capacity_];
    }

    /**
     * @brief Returns the initial slot for a given key `k`
     *
     * To be used for Cooperative Group based probing.
     *
     * @tparam CG Cooperative Group type
     * @tparam ProbeKey Probe key type
     * @tparam Hash Unary callable type
     *
     * @param g the Cooperative Group for which the initial slot is needed
     * @param k The key to get the slot for
     * @param hash The unary callable used to hash the key
     * @return Pointer to the initial slot for `k`
     */
    template <typename CG, typename ProbeKey, typename Hash>
    __device__ const_iterator initial_slot(CG g, ProbeKey const& k, Hash hash) const noexcept
    {
      return &slots_[(hash(k) + g.thread_rank()) % capacity_];
    }

    /**
     * @brief Given a slot `s`, returns the next slot.
     *
     * If `s` is the last slot, wraps back around to the first slot.
     *
     * @param s The slot to advance
     * @return The next slot after `s`
     */
    __device__ iterator next_slot(iterator s) noexcept { return (++s < end()) ? s : begin_slot(); }

    /**
     * @brief Given a slot `s`, returns the next slot.
     *
     * If `s` is the last slot, wraps back around to the first slot.
     *
     * @param s The slot to advance
     * @return The next slot after `s`
     */
    __device__ const_iterator next_slot(const_iterator s) const noexcept
    {
      return (++s < end()) ? s : begin_slot();
    }

    /**
     * @brief Given a slot `s`, returns the next slot.
     *
     * If `s` is the last slot, wraps back around to the first slot. To
     * be used for Cooperative Group based probing.
     *
     * @tparam CG The Cooperative Group type
     * @param g The Cooperative Group for which the next slot is needed
     * @param s The slot to advance
     * @return The next slot after `s`
     */
    template <typename CG>
    __device__ iterator next_slot(CG g, iterator s) noexcept
    {
      uint32_t index = s - slots_;
      return &slots_[(index + g.size()) % capacity_];
    }

    /**
     * @brief Given a slot `s`, returns the next slot.
     *
     * If `s` is the last slot, wraps back around to the first slot. To
     * be used for Cooperative Group based probing.
     *
     * @tparam CG The Cooperative Group type
     * @param g The Cooperative Group for which the next slot is needed
     * @param s The slot to advance
     * @return The next slot after `s`
     */
    template <typename CG>
    __device__ const_iterator next_slot(CG g, const_iterator s) const noexcept
    {
      uint32_t index = s - slots_;
      return &slots_[(index + g.size()) % capacity_];
    }

    /**
     * @brief Initializes the given array of slots to the specified values given by `k` and `v`
     * using the threads in the group `g`.
     *
     * @note This function synchronizes the group `g`.
     *
     * @tparam CG The type of the cooperative thread group
     * @param g The cooperative thread group used to initialize the slots
     * @param slots Pointer to the array of slots to initialize
     * @param num_slots Number of slots to initialize
     * @param k The desired key value for each slot
     * @param v The desired mapped value for each slot
     */

    template <typename CG>
    __device__ static void initialize_slots(
      CG g, pair_atomic_type* slots, std::size_t num_slots, Key k, Value v)
    {
      auto tid = g.thread_rank();
      while (tid < num_slots) {
        new (&slots[tid].first) atomic_key_type{k};
        new (&slots[tid].second) atomic_mapped_type{v};
        tid += g.size();
      }
      g.sync();
    }

   public:
    /**
     * @brief Gets slots array.
     *
     * @return Slots array
     */
    __host__ __device__ pair_atomic_type* get_slots() noexcept { return slots_; }

    /**
     * @brief Gets slots array.
     *
     * @return Slots array
     */
    __host__ __device__ pair_atomic_type const* get_slots() const noexcept { return slots_; }

    /**
     * @brief Gets the maximum number of elements the hash map can hold.
     *
     * @return The maximum number of elements the hash map can hold
     */
    __host__ __device__ std::size_t get_capacity() const noexcept { return capacity_; }

    /**
     * @brief Gets the sentinel value used to represent an empty key slot.
     *
     * @return The sentinel value used to represent an empty key slot
     */
    __host__ __device__ Key get_empty_key_sentinel() const noexcept { return empty_key_sentinel_; }

    /**
     * @brief Gets the sentinel value used to represent an empty value slot.
     *
     * @return The sentinel value used to represent an empty value slot
     */
    __host__ __device__ Value get_empty_value_sentinel() const noexcept
    {
      return empty_value_sentinel_;
    }

    __host__ __device__ Key get_erased_key_sentinel() const noexcept
    {
      return erased_key_sentinel_;
    }

    /**
     * @brief Returns iterator to the first slot.
     *
     * @note Unlike `std::map::begin()`, the `begin_slot()` iterator does _not_ point to the first
     * occupied slot. Instead, it refers to the first slot in the array of contiguous slot storage.
     * Iterating from `begin_slot()` to `end_slot()` will iterate over all slots, including those
     * both empty and filled.
     *
     * There is no `begin()` iterator to avoid confusion as it is not possible to provide an
     * iterator over only the filled slots.
     *
     * @return Iterator to the first slot
     */
    __device__ iterator begin_slot() noexcept { return slots_; }

    /**
     * @brief Returns iterator to the first slot.
     *
     * @note Unlike `std::map::begin()`, the `begin_slot()` iterator does _not_ point to the first
     * occupied slot. Instead, it refers to the first slot in the array of contiguous slot storage.
     * Iterating from `begin_slot()` to `end_slot()` will iterate over all slots, including those
     * both empty and filled.
     *
     * There is no `begin()` iterator to avoid confusion as it is not possible to provide an
     * iterator over only the filled slots.
     *
     * @return Iterator to the first slot
     */
    __device__ const_iterator begin_slot() const noexcept { return slots_; }

    /**
     * @brief Returns a const_iterator to one past the last slot.
     *
     * @return A const_iterator to one past the last slot
     */
    __host__ __device__ const_iterator end_slot() const noexcept { return slots_ + capacity_; }

    /**
     * @brief Returns an iterator to one past the last slot.
     *
     * @return An iterator to one past the last slot
     */
    __host__ __device__ iterator end_slot() noexcept { return slots_ + capacity_; }

    /**
     * @brief Returns a const_iterator to one past the last slot.
     *
     * `end()` calls `end_slot()` and is provided for convenience for those familiar with checking
     * an iterator returned from `find()` against the `end()` iterator.
     *
     * @return A const_iterator to one past the last slot
     */
    __host__ __device__ const_iterator end() const noexcept { return end_slot(); }

    /**
     * @brief Returns an iterator to one past the last slot.
     *
     * `end()` calls `end_slot()` and is provided for convenience for those familiar with checking
     * an iterator returned from `find()` against the `end()` iterator.
     *
     * @return An iterator to one past the last slot
     */
    __host__ __device__ iterator end() noexcept { return end_slot(); }
  };

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
   * cuco::static_map<int,int> m{100'000, -1, -1};
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
  class device_mutable_view : public device_view_base {
   public:
    using value_type  = typename device_view_base::value_type;   ///< Type of key/value pairs
    using key_type    = typename device_view_base::key_type;     ///< Key type
    using mapped_type = typename device_view_base::mapped_type;  ///< Type of the mapped values
    using iterator =
      typename device_view_base::iterator;  ///< Type of the forward iterator to `value_type`
    using const_iterator =
      typename device_view_base::const_iterator;  ///< Type of the forward iterator to `const
                                                  ///< value_type`
    using slot_type = typename device_view_base::slot_type;  ///< Type of hash map slots

    /**
     * @brief Constructs a mutable view of the first `capacity` slots of the
     * slots array pointed to by `slots`.
     *
     * @param slots Pointer to beginning of initialized slots array
     * @param capacity The number of slots viewed by this object
     * @param empty_key_sentinel The reserved value for keys to represent empty slots
     * @param empty_value_sentinel The reserved value for mapped values to
     * represent empty slots
     */
    __host__ __device__ device_mutable_view(pair_atomic_type* slots,
                                            std::size_t capacity,
                                            empty_key<Key> empty_key_sentinel,
                                            empty_value<Value> empty_value_sentinel) noexcept
      : device_view_base{slots, capacity, empty_key_sentinel, empty_value_sentinel}
    {
    }

    /**
     * @brief Constructs a mutable view of the first `capacity` slots of the
     * slots array pointed to by `slots`.
     *
     * @param slots Pointer to beginning of initialized slots array
     * @param capacity The number of slots viewed by this object
     * @param empty_key_sentinel The reserved value for keys to represent empty slots
     * @param empty_value_sentinel The reserved value for mapped values to represent empty slots
     * @param erased_key_sentinel The reserved value for keys to represent erased slots
     */
    __host__ __device__ device_mutable_view(pair_atomic_type* slots,
                                            std::size_t capacity,
                                            empty_key<Key> empty_key_sentinel,
                                            empty_value<Value> empty_value_sentinel,
                                            erased_key<Key> erased_key_sentinel) noexcept
      : device_view_base{
          slots, capacity, empty_key_sentinel, empty_value_sentinel, erased_key_sentinel}
    {
    }

   private:
    /**
     * @brief Enumeration of the possible results of attempting to insert into a hash slot.
     */
    enum class insert_result {
      CONTINUE,  ///< Insert did not succeed, continue trying to insert
      SUCCESS,   ///< New pair inserted successfully
      DUPLICATE  ///< Insert did not succeed, key is already present
    };

    /**
     * @brief Inserts the specified key/value pair with one single CAS operation.
     *
     * @tparam KeyEqual Binary callable type
     * @param current_slot The slot to insert
     * @param insert_pair The pair to insert
     * @param key_equal The binary callable used to compare two keys for
     * equality
     * @param expected_key The expected value of the key in the target slot
     * @return An insert result from the `insert_resullt` enumeration.
     */
    template <typename KeyEqual>
    __device__ insert_result packed_cas(iterator current_slot,
                                        value_type const& insert_pair,
                                        KeyEqual key_equal,
                                        Key expected_key) noexcept;

    /**
     * @brief Inserts the specified key/value pair with two back-to-back CAS operations.
     *
     * @tparam KeyEqual Binary callable type
     * @param current_slot The slot to insert
     * @param insert_pair The pair to insert
     * @param key_equal The binary callable used to compare two keys for
     * equality
     * @param expected_key The expected value of the key in the target slot
     * @return An insert result from the `insert_resullt` enumeration.
     */
    template <typename KeyEqual>
    __device__ insert_result back_to_back_cas(iterator current_slot,
                                              value_type const& insert_pair,
                                              KeyEqual key_equal,
                                              Key expected_key) noexcept;

    /**
     * @brief Inserts the specified key/value pair with a CAS of the key and a dependent write of
     * the value.
     *
     * @tparam KeyEqual Binary callable type
     * @param current_slot The slot to insert
     * @param insert_pair The pair to insert
     * @param key_equal The binary callable used to compare two keys for
     * equality
     * @param expected_key The expected value of the key in the target slot
     * @return An insert result from the `insert_resullt` enumeration.
     */
    template <typename KeyEqual>
    __device__ insert_result cas_dependent_write(iterator current_slot,
                                                 value_type const& insert_pair,
                                                 KeyEqual key_equal,
                                                 Key expected_key) noexcept;

   public:
    /**
     * @brief Given a slot pointer `slots`, initializes the first `capacity` slots with the given
     * sentinel values and returns a `device_mutable_view` object of those slots.
     *
     * @tparam CG The type of the cooperative thread group
     *
     * @param g The cooperative thread group used to copy the slots
     * @param slots Pointer to the hash map slots
     * @param capacity The total number of slots in the map
     * @param empty_key_sentinel The reserved value for keys to represent empty slots
     * @param empty_value_sentinel The reserved value for mapped values to represent empty slots
     * @return A device_mutable_view object based on the given parameters
     */
    template <typename CG>
    __device__ static device_mutable_view make_from_uninitialized_slots(
      CG g,
      pair_atomic_type* slots,
      std::size_t capacity,
      empty_key<Key> empty_key_sentinel,
      empty_value<Value> empty_value_sentinel) noexcept
    {
      device_view_base::initialize_slots(
        g, slots, capacity, empty_key_sentinel.value, empty_value_sentinel.value);
      return device_mutable_view{slots,
                                 capacity,
                                 empty_key_sentinel,
                                 empty_value_sentinel,
                                 erased_key<Key>{empty_key_sentinel.value}};
    }

    /**
     * @brief Given a slot pointer `slots`, initializes the first `capacity` slots with the given
     * sentinel values and returns a `device_mutable_view` object of those slots.
     *
     * @tparam CG The type of the cooperative thread group
     *
     * @param g The cooperative thread group used to copy the slots
     * @param slots Pointer to the hash map slots
     * @param capacity The total number of slots in the map
     * @param empty_key_sentinel The reserved value for keys to represent empty slots
     * @param empty_value_sentinel The reserved value for mapped values to represent empty slots
     * @param erased_key_sentinel The reserved value for keys to represent erased slots
     * @return A device_mutable_view object based on the given parameters
     */
    template <typename CG>
    __device__ static device_mutable_view make_from_uninitialized_slots(
      CG g,
      pair_atomic_type* slots,
      std::size_t capacity,
      empty_key<Key> empty_key_sentinel,
      empty_value<Value> empty_value_sentinel,
      erased_key<Key> erased_key_sentinel) noexcept
    {
      device_view_base::initialize_slots(
        g, slots, capacity, empty_key_sentinel, empty_value_sentinel);
      return device_mutable_view{
        slots, capacity, empty_key_sentinel, empty_value_sentinel, erased_key_sentinel};
    }

    /**
     * @brief Inserts the specified key/value pair into the map.
     *
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     * @param insert_pair The pair to insert
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys for
     * equality
     * @return `true` if the insert was successful, `false` otherwise.
     */
    template <typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ bool insert(value_type const& insert_pair,
                           Hash hash          = Hash{},
                           KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Inserts the specified key/value pair into the map.
     *
     * Returns a pair consisting of an iterator to the inserted element (or to
     * the element that prevented the insertion) and a `bool` denoting whether
     * the insertion took place.
     *
     * Note: In order to guarantee the validity of the returned iterator,
     * `insert_and_find` may be less efficient than `insert` in some situations.
     * Prefer using `insert` unless the returned iterator is required.
     *
     * Note: `insert_and_find` may only be used concurrently with `insert`,
     * `find`, and `erase` when `supports_concurrent_insert_find()` returns
     * true.
     *
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     *
     * @param insert_pair The pair to insert
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys for
     * equality
     * @return a pair consisting of an iterator to the element and a bool,
     * either `true` if the insert was successful, `false` otherwise.
     */
    template <typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ cuda::std::pair<iterator, bool> insert_and_find(
      value_type const& insert_pair, Hash hash = Hash{}, KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Inserts the specified key/value pair into the map.
     *
     * Uses the CUDA Cooperative Groups API to to leverage multiple threads to
     * perform a single insert. This provides a significant boost in throughput
     * compared to the non Cooperative Group `insert` at moderate to high load
     * factors.
     *
     * @tparam CG Cooperative Group type
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     *
     * @param g The Cooperative Group that performs the insert
     * @param insert_pair The pair to insert
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys for
     * equality
     * @return `true` if the insert was successful, `false` otherwise.
     */
    template <typename CG,
              typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ bool insert(CG g,
                           value_type const& insert_pair,
                           Hash hash          = Hash{},
                           KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Erases the specified key across the map.
     *
     * Behavior is undefined if `empty_key_sentinel_` equals to `erased_key_sentinel_`.
     *
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     *
     * @param k The key to be erased
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys for
     * equality
     * @return `true` if the erasure was successful, `false` otherwise.
     */
    template <typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ bool erase(key_type const& k,
                          Hash hash          = Hash{},
                          KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Erases the specified key across the map.
     *
     * Behavior is undefined if `empty_key_sentinel_` equals to `erased_key_sentinel_`.
     *
     * @tparam CG Cooperative Group type
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     *
     * @param g The Cooperative Group that performs the erasure
     * @param k The key to be erased
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys for
     * equality
     * @return `true` if the erasure was successful, `false` otherwise.
     */
    template <typename CG,
              typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ bool erase(CG g,
                          key_type const& k,
                          Hash hash          = Hash{},
                          KeyEqual key_equal = KeyEqual{}) noexcept;

  };  // class device mutable view

  /**
   * @brief Non-owning view-type that may be used in device code to
   * perform singular find and contains operations for the map.
   *
   * `device_view` is trivially-copyable and is intended to be passed by
   * value.
   *
   */
  class device_view : public device_view_base {
   public:
    using value_type  = typename device_view_base::value_type;   ///< Type of key/value pairs
    using key_type    = typename device_view_base::key_type;     ///< Key type
    using mapped_type = typename device_view_base::mapped_type;  ///< Type of the mapped values
    using iterator =
      typename device_view_base::iterator;  ///< Type of the forward iterator to `value_type`
    using const_iterator =
      typename device_view_base::const_iterator;  ///< Type of the forward iterator to `const
                                                  ///< value_type`
    using slot_type = typename device_view_base::slot_type;  ///< Type of hash map slots

    /**
     * @brief Construct a view of the first `capacity` slots of the
     * slots array pointed to by `slots`.
     *
     * @param slots Pointer to beginning of initialized slots array
     * @param capacity The number of slots viewed by this object
     * @param empty_key_sentinel The reserved value for keys to represent empty slots
     * @param empty_value_sentinel The reserved value for mapped values to represent empty slots
     */
    __host__ __device__ device_view(pair_atomic_type* slots,
                                    std::size_t capacity,
                                    empty_key<Key> empty_key_sentinel,
                                    empty_value<Value> empty_value_sentinel) noexcept
      : device_view_base{slots, capacity, empty_key_sentinel, empty_value_sentinel}
    {
    }

    /**
     * @brief Construct a view of the first `capacity` slots of the
     * slots array pointed to by `slots`.
     *
     * @param slots Pointer to beginning of initialized slots array
     * @param capacity The number of slots viewed by this object
     * @param empty_key_sentinel The reserved value for keys to represent empty slots
     * @param empty_value_sentinel The reserved value for mapped values to represent empty slots
     * @param erased_key_sentinel The reserved value for keys to represent erased slots
     */
    __host__ __device__ device_view(pair_atomic_type* slots,
                                    std::size_t capacity,
                                    empty_key<Key> empty_key_sentinel,
                                    empty_value<Value> empty_value_sentinel,
                                    erased_key<Key> erased_key_sentinel) noexcept
      : device_view_base{
          slots, capacity, empty_key_sentinel, empty_value_sentinel, erased_key_sentinel}
    {
    }

    /**
     * @brief Construct a `device_view` from a `device_mutable_view` object
     *
     * @param mutable_map object of type `device_mutable_view`
     */
    __host__ __device__ explicit device_view(device_mutable_view mutable_map)
      : device_view_base{mutable_map.get_slots(),
                         mutable_map.get_capacity(),
                         empty_key<Key>{mutable_map.get_empty_key_sentinel()},
                         empty_value<Value>{mutable_map.get_empty_value_sentinel()},
                         erased_key<Key>{mutable_map.get_erased_key_sentinel()}}
    {
    }

    /**
     * @brief Makes a copy of given `device_view` using non-owned memory.
     *
     * This function is intended to be used to create shared memory copies of small static maps,
     * although global memory can be used as well.
     *
     * Example:
     * @code{.cpp}
     * template <typename MapType, int CAPACITY>
     * __global__ void use_device_view(const typename MapType::device_view device_view,
     *                                 map_key_t const* const keys_to_search,
     *                                 map_value_t* const values_found,
     *                                 const size_t number_of_elements)
     * {
     *     const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
     *
     *     __shared__ typename MapType::pair_atomic_type sm_buffer[CAPACITY];
     *
     *     auto g = cg::this_thread_block();
     *
     *     const map_t::device_view sm_static_map = device_view.make_copy(g,
     *                                                                    sm_buffer);
     *
     *     for (size_t i = g.thread_rank(); i < number_of_elements; i += g.size())
     *     {
     *         values_found[i] = sm_static_map.find(keys_to_search[i])->second;
     *     }
     * }
     * @endcode
     *
     * @tparam CG The type of the cooperative thread group
     * @param g The cooperative thread group used to copy the slots
     * @param source_device_view `device_view` to copy from
     * @param memory_to_use Array large enough to support `capacity` elements. Object does not take
     * the ownership of the memory
     * @return Copy of passed `device_view`
     */
    template <typename CG>
    __device__ static device_view make_copy(CG g,
                                            pair_atomic_type* const memory_to_use,
                                            device_view source_device_view) noexcept
    {
#if defined(CUCO_HAS_CUDA_BARRIER)
      __shared__ cuda::barrier<cuda::thread_scope::thread_scope_block> barrier;
      if (g.thread_rank() == 0) { init(&barrier, g.size()); }
      g.sync();

      cuda::memcpy_async(g,
                         memory_to_use,
                         source_device_view.get_slots(),
                         sizeof(pair_atomic_type) * source_device_view.get_capacity(),
                         barrier);

      barrier.arrive_and_wait();
#else
      pair_atomic_type const* const slots_ptr = source_device_view.get_slots();
      for (std::size_t i = g.thread_rank(); i < source_device_view.get_capacity(); i += g.size()) {
        new (&memory_to_use[i].first)
          atomic_key_type{slots_ptr[i].first.load(cuda::memory_order_relaxed)};
        new (&memory_to_use[i].second)
          atomic_mapped_type{slots_ptr[i].second.load(cuda::memory_order_relaxed)};
      }
      g.sync();
#endif

      return device_view(memory_to_use,
                         source_device_view.get_capacity(),
                         empty_key<Key>{source_device_view.get_empty_key_sentinel()},
                         empty_value<Value>{source_device_view.get_empty_value_sentinel()},
                         erased_key<Key>{source_device_view.get_erased_key_sentinel()});
    }

    /**
     * @brief Finds the value corresponding to the key `k`.
     *
     * Returns an iterator to the pair whose key is equivalent to `k`.
     * If no such pair exists, returns `end()`.
     *
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     * @param k The key to search for
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return An iterator to the position at which the key/value pair
     * containing `k` was inserted
     */
    template <typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ iterator find(Key const& k,
                             Hash hash          = Hash{},
                             KeyEqual key_equal = KeyEqual{}) noexcept;

    /** @brief Finds the value corresponding to the key `k`.
     *
     * Returns a const_iterator to the pair whose key is equivalent to `k`.
     * If no such pair exists, returns `end()`.
     *
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     * @param k The key to search for
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return An iterator to the position at which the key/value pair
     * containing `k` was inserted
     */
    template <typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ const_iterator find(Key const& k,
                                   Hash hash          = Hash{},
                                   KeyEqual key_equal = KeyEqual{}) const noexcept;

    /**
     * @brief Finds the value corresponding to the key `k`.
     *
     * Returns an iterator to the pair whose key is equivalent to `k`.
     * If no such pair exists, returns `end()`. Uses the CUDA Cooperative Groups API to
     * to leverage multiple threads to perform a single find. This provides a
     * significant boost in throughput compared to the non Cooperative Group
     * `find` at moderate to high load factors.
     *
     * @tparam CG Cooperative Group type
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     * @param g The Cooperative Group used to perform the find
     * @param k The key to search for
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return An iterator to the position at which the key/value pair
     * containing `k` was inserted
     */
    template <typename CG,
              typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ iterator
    find(CG g, Key const& k, Hash hash = Hash{}, KeyEqual key_equal = KeyEqual{}) noexcept;

    /**
     * @brief Finds the value corresponding to the key `k`.
     *
     * Returns a const_iterator to the pair whose key is equivalent to `k`.
     * If no such pair exists, returns `end()`. Uses the CUDA Cooperative Groups API to
     * to leverage multiple threads to perform a single find. This provides a
     * significant boost in throughput compared to the non Cooperative Group
     * `find` at moderate to high load factors.
     *
     * @tparam CG Cooperative Group type
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     * @param g The Cooperative Group used to perform the find
     * @param k The key to search for
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return An iterator to the position at which the key/value pair
     * containing `k` was inserted
     */
    template <typename CG,
              typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ const_iterator
    find(CG g, Key const& k, Hash hash = Hash{}, KeyEqual key_equal = KeyEqual{}) const noexcept;

    /**
     * @brief Indicates whether the key `k` was inserted into the map.
     *
     * If the key `k` was inserted into the map, find returns
     * true. Otherwise, it returns false.
     *
     * Hash should be callable with both ProbeKey and Key type. `std::invoke_result<KeyEqual,
     * ProbeKey, Key>` must be well-formed.
     *
     * If `key_equal(probe_key, slot_key)` returns true, `hash(probe_key) == hash(slot_key)` must
     * also be true.
     *
     * @tparam ProbeKey Probe key type
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     *
     * @param k The key to search for
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return A boolean indicating whether the key/value pair
     * containing `k` was inserted
     */
    template <typename ProbeKey,
              typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ bool contains(ProbeKey const& k,
                             Hash hash          = Hash{},
                             KeyEqual key_equal = KeyEqual{}) const noexcept;

    /**
     * @brief Indicates whether the key `k` was inserted into the map.
     *
     * If the key `k` was inserted into the map, find returns true. Otherwise, it returns false.
     * Uses the CUDA Cooperative Groups API to to leverage multiple threads to perform a single
     * contains operation. This provides a significant boost in throughput compared to the non
     * Cooperative Group `contains` at moderate to high load factors.
     *
     * Hash should be callable with both ProbeKey and Key type. `std::invoke_result<KeyEqual,
     * ProbeKey, Key>` must be well-formed.
     *
     * If `key_equal(probe_key, slot_key)` returns true, `hash(probe_key) == hash(slot_key)` must
     * also be true.
     *
     * @tparam CG Cooperative Group type
     * @tparam ProbeKey Probe key type
     * @tparam Hash Unary callable type
     * @tparam KeyEqual Binary callable type
     *
     * @param g The Cooperative Group used to perform the contains operation
     * @param k The key to search for
     * @param hash The unary callable used to hash the key
     * @param key_equal The binary callable used to compare two keys
     * for equality
     * @return A boolean indicating whether the key/value pair
     * containing `k` was inserted
     */
    template <typename CG,
              typename ProbeKey,
              typename Hash     = cuco::default_hash_function<key_type>,
              typename KeyEqual = cuda::std::equal_to<key_type>>
    __device__ cuda::std::enable_if_t<std::is_invocable_v<KeyEqual, ProbeKey, Key>, bool> contains(
      CG g, ProbeKey const& k, Hash hash = Hash{}, KeyEqual key_equal = KeyEqual{}) const noexcept;
  };  // class device_view

  /**
   * @brief Gets the maximum number of elements the hash map can hold.
   *
   * @return The maximum number of elements the hash map can hold
   */
  std::size_t get_capacity() const noexcept { return capacity_; }

  /**
   * @brief Gets the number of elements in the hash map.
   *
   * @return The number of elements in the map
   */
  std::size_t get_size() const noexcept { return size_; }

  /**
   * @brief Gets the load factor of the hash map.
   *
   * @return The load factor of the hash map
   */
  float get_load_factor() const noexcept { return static_cast<float>(size_) / capacity_; }

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
   * @brief Gets the sentinel value used to represent an erased value slot.
   *
   * @return The sentinel value used to represent an erased value slot
   */
  Key get_erased_key_sentinel() const noexcept { return erased_key_sentinel_; }

  /**
   * @brief Constructs a device_view object based on the members of the `static_map` object.
   *
   * @return A device_view object based on the members of the `static_map` object
   */
  device_view get_device_view() const noexcept
  {
    return device_view(slots_,
                       capacity_,
                       empty_key<Key>{empty_key_sentinel_},
                       empty_value<Value>{empty_value_sentinel_},
                       erased_key<Key>{erased_key_sentinel_});
  }

  /**
   * @brief Constructs a device_mutable_view object based on the members of the `static_map` object
   *
   * @return A device_mutable_view object based on the members of the `static_map` object
   */
  device_mutable_view get_device_mutable_view() const noexcept
  {
    return device_mutable_view(slots_,
                               capacity_,
                               empty_key<Key>{empty_key_sentinel_},
                               empty_value<Value>{empty_value_sentinel_},
                               erased_key<Key>{erased_key_sentinel_});
  }

 private:
  pair_atomic_type* slots_{};                   ///< Pointer to flat slots storage
  std::size_t capacity_{};                      ///< Total number of slots
  std::size_t size_{};                          ///< Number of keys in map
  Key empty_key_sentinel_{};                    ///< Key value that represents an empty slot
  Value empty_value_sentinel_{};                ///< Initial value of empty slot
  Key erased_key_sentinel_{};                   ///< Key value that represents an erased slot
  atomic_ctr_type* num_successes_{};            ///< Number of successfully inserted keys on insert
  slot_allocator_type slot_allocator_{};        ///< Allocator used to allocate slots
  counter_allocator_type counter_allocator_{};  ///< Allocator used to allocate `num_successes_`
};
}  // namespace legacy
}  // namespace cuco

#include <cuco/detail/static_map.inl>
#include <cuco/detail/static_map/static_map.inl>
