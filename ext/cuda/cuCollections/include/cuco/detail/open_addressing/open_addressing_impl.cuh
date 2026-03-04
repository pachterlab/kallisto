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

#include <cuco/detail/__config>
#include <cuco/detail/open_addressing/functors.cuh>
#include <cuco/detail/open_addressing/kernels.cuh>
#include <cuco/detail/storage/counter_storage.cuh>
#include <cuco/detail/utility/cuda.hpp>
#include <cuco/detail/utils.hpp>
#include <cuco/extent.cuh>
#include <cuco/operator.hpp>
#include <cuco/probing_scheme.cuh>
#include <cuco/storage.cuh>
#include <cuco/utility/traits.hpp>

#include <cub/device/device_for.cuh>
#include <cub/device/device_select.cuh>
#include <cuda/atomic>
#include <cuda/std/functional>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <cmath>
#include <cstdint>

namespace cuco {
namespace detail {
/**
 * @brief An open addressing impl class.
 *
 * @note This class should NOT be used directly.
 *
 * @throw If the size of the given key type is larger than 8 bytes
 * @throw If the size of the given slot type is larger than 16 bytes
 * @throw If the given key type doesn't have unique object representations, i.e.,
 * `cuco::bitwise_comparable_v<Key> == false`
 * @throw If the probing scheme type is not inherited from `cuco::detail::probing_scheme_base`
 *
 * @tparam Key Type used for keys. Requires `cuco::is_bitwise_comparable_v<Key>`
 * @tparam Value Type used for storage values.
 * @tparam Extent Data structure size type
 * @tparam Scope The scope in which operations will be performed by individual threads.
 * @tparam KeyEqual Binary callable type used to compare two keys for equality
 * @tparam ProbingScheme Probing scheme (see `include/cuco/probing_scheme.cuh` for choices)
 * @tparam Allocator Type of allocator used for device storage
 * @tparam Storage Slot bucket storage type
 */
template <class Key,
          class Value,
          class Extent,
          cuda::thread_scope Scope,
          class KeyEqual,
          class ProbingScheme,
          class Allocator,
          class Storage>
class open_addressing_impl {
  static_assert(sizeof(Key) <= 8, "Container does not support key types larger than 8 bytes.");

  static_assert(sizeof(Value) <= 16, "Container does not support slot types larger than 16 bytes.");

  static_assert(
    cuco::is_bitwise_comparable_v<Key>,
    "Key type must have unique object representations or have been explicitly declared as safe for "
    "bitwise comparison via specialization of cuco::is_bitwise_comparable_v<Key>.");

  static_assert(cuda::std::is_base_of_v<cuco::detail::probing_scheme_base<ProbingScheme::cg_size>,
                                        ProbingScheme>,
                "ProbingScheme must inherit from cuco::detail::probing_scheme_base");

  /// Determines if the container is a key/value or key-only store
  static constexpr auto has_payload = not cuda::std::is_same_v<Key, Value>;

 public:
  static constexpr auto cg_size      = ProbingScheme::cg_size;  ///< CG size used for probing
  static constexpr auto bucket_size  = Storage::bucket_size;    ///< Bucket size used for probing
  static constexpr auto thread_scope = Scope;                   ///< CUDA thread scope

  using key_type            = Key;            ///< Key type
  using value_type          = Value;          ///< The storage value type, NOT payload type
  using probing_scheme_type = ProbingScheme;  ///< Probe scheme type
  using hasher              = typename probing_scheme_type::hasher;  ///< Hash function type
  /// Extent type
  using extent_type =
    decltype(make_valid_extent<probing_scheme_type, Storage>(std::declval<Extent>()));
  using size_type = typename extent_type::value_type;  ///< Size type
  using key_equal = KeyEqual;                          ///< Key equality comparator type
  using storage_type =
    detail::storage<Storage, value_type, extent_type, Allocator>;  ///< Storage type
  using allocator_type = typename storage_type::allocator_type;    ///< Allocator type

  using storage_ref_type = typename storage_type::ref_type;  ///< Non-owning bucket storage ref type

  /**
   * @brief Constructs a statically-sized open addressing data structure with the specified initial
   * capacity, sentinel values and CUDA stream.
   *
   * @note The actual capacity depends on the given `capacity`, the probing scheme, CG size, and the
   * bucket size and it is computed via the `make_valid_extent` factory. Insert operations will not
   * automatically grow the container. Attempting to insert more unique keys than the capacity of
   * the container results in undefined behavior.
   * @note Any `*_sentinel`s are reserved and behavior is undefined when attempting to insert
   * this sentinel value.
   * @note If a non-default CUDA stream is provided, the caller is responsible for synchronizing the
   * stream before the object is first used.
   *
   * @param capacity The requested lower-bound size
   * @param empty_slot_sentinel The reserved slot value for empty slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the data structure
   */
  constexpr open_addressing_impl(Extent capacity,
                                 Value empty_slot_sentinel,
                                 KeyEqual const& pred,
                                 ProbingScheme const& probing_scheme,
                                 Allocator const& alloc,
                                 cuda::stream_ref stream)
    : empty_slot_sentinel_{empty_slot_sentinel},
      erased_key_sentinel_{this->extract_key(empty_slot_sentinel)},
      predicate_{pred},
      probing_scheme_{probing_scheme},
      storage_{make_valid_extent<probing_scheme_type, Storage>(capacity), alloc, stream}
  {
    this->clear_async(stream);
  }

  /**
   * @brief Constructs a statically-sized open addressing data structure with the number of elements
   * to insert `n`, the desired load factor, etc.
   *
   * @note This constructor helps users create a data structure based on the number of elements to
   * insert and the desired load factor without manually computing the desired capacity. The actual
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
   * @param empty_slot_sentinel The reserved slot value for empty slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the data structure
   */
  constexpr open_addressing_impl(Extent n,
                                 double desired_load_factor,
                                 Value empty_slot_sentinel,
                                 KeyEqual const& pred,
                                 ProbingScheme const& probing_scheme,
                                 Allocator const& alloc,
                                 cuda::stream_ref stream)
    : empty_slot_sentinel_{empty_slot_sentinel},
      erased_key_sentinel_{this->extract_key(empty_slot_sentinel)},
      predicate_{pred},
      probing_scheme_{probing_scheme},
      storage_{
        make_valid_extent<probing_scheme_type, Storage>(n, desired_load_factor), alloc, stream}
  {
    this->clear_async(stream);
  }

  /**
   * @brief Constructs a statically-sized open addressing data structure with the specified initial
   * capacity, sentinel values and CUDA stream.
   *
   * @note The actual capacity depends on the given `capacity`, the probing scheme, CG size, and the
   * bucket size and it is computed via the `make_valid_extent` factory. Insert operations will not
   * automatically grow the container. Attempting to insert more unique keys than the capacity of
   * the container results in undefined behavior.
   * @note Any `*_sentinel`s are reserved and behavior is undefined when attempting to insert
   * this sentinel value.
   * @note If a non-default CUDA stream is provided, the caller is responsible for synchronizing the
   * stream before the object is first used.
   *
   * @param capacity The requested lower-bound size
   * @param empty_slot_sentinel The reserved slot value for empty slots
   * @param erased_key_sentinel The reserved key value for erased slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the data structure
   */
  constexpr open_addressing_impl(Extent capacity,
                                 Value empty_slot_sentinel,
                                 Key erased_key_sentinel,
                                 KeyEqual const& pred,
                                 ProbingScheme const& probing_scheme,
                                 Allocator const& alloc,
                                 cuda::stream_ref stream)
    : empty_slot_sentinel_{empty_slot_sentinel},
      erased_key_sentinel_{erased_key_sentinel},
      predicate_{pred},
      probing_scheme_{probing_scheme},
      storage_{make_valid_extent<probing_scheme_type, Storage>(capacity), alloc, stream}
  {
    CUCO_EXPECTS(this->empty_key_sentinel() != this->erased_key_sentinel(),
                 "The empty key sentinel and erased key sentinel cannot be the same value.",
                 std::logic_error);

    this->clear_async(stream);
  }

  /**
   * @brief Erases all elements from the container. After this call, `size()` returns zero.
   * Invalidates any references, pointers, or iterators referring to contained elements.
   *
   * @param stream CUDA stream this operation is executed in
   */
  void clear(cuda::stream_ref stream) { storage_.initialize(empty_slot_sentinel_, stream); }

  /**
   * @brief Asynchronously erases all elements from the container. After this call, `size()` returns
   * zero. Invalidates any references, pointers, or iterators referring to contained elements.
   *
   * @param stream CUDA stream this operation is executed in
   */
  void clear_async(cuda::stream_ref stream) noexcept
  {
    storage_.initialize_async(empty_slot_sentinel_, stream);
  }

  /**
   * @brief Inserts all keys in the range `[first, last)` and returns the number of successful
   * insertions.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `insert_async`.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * open_addressing_impl::value_type></tt> is `true`
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream CUDA stream used for insert
   *
   * @return Number of successfully inserted keys
   */
  template <typename InputIt, typename Ref>
  size_type insert(InputIt first, InputIt last, Ref container_ref, cuda::stream_ref stream)
  {
    auto const always_true = thrust::constant_iterator<bool>{true};
    return this->insert_if(first, last, always_true, cuda::std::identity{}, container_ref, stream);
  }

  /**
   * @brief Asynchronously inserts all keys in the range `[first, last)`.
   *
   * @tparam InputIt Device accessible random access input iterator where
   * <tt>std::is_convertible<std::iterator_traits<InputIt>::value_type,
   * open_addressing_impl::value_type></tt> is `true`
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream CUDA stream used for insert
   */
  template <typename InputIt, typename Ref>
  void insert_async(InputIt first,
                    InputIt last,
                    Ref container_ref,
                    cuda::stream_ref stream) noexcept
  {
    auto const always_true = thrust::constant_iterator<bool>{true};
    this->insert_if_async(first, last, always_true, cuda::std::identity{}, container_ref, stream);
  }

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
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream CUDA stream used for the operation
   *
   * @return Number of successfully inserted keys
   */
  template <typename InputIt, typename StencilIt, typename Predicate, typename Ref>
  size_type insert_if(InputIt first,
                      InputIt last,
                      StencilIt stencil,
                      Predicate pred,
                      Ref container_ref,
                      cuda::stream_ref stream)
  {
    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return 0; }

    auto counter =
      detail::counter_storage<size_type, thread_scope, allocator_type>{this->allocator(), stream};
    counter.reset(stream);

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::insert_if_n<cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, stencil, pred, counter.data(), container_ref);

    return counter.load_to_host(stream);
  }

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
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream CUDA stream used for the operation
   */
  template <typename InputIt, typename StencilIt, typename Predicate, typename Ref>
  void insert_if_async(InputIt first,
                       InputIt last,
                       StencilIt stencil,
                       Predicate pred,
                       Ref container_ref,
                       cuda::stream_ref stream) noexcept
  {
    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return; }

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::insert_if_n<cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, stencil, pred, container_ref);
  }

  /**
   * @brief Asynchronously inserts all elements in the range `[first, last)`.
   *
   * @note: For a given element `*(first + i)`, if the container doesn't already contain an element
   * with an equivalent key, inserts the element at a location pointed by `iter` and writes
   * `iter` to `found_begin + i` and writes `true` to `inserted_begin + i`. Otherwise, finds the
   * location of the equivalent element, `iter` and writes `iter` to `found_begin + i` and writes
   * `false` to `inserted_begin + i`.
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the `value_type` of the data structure
   * @tparam FoundIt Device accessible random access output iterator whose `value_type`
   * is constructible from `map::iterator` type
   * @tparam InsertedIt Device accessible random access output iterator whose `value_type`
   * is constructible from `bool`
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of input elements
   * @param last End of the sequence of elements
   * @param found_begin Beginning of the sequence of elements found for each key
   * @param inserted_begin Beginning of the sequence of booleans for the presence of each key
   * @param ref Non-owning container device ref used to access the slot storage
   * @param stream CUDA stream used for the operation
   */
  template <typename InputIt, typename FoundIt, typename InsertedIt, typename Ref>
  void insert_and_find_async(InputIt first,
                             InputIt last,
                             FoundIt found_begin,
                             InsertedIt inserted_begin,
                             Ref container_ref,
                             cuda::stream_ref stream) noexcept
  {
    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return; }

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::insert_and_find<cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, found_begin, inserted_begin, container_ref);
  }

  /**
   * @brief Asynchronously erases keys in the range `[first, last)`.
   *
   * @note For each key `k` in `[first, last)`, if contains(k) returns true, removes `k` and it's
   * associated value from the container. Else, no effect.
   *
   * @note Side-effects:
   *  - `contains(k) == false`
   *  - `find(k) == end()`
   *  - `insert({k,v}) == true`
   *  - `size()` is reduced by the total number of erased keys
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream Stream used for executing the kernels
   *
   * @throw std::runtime_error if a unique erased key sentinel value was not
   * provided at construction
   */
  template <typename InputIt, typename Ref>
  void erase_async(InputIt first,
                   InputIt last,
                   Ref container_ref,
                   cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}})
  {
    CUCO_EXPECTS(this->empty_key_sentinel() != this->erased_key_sentinel(),
                 "The empty key sentinel and erased key sentinel cannot be the same value.",
                 std::logic_error);

    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return; }

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::erase<cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, container_ref);
  }

  /**
   * @brief Asynchronously indicates whether the keys in the range `[first, last)` are contained in
   * the container.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename OutputIt, typename Ref>
  void contains_async(InputIt first,
                      InputIt last,
                      OutputIt output_begin,
                      Ref container_ref,
                      cuda::stream_ref stream) const noexcept
  {
    auto const always_true = thrust::constant_iterator<bool>{true};
    this->contains_if_async(
      first, last, always_true, cuda::std::identity{}, output_begin, container_ref, stream);
  }

  /**
   * @brief Asynchronously indicates whether the keys in the range `[first, last)` are contained in
   * the container if `pred` of the corresponding stencil returns true.
   *
   * @note If `pred( *(stencil + i) )` is true, stores `true` or `false` to `(output_begin + i)`
   * indicating if the key `*(first + i)` is present int the container. If `pred( *(stencil + i) )`
   * is false, stores false to `(output_begin + i)`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam OutputIt Device accessible output iterator assignable from `bool`
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param output_begin Beginning of the sequence of booleans for the presence of each key
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt,
            typename StencilIt,
            typename Predicate,
            typename OutputIt,
            typename Ref>
  void contains_if_async(InputIt first,
                         InputIt last,
                         StencilIt stencil,
                         Predicate pred,
                         OutputIt output_begin,
                         Ref container_ref,
                         cuda::stream_ref stream) const noexcept
  {
    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return; }

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::contains_if_n<cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, stencil, pred, output_begin, container_ref);
  }

  /**
   * @brief For all keys in the range `[first, last)`, asynchronously finds
   * a match with its key equivalent to the query key.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param output_begin Beginning of the sequence of matches retrieved for each key
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt, typename OutputIt, typename Ref>
  void find_async(InputIt first,
                  InputIt last,
                  OutputIt output_begin,
                  Ref container_ref,
                  cuda::stream_ref stream) const noexcept
  {
    auto const always_true = thrust::constant_iterator<bool>{true};

    this->find_if_async(
      first, last, always_true, cuda::std::identity{}, output_begin, container_ref, stream);
  }

  /**
   * @brief For all keys in the range `[first, last)`, asynchronously finds
   * a match with its key equivalent to the query key.
   *
   * @note If `pred( *(stencil + i) )` is true, stores the payload of the
   * matched key or the `empty_value_sentinel` to `(output_begin + i)`. If `pred( *(stencil + i) )`
   * is false, stores `empty_value_sentinel` to `(output_begin + i)`.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam StencilIt Device accessible random access iterator whose value_type is
   * convertible to Predicate's argument type
   * @tparam Predicate Unary predicate callable whose return type must be convertible to `bool` and
   * argument type is convertible from <tt>std::iterator_traits<StencilIt>::value_type</tt>
   * @tparam OutputIt Device accessible output iterator
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param stencil Beginning of the stencil sequence
   * @param pred Predicate to test on every element in the range `[stencil, stencil +
   * std::distance(first, last))`
   * @param output_begin Beginning of the sequence of matches retrieved for each key
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt,
            typename StencilIt,
            typename Predicate,
            typename OutputIt,
            typename Ref>
  void find_if_async(InputIt first,
                     InputIt last,
                     StencilIt stencil,
                     Predicate pred,
                     OutputIt output_begin,
                     Ref container_ref,
                     cuda::stream_ref stream) const noexcept
  {
    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return; }

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::find_if_n<cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, stencil, pred, output_begin, container_ref);
  }

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
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the input sequence of keys
   * @param last End of the input sequence of keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param container_ref Non-owning device reference to the container
   * @param stream CUDA stream this operation is executed in
   *
   * @return Iterator pair indicating the the end of the output sequences
   */
  template <class InputProbeIt, class OutputProbeIt, class OutputMatchIt, class Ref>
  std::pair<OutputProbeIt, OutputMatchIt> retrieve(InputProbeIt first,
                                                   InputProbeIt last,
                                                   OutputProbeIt output_probe,
                                                   OutputMatchIt output_match,
                                                   Ref container_ref,
                                                   cuda::stream_ref stream) const
  {
    auto constexpr is_outer = false;
    return this->retrieve_impl<is_outer>(
      first, last, output_probe, output_match, container_ref, stream);
  }

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
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the `InputProbeIt`'s `value_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the input sequence of keys
   * @param last End of the input sequence of keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param container_ref Non-owning device reference to the container
   * @param stream CUDA stream this operation is executed in
   *
   * @return Iterator pair indicating the the end of the output sequences
   */
  template <class InputProbeIt, class OutputProbeIt, class OutputMatchIt, class Ref>
  std::pair<OutputProbeIt, OutputMatchIt> retrieve_outer(InputProbeIt first,
                                                         InputProbeIt last,
                                                         OutputProbeIt output_probe,
                                                         OutputMatchIt output_match,
                                                         Ref container_ref,
                                                         cuda::stream_ref stream) const
  {
    auto constexpr is_outer = true;
    return this->retrieve_impl<is_outer>(
      first, last, output_probe, output_match, container_ref, stream);
  }

  /**
   * @brief Counts the occurrences of keys in `[first, last)` contained in the container
   *
   * @tparam Input Device accessible input iterator
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param stream CUDA stream used for count
   *
   * @return The sum of total occurrences of all keys in `[first, last)`
   */
  template <typename InputIt, typename Ref>
  [[nodiscard]] size_type count(InputIt first,
                                InputIt last,
                                Ref container_ref,
                                cuda::stream_ref stream) const
  {
    auto constexpr is_outer = false;
    return this->count<is_outer>(first, last, container_ref, stream);
  }

  /**
   * @brief Counts the occurrences of keys in `[first, last)` contained in the container
   *
   * @note If a given key has no matches, its occurrence is 1.
   *
   * @tparam Input Device accessible input iterator
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param stream CUDA stream used for count
   *
   * @return The sum of total occurrences of all keys in `[first, last)`
   */
  template <typename InputIt, typename Ref>
  [[nodiscard]] size_type count_outer(InputIt first,
                                      InputIt last,
                                      Ref container_ref,
                                      cuda::stream_ref stream) const noexcept
  {
    auto constexpr is_outer = true;
    return this->count<is_outer>(first, last, container_ref, stream);
  }

  /**
   * @brief Counts the number of occurrences of each query key in the container
   *
   * For each key in the input range `[first, last)`, this function computes the number of matching
   * elements in the container and writes the result to the corresponding position in the output
   * range starting at `output_begin`.
   *
   * @note The input and output ranges must be device-accessible and of the same length.
   * @note The behavior is undefined if the input and output ranges overlap.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param output_begin Beginning of the output range where per-key counts will be stored
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream CUDA stream used for this operation
   */
  template <typename InputIt, typename OutputIt, typename Ref>
  void count_each(InputIt first,
                  InputIt last,
                  OutputIt output_begin,
                  Ref container_ref,
                  cuda::stream_ref stream) const noexcept
  {
    auto constexpr is_outer = false;
    this->count_each<is_outer>(first, last, output_begin, container_ref, stream);
  }

  /**
   * @brief Counts the number of occurrences of each query key in the container with outer semantics
   *
   * For each key in the input range `[first, last)`, this function computes the number of matching
   * elements in the container and writes the result to the corresponding position in the output
   * range starting at `output_begin`.
   *
   * If a query key has no matches in the container, the result for that key will be 1 instead of 0.
   *
   * @note The input and output ranges must be device-accessible and of the same length.
   * @note The behavior is undefined if the input and output ranges overlap.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param output_begin Beginning of the output range where per-key counts will be stored
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream CUDA stream used for this operation
   */
  template <typename InputIt, typename OutputIt, typename Ref>
  void count_each_outer(InputIt first,
                        InputIt last,
                        OutputIt output_begin,
                        Ref container_ref,
                        cuda::stream_ref stream) const noexcept
  {
    auto constexpr is_outer = true;
    this->count_each<is_outer>(first, last, output_begin, container_ref, stream);
  }

  /**
   * @brief Retrieves all keys contained in the container.
   *
   * @note This API synchronizes the given stream.
   * @note The order in which keys are returned is implementation defined and not guaranteed to be
   * consistent between subsequent calls to `retrieve_all`.
   * @note Behavior is undefined if the range beginning at `output_begin` is smaller than the return
   * value of `size()`.
   *
   * @tparam OutputIt Device accessible random access output iterator whose `value_type` is
   * convertible from the container's `value_type`
   *
   * @param output_begin Beginning output iterator for keys
   * @param stream CUDA stream used for this operation
   *
   * @return Iterator indicating the end of the output
   */
  template <typename OutputIt>
  [[nodiscard]] OutputIt retrieve_all(OutputIt output_begin, cuda::stream_ref stream) const
  {
    using temp_allocator_type =
      typename std::allocator_traits<allocator_type>::template rebind_alloc<char>;

    cuco::detail::index_type constexpr stride = std::numeric_limits<std::int32_t>::max();

    cuco::detail::index_type h_num_out{0};
    auto temp_allocator = temp_allocator_type{this->allocator()};
    auto d_num_out =
      reinterpret_cast<size_type*>(temp_allocator.allocate(sizeof(size_type), stream));

    // TODO: PR #580 to be reverted once https://github.com/NVIDIA/cccl/issues/1422 is resolved
    for (cuco::detail::index_type offset = 0;
         offset < static_cast<cuco::detail::index_type>(this->capacity());
         offset += stride) {
      auto const num_items =
        std::min(static_cast<cuco::detail::index_type>(this->capacity()) - offset, stride);
      auto const begin = thrust::make_transform_iterator(
        thrust::counting_iterator{static_cast<size_type>(offset)},
        detail::open_addressing_ns::get_slot<has_payload, storage_ref_type>(this->storage_ref()));
      auto const is_filled = detail::open_addressing_ns::slot_is_filled<has_payload, key_type>{
        this->empty_key_sentinel(), this->erased_key_sentinel()};

      std::size_t temp_storage_bytes = 0;

      CUCO_CUDA_TRY(cub::DeviceSelect::If(nullptr,
                                          temp_storage_bytes,
                                          begin,
                                          output_begin + h_num_out,
                                          d_num_out,
                                          static_cast<std::int32_t>(num_items),
                                          is_filled,
                                          stream.get()));

      // Allocate temporary storage
      auto d_temp_storage = temp_allocator.allocate(temp_storage_bytes, stream);

      CUCO_CUDA_TRY(cub::DeviceSelect::If(d_temp_storage,
                                          temp_storage_bytes,
                                          begin,
                                          output_begin + h_num_out,
                                          d_num_out,
                                          static_cast<std::int32_t>(num_items),
                                          is_filled,
                                          stream.get()));

      size_type temp_count;
      CUCO_CUDA_TRY(cudaMemcpyAsync(
        &temp_count, d_num_out, sizeof(size_type), cudaMemcpyDeviceToHost, stream.get()));
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
      stream.sync();
#else
      stream.wait();
#endif
      h_num_out += temp_count;
      temp_allocator.deallocate(d_temp_storage, temp_storage_bytes, stream);
    }

    temp_allocator.deallocate(reinterpret_cast<char*>(d_num_out), sizeof(size_type), stream);

    return output_begin + h_num_out;
  }

  /**
   * @brief Asynchronously applies the given function object `callback_op` to the copy of every
   * filled slot in the container
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @tparam CallbackOp Type of unary callback function object
   *
   * @param callback_op Function to call on every filled slot in the container
   * @param stream CUDA stream used for this operation
   */
  template <typename CallbackOp>
  void for_each_async(CallbackOp&& callback_op, cuda::stream_ref stream) const
  {
    auto const is_filled = detail::open_addressing_ns::slot_is_filled<has_payload, key_type>{
      this->empty_key_sentinel(), this->erased_key_sentinel()};

    auto storage_ref = this->storage_ref();
    auto const op    = [callback_op, is_filled, storage_ref] __device__(auto const slot_or_bucket) {
      if (is_filled(slot_or_bucket)) { callback_op(slot_or_bucket); }
    };

    CUCO_CUDA_TRY(cub::DeviceFor::ForEachCopyN(
      storage_ref.data(), static_cast<size_type>(storage_ref.extent()), op, stream.get()));
  }

  /**
   * @brief For each key in the range [first, last), asynchronously applies the function object
   * `callback_op` to the copy of all corresponding matches found in the container.
   *
   * @note The return value of `callback_op`, if any, is ignored.
   *
   * @tparam InputIt Device accessible random access input iterator
   * @tparam CallbackOp Type of unary callback function object
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys
   * @param last End of the sequence of keys
   * @param callback_op Function to call on every match found in the container
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream CUDA stream used for this operation
   */
  template <typename InputIt, typename CallbackOp, typename Ref>
  void for_each_async(InputIt first,
                      InputIt last,
                      CallbackOp&& callback_op,
                      Ref container_ref,
                      cuda::stream_ref stream) const noexcept
  {
    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return; }

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::for_each_n<cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, std::forward<CallbackOp>(callback_op), container_ref);
  }

  /**
   * @brief Gets the number of elements in the container
   *
   * @note This function synchronizes the given stream.
   *
   * @param stream CUDA stream used to get the number of inserted elements
   *
   * @return The number of elements in the container
   */
  [[nodiscard]] size_type size(cuda::stream_ref stream) const
  {
    auto counter =
      detail::counter_storage<size_type, thread_scope, allocator_type>{this->allocator(), stream};
    counter.reset(stream);

    auto const grid_size = cuco::detail::grid_size(this->capacity());
    auto const is_filled = detail::open_addressing_ns::slot_is_filled<has_payload, key_type>{
      this->empty_key_sentinel(), this->erased_key_sentinel()};

    // TODO: custom kernel to be replaced by cub::DeviceReduce::Sum when cub version is bumped to
    // v2.1.0
    detail::open_addressing_ns::size<cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        storage_.ref(), is_filled, counter.data());

    return counter.load_to_host(stream);
  }

  /**
   * @brief Regenerates the container
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `rehash_async`.
   *
   * @tparam Container The container type this function operates on
   *
   * @param extent The container's new `bucket_extent` after this operation took place
   * @param container The container to be rehashed
   * @param stream CUDA stream used for this operation
   */
  template <typename Container>
  void rehash(Container const& container, cuda::stream_ref stream)
  {
    this->rehash_async(container, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
    stream.sync();
#else
    stream.wait();
#endif
  }

  /**
   * @brief Asynchronously reserves at least the specified number of slots and regenerates the
   * container
   *
   * @note Changes the number of buckets to a value that is not less than `extent`, then
   * rehashes the container, i.e. puts the elements into appropriate slots considering
   * that the total number of slots has changed.
   *
   * @note This function synchronizes the given stream. For asynchronous execution use
   * `rehash_async`.
   *
   * @note Behavior is undefined if the desired `extent` is insufficient to store all of the
   * contained elements.
   *
   * @note This function is not available if the container's `extent_type` is static.
   *
   * @tparam Container The container type this function operates on
   *
   * @param extent The container's new `bucket_extent` after this operation took place
   * @param container The container to be rehashed
   * @param stream CUDA stream used for this operation
   */
  template <typename Container>
  void rehash(extent_type extent, Container const& container, cuda::stream_ref stream)
  {
    this->rehash_async(extent, container, stream);
#if CCCL_MAJOR_VERSION > 3 || (CCCL_MAJOR_VERSION == 3 && CCCL_MINOR_VERSION >= 1)
    stream.sync();
#else
    stream.wait();
#endif
  }

  /**
   * @brief Asynchronously regenerates the container
   *
   * @tparam Container The container type this function operates on
   *
   * @param extent The container's new `bucket_extent` after this operation took place
   * @param container The container to be rehashed
   * @param stream CUDA stream used for this operation
   */
  template <typename Container>
  void rehash_async(Container const& container, cuda::stream_ref stream)
  {
    this->rehash_async(storage_.extent(), container, stream);
  }

  /**
   * @brief Asynchronously reserves at least the specified number of slots and regenerates the
   * container
   *
   * @note Changes the number of buckets to a value that is not less than `extent`, then
   * rehashes the container, i.e. puts the elements into appropriate slots considering
   * that the total number of slots has changed.
   *
   * @note Behavior is undefined if the desired `extent` is insufficient to store all of the
   * contained elements.
   *
   * @note This function is not available if the container's `extent_type` is static.
   *
   * @tparam Container The container type this function operates on
   *
   * @param extent The container's new `bucket_extent` after this operation took place
   * @param container The container to be rehashed
   * @param stream CUDA stream used for this operation
   */
  template <typename Container>
  void rehash_async(extent_type extent, Container const& container, cuda::stream_ref stream)
  {
    auto const old_storage = std::move(this->storage_);
    new (&storage_) storage_type{extent, this->allocator(), stream};
    this->clear_async(stream);

    auto const num_buckets = old_storage.num_buckets();
    if (num_buckets == 0) { return; }

    auto constexpr block_size = cuco::detail::default_block_size();
    auto constexpr stride     = cuco::detail::default_stride();
    auto const grid_size      = cuco::detail::grid_size(num_buckets, 1, stride, block_size);
    auto const is_filled      = detail::open_addressing_ns::slot_is_filled<has_payload, key_type>{
      this->empty_key_sentinel(), this->erased_key_sentinel()};

    detail::open_addressing_ns::rehash<block_size><<<grid_size, block_size, 0, stream.get()>>>(
      old_storage.ref(), container.ref(op::insert), is_filled);
  }

  /**
   * @brief Gets the maximum number of elements the container can hold.
   *
   * @return The maximum number of elements the container can hold
   */
  [[nodiscard]] constexpr auto capacity() const noexcept { return storage_.capacity(); }

  /**
   * @brief Gets a pointer to the underlying slot storage.
   *
   * @return Pointer to the underlying slot storage
   */
  [[nodiscard]] __host__ value_type* data() const { return storage_.data(); }

  /**
   * @brief Gets the sentinel value used to represent an empty key slot.
   *
   * @return The sentinel value used to represent an empty key slot
   */
  [[nodiscard]] constexpr key_type empty_key_sentinel() const noexcept
  {
    return this->extract_key(this->empty_slot_sentinel_);
  }

  /**
   * @brief Gets the sentinel value used to represent an erased key slot.
   *
   * @return The sentinel value used to represent an erased key slot
   */
  [[nodiscard]] constexpr key_type erased_key_sentinel() const noexcept
  {
    return erased_key_sentinel_;
  }

  /**
   * @brief Gets the key comparator.
   *
   * @return The comparator used to compare keys
   */
  [[nodiscard]] constexpr key_equal key_eq() const noexcept { return predicate_; }

  /**
   * @brief Gets the probing scheme.
   *
   * @return The probing scheme used for the container
   */
  [[nodiscard]] constexpr probing_scheme_type const& probing_scheme() const noexcept
  {
    return probing_scheme_;
  }

  /**
   * @brief Gets the function(s) used to hash keys
   *
   * @return The function(s) used to hash keys
   */
  [[nodiscard]] constexpr hasher hash_function() const noexcept
  {
    return this->probing_scheme().hash_function();
  }

  /**
   * @brief Gets the container allocator.
   *
   * @return The container allocator
   */
  [[nodiscard]] constexpr allocator_type allocator() const noexcept { return storage_.allocator(); }

  /**
   * @brief Gets the non-owning storage ref.
   *
   * @return The non-owning storage ref of the container
   */
  [[nodiscard]] constexpr storage_ref_type storage_ref() const noexcept { return storage_.ref(); }

 private:
  /**
   * @brief Counts the occurrences of keys in `[first, last)` contained in the container
   *
   * @note If `IsOuter` is `true`, the occurrence of a non-match key is 1. Else, it's 0.
   *
   * @tparam IsOuter Flag indicating whether it's an outer count or not
   * @tparam Input Device accessible input iterator
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param stream CUDA stream used for count
   *
   * @return The sum of total occurrences of all keys in `[first, last)`
   */
  template <bool IsOuter, typename InputIt, typename Ref>
  [[nodiscard]] size_type count(InputIt first,
                                InputIt last,
                                Ref container_ref,
                                cuda::stream_ref stream) const noexcept
  {
    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return 0; }

    auto counter =
      detail::counter_storage<size_type, thread_scope, allocator_type>{this->allocator(), stream};
    counter.reset(stream);

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::count<IsOuter, cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, counter.data(), container_ref);

    return counter.load_to_host(stream);
  }

  /**
   * @brief Counts the number of occurrences of each query key in the container
   *
   * For each key in the input range `[first, last)`, this function computes the number of matching
   * elements in the container and writes the result to the corresponding position in the output
   * range starting at `output_begin`.
   *
   * If `IsOuter` is `true` and a query key has no matches in the container, the result for that key
   * will be 1 instead of 0. Otherwise, the actual number of matches is returned.
   *
   * @note The input and output ranges must be device-accessible and of the same length.
   * @note The behavior is undefined if the input and output ranges overlap.
   *
   * @tparam IsOuter Flag indicating whether to use outer semantics (non-matches count as 1)
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the sequence of keys to count
   * @param last End of the sequence of keys to count
   * @param output_begin Beginning of the output range where per-key counts will be stored
   * @param container_ref Non-owning device container ref used to access the slot storage
   * @param stream CUDA stream used for this operation
   */
  template <bool IsOuter, typename InputIt, typename OutputIt, typename Ref>
  void count_each(InputIt first,
                  InputIt last,
                  OutputIt output_begin,
                  Ref container_ref,
                  cuda::stream_ref stream) const noexcept
  {
    auto const num_keys = cuco::detail::distance(first, last);
    if (num_keys == 0) { return; }

    auto const grid_size = cuco::detail::grid_size(num_keys, cg_size);

    detail::open_addressing_ns::count_each<IsOuter, cg_size, cuco::detail::default_block_size()>
      <<<grid_size, cuco::detail::default_block_size(), 0, stream.get()>>>(
        first, num_keys, output_begin, container_ref);
  }

  /**
   * @brief Retrieves all the slots corresponding to all keys in the range `[first, last)`.
   *
   * If key `k = *(first + i)` exists in the container, copies `k` to `output_probe` and associated
   * slot contents to `output_match`, respectively. The output order is unspecified.
   *
   * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
   * Use `count()/count_outer()` to determine the size of the output range.
   *
   * If `IsOuter == true` and a key `k` has no matches in the container, then `{key,
   * empty_slot_sentinel}` will be added to the output sequence.
   *
   * This function synchronizes the given CUDA stream.
   *
   * @tparam IsOuter Flag indicating if an inner or outer retrieve operation should be performed
   * @tparam InputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   * @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `key_type`
   * @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
   * convertible to the container's `value_type`
   * @tparam Ref Type of non-owning device container ref allowing access to storage
   *
   * @param first Beginning of the input sequence of keys
   * @param last End of the input sequence of keys
   * @param output_probe Beginning of the sequence of keys corresponding to matching elements in
   * `output_match`
   * @param output_match Beginning of the sequence of matching elements
   * @param container_ref Non-owning device reference to the container
   * @param stream CUDA stream this operation is executed in
   *
   * @return Iterator pair indicating the the end of the output sequences
   */
  template <bool IsOuter, class InputProbeIt, class OutputProbeIt, class OutputMatchIt, class Ref>
  std::pair<OutputProbeIt, OutputMatchIt> retrieve_impl(InputProbeIt first,
                                                        InputProbeIt last,
                                                        OutputProbeIt output_probe,
                                                        OutputMatchIt output_match,
                                                        Ref container_ref,
                                                        cuda::stream_ref stream) const
  {
    auto const n = detail::distance(first, last);
    if (n == 0) { return {output_probe, output_match}; }

    using counter_type = detail::counter_storage<size_type, thread_scope, allocator_type>;
    auto counter       = counter_type{this->allocator(), stream};
    counter.reset(stream.get());

    auto constexpr block_size = cuco::detail::default_block_size();

    auto constexpr grid_stride = 1;
    auto const grid_size       = cuco::detail::grid_size(n, cg_size, grid_stride, block_size);

    detail::open_addressing_ns::retrieve<IsOuter, block_size>
      <<<grid_size, block_size, 0, stream.get()>>>(
        first, n, output_probe, output_match, counter.data(), container_ref);

    auto const num_retrieved = counter.load_to_host(stream.get());

    return {output_probe + num_retrieved, output_match + num_retrieved};
  }

  /**
   * @brief Extracts the key from a given slot.
   *
   * @param value The input value
   *
   * @return The key
   */
  [[nodiscard]] constexpr key_type const& extract_key(value_type const& slot) const noexcept
  {
    if constexpr (has_payload) {
      return slot.first;
    } else {
      return slot;
    }
  }

 protected:
  // TODO: cleanup by using equal wrapper as a data member
  value_type empty_slot_sentinel_;      ///< Slot value that represents an empty slot
  key_type erased_key_sentinel_;        ///< Key value that represents an erased slot
  key_equal predicate_;                 ///< Key equality binary predicate
  probing_scheme_type probing_scheme_;  ///< Probing scheme
  storage_type storage_;                ///< Slot bucket storage
};

}  // namespace detail
}  // namespace cuco
