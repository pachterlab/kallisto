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

#include <cuco/detail/dynamic_map_kernels.cuh>
#include <cuco/hash_functions.cuh>
#include <cuco/static_map.cuh>
#include <cuco/types.cuh>

#include <cuda/std/atomic>
#include <cuda/std/functional>
#include <thrust/device_vector.h>

#include <cstddef>
#include <memory>
#include <numeric>
#include <type_traits>
#include <vector>

namespace cuco {

namespace experimental {
/**
 * @brief A GPU-accelerated, unordered, associative container of key-value
 * pairs with unique keys.
 *
 * This container automatically grows its capacity as necessary until device memory runs out.
 *
 * @tparam Key The type of the keys.
 * @tparam T The type of the mapped values.
 * @tparam Extent The type representing the extent of the container.
 * @tparam Scope The thread scope for the container's operations.
 * @tparam KeyEqual The equality comparison function for keys.
 * @tparam ProbingScheme The probing scheme for resolving hash collisions.
 * @tparam Allocator The allocator used for memory management.
 * @tparam Storage The storage policy for the container.
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
class dynamic_map {
  using map_type = static_map<Key, T, Extent, Scope, KeyEqual, ProbingScheme, Allocator, Storage>;

 public:
  static constexpr auto thread_scope = map_type::thread_scope;  ///< CUDA thread scope

  using key_type    = typename map_type::key_type;    ///< Key type
  using value_type  = typename map_type::value_type;  ///< Key-value pair type
  using size_type   = typename map_type::size_type;   ///< Size type
  using key_equal   = typename map_type::key_equal;   ///< Key equality comparator type
  using mapped_type = T;                              ///< Payload type

  dynamic_map(dynamic_map const&)            = delete;
  dynamic_map& operator=(dynamic_map const&) = delete;

  dynamic_map(dynamic_map&&) = default;  ///< Move constructor

  /**
   * @brief Replaces the contents of the container with another container.
   *
   * @return Reference of the current map object
   */
  dynamic_map& operator=(dynamic_map&&) = default;
  ~dynamic_map()                        = default;

  /**
   * @brief Constructs a dynamically-sized map with erase capability.
   *
   * The capacity of the map will automatically increase as the user adds key/value pairs using
   * `insert`.
   *
   * Capacity increases by a factor of growth_factor each time the size of the map exceeds a
   * threshold occupancy. The performance of `find` and `contains` gradually decreases each time the
   * map's capacity grows.
   *
   * @param initial_capacity The initial number of slots in the map
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param empty_value_sentinel The reserved mapped value for empty slots
   * @param pred Key equality binary predicate
   * @param probing_scheme Probing scheme
   * @param scope The scope in which operations will be performed
   * @param storage Kind of storage to use
   * @param alloc Allocator used for allocating device storage
   * @param stream CUDA stream used to initialize the map
   */
  constexpr dynamic_map(Extent initial_capacity,
                        empty_key<Key> empty_key_sentinel,
                        empty_value<T> empty_value_sentinel,
                        KeyEqual const& pred                = {},
                        ProbingScheme const& probing_scheme = {},
                        cuda_thread_scope<Scope> scope      = {},
                        Storage storage                     = {},
                        Allocator const& alloc              = {},
                        cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Grows the capacity of the map so there is enough space for `n` key/value pairs.
   *
   * If there is already enough space for `n` key/value pairs, the capacity remains the same.
   *
   * @param n The number of key value pairs for which there must be space
   * @param stream Stream used for executing the kernels
   */
  void reserve(size_type n, cuda::stream_ref stream);

  /**
   * @brief Inserts all key/value pairs in the range `[first, last)`.
   *
   * If multiple keys in `[first, last)` compare equal, it is unspecified which
   * element is inserted.
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `value_type`
   * @param first Beginning of the sequence of key/value pairs
   * @param last End of the sequence of key/value pairs
   * @param stream Stream used for executing the kernels
   */
  template <typename InputIt>
  void insert(InputIt first,
              InputIt last,
              cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Indicates whether the keys in the range `[first, last)` are contained in the map.
   *
   * Writes a `bool` to `(output + i)` indicating if the key `*(first + i)` exists in the map.
   *
   * @tparam InputIt Device accessible input iterator
   * @tparam OutputIt Device accessible output iterator whose `value_type` is
   * convertible to the map's `mapped_type`
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

 private:
  size_type size_{};      ///< Number of keys in the map
  size_type capacity_{};  ///< Maximum number of keys that can be inserted

  std::vector<std::unique_ptr<map_type>> submaps_;  ///< vector of pointers to each submap
  size_type min_insert_size_{};                     ///< min remaining capacity of submap for insert
  float max_load_factor_{};                         ///< Maximum load factor
  Allocator alloc_{};  ///< Allocator passed to submaps to allocate their device storage
};

}  // namespace experimental

/**
 * @brief A GPU-accelerated, unordered, associative container of key-value
 * pairs with unique keys
 *
 * Automatically grows capacity as necessary until device memory runs out.
 *
 * Allows constant time concurrent inserts or concurrent find operations (not
 * concurrent insert and find) from threads in device code.
 *
 * Current limitations:
 * - Requires keys and values that where `cuco::is_bitwise_comparable_v<T>` is true
 *    - Comparisons against the "sentinel" values will always be done with bitwise comparisons.
 * - Capacity does not shrink automatically
 * - Requires the user to specify sentinel values for both key and mapped value
 *   to indicate empty slots
 * - Does not support concurrent insert and find operations
 *
 * The `dynamic_map` supports host-side "bulk" operations which include `insert`, `find`
 * and `contains`. These are to be used when there are a large number of keys to insert
 * or lookup in the map. For example, given a range of keys specified by device-accessible
 * iterators, the bulk `insert` function will insert all keys into the map.
 *
 * Example:
 * \code{.cpp}
 * int empty_key_sentinel = -1;
 * int empty_value_sentinel = -1;
 *
 * // Constructs a map with 100,000 initial slots using -1 and -1 as the empty key/value
 * // sentinels. Performs one bulk insert of 50,000 keys and a second bulk insert of
 * // 100,000 keys. The map automatically increases capacity to accomodate the excess keys
 * // within the second insert.
 *
 * dynamic_map<int, int> m{100'000,
 *                         empty_key<int>{empty_key_sentinel},
 *                         empty_value<int>{empty_value_sentinel}};
 *
 * // Create a sequence of pairs {{0,0}, {1,1}, ... {i,i}}
 * thrust::device_vector<cuda::std::pair<int,int>> pairs_0(50'000);
 * thrust::transform(thrust::make_counting_iterator(0),
 *                   thrust::make_counting_iterator(pairs_0.size()),
 *                   pairs_0.begin(),
 *                   []__device__(auto i){ return cuco::pair{i,i}; };
 *
 * thrust::device_vector<cuda::std::pair<int,int>> pairs_1(100'000);
 * thrust::transform(thrust::make_counting_iterator(50'000),
 *                   thrust::make_counting_iterator(pairs_1.size()),
 *                   pairs_1.begin(),
 *                   []__device__(auto i){ return cuco::pair{i,i}; };
 *
 * // Inserts all pairs into the map
 * m.insert(pairs_0.begin(), pairs_0.end());
 * m.insert(pairs_1.begin(), pairs_1.end());
 * \endcode
 *
 * @tparam Key Arithmetic type used for key
 * @tparam Value Type of the mapped values
 * @tparam Scope The scope in which insert/find/contains will be performed by
 * individual threads.
 * @tparam Allocator Type of allocator used to allocate submap device storage
 */
template <typename Key,
          typename Value,
          cuda::thread_scope Scope = cuda::thread_scope_device,
          typename Allocator       = cuco::cuda_allocator<char>>
class dynamic_map {
  static_assert(std::is_arithmetic<Key>::value, "Unsupported, non-arithmetic key type.");

 public:
  using value_type      = cuco::pair<Key, Value>;            ///< Type of key/value pairs
  using key_type        = Key;                               ///< Key type
  using mapped_type     = Value;                             ///< Type of mapped values
  using atomic_ctr_type = cuda::atomic<std::size_t, Scope>;  ///< Atomic counter type
  using view_type =
    typename cuco::legacy::static_map<Key, Value, Scope>::device_view;  ///< Type for submap device
                                                                        ///< view
  using mutable_view_type =
    typename cuco::legacy::static_map<Key, Value, Scope>::device_mutable_view;  ///< Type for submap
                                                                                ///< mutable device
                                                                                ///< view

  dynamic_map(dynamic_map const&) = delete;
  dynamic_map(dynamic_map&&)      = delete;

  dynamic_map& operator=(dynamic_map const&) = delete;
  dynamic_map& operator=(dynamic_map&&)      = delete;

  /**
   * @brief Constructs a dynamically-sized map with the specified initial capacity, growth factor
   * and sentinel values.
   *
   * The capacity of the map will automatically increase as the user adds key/value pairs using
   * `insert`.
   *
   * Capacity increases by a factor of growth_factor each time the size of the map exceeds a
   * threshold occupancy. The performance of `find` and `contains` decreases somewhat each time the
   * map's capacity grows.
   *
   * The `empty_key_sentinel` and `empty_value_sentinel` values are reserved and
   * undefined behavior results from attempting to insert any key/value pair
   * that contains either.
   *
   * @param initial_capacity The initial number of slots in the map
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param empty_value_sentinel The reserved mapped value for empty slots
   * @param alloc Allocator used to allocate submap device storage
   * @param stream Stream used for executing the kernels
   */
  dynamic_map(std::size_t initial_capacity,
              empty_key<Key> empty_key_sentinel,
              empty_value<Value> empty_value_sentinel,
              Allocator const& alloc = Allocator{},
              cudaStream_t stream    = nullptr);

  /**
   * @brief Constructs a dynamically-sized map with erase capability.
   *
   * The capacity of the map will automatically increase as the user adds key/value pairs using
   * `insert`.
   *
   * Capacity increases by a factor of growth_factor each time the size of the map exceeds a
   * threshold occupancy. The performance of `find` and `contains` decreases somewhat each time the
   * map's capacity grows.
   *
   * The `empty_key_sentinel` and `empty_value_sentinel` values are reserved and
   * undefined behavior results from attempting to insert any key/value pair
   * that contains either.
   *
   * @param initial_capacity The initial number of slots in the map
   * @param empty_key_sentinel The reserved key value for empty slots
   * @param empty_value_sentinel The reserved mapped value for empty slots
   * @param erased_key_sentinel The reserved key value for erased slots
   * @param alloc Allocator used to allocate submap device storage
   * @param stream Stream used for executing the kernels
   *
   * @throw std::runtime error if the empty key sentinel and erased key sentinel
   * are the same value
   */
  dynamic_map(std::size_t initial_capacity,
              empty_key<Key> empty_key_sentinel,
              empty_value<Value> empty_value_sentinel,
              erased_key<Key> erased_key_sentinel,
              Allocator const& alloc = Allocator{},
              cudaStream_t stream    = nullptr);

  /**
   * @brief Destroys the map and frees its contents
   *
   */
  ~dynamic_map() {}

  /**
   * @brief Grows the capacity of the map so there is enough space for `n` key/value pairs.
   *
   * If there is already enough space for `n` key/value pairs, the capacity remains the same.
   *
   * @param n The number of key value pairs for which there must be space
   * @param stream Stream used for executing the kernels
   */
  void reserve(std::size_t n, cudaStream_t stream = nullptr);

  /**
   * @brief Inserts all key/value pairs in the range `[first, last)`.
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
              cudaStream_t stream = nullptr);

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
   * Keep in mind that `erase` does not cause the map to shrink its memory allocation.
   *
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `value_type`
   * @tparam Hash Unary callable type
   * @tparam KeyEqual Binary callable type
   *
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
             cudaStream_t stream = nullptr);

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
   *
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
            cudaStream_t stream = nullptr);

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
   * @tparam InputIt Device accessible input iterator whose `value_type` is
   * convertible to the map's `key_type`
   * @tparam OutputIt Device accessible output iterator whose `value_type` is
   * convertible to the map's `mapped_type`
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
                cudaStream_t stream = nullptr);

  /**
   * @brief Gets the current number of elements in the map
   *
   * @return The current number of elements in the map
   */
  std::size_t get_size() const noexcept { return size_; }

  /**
   * @brief Gets the maximum number of elements the hash map can hold.
   *
   * @return The maximum number of elements the hash map can hold
   */
  std::size_t get_capacity() const noexcept { return capacity_; }

  /**
   * @brief Gets the load factor of the hash map.
   *
   * @return The load factor of the hash map
   */
  float get_load_factor() const noexcept { return static_cast<float>(size_) / capacity_; }

 private:
  key_type empty_key_sentinel_{};       ///< Key value that represents an empty slot
  mapped_type empty_value_sentinel_{};  ///< Initial value of empty slot
  key_type erased_key_sentinel_{};      ///< Key value that represents an erased slot

  // TODO: initialize this
  std::size_t size_{};       ///< Number of keys in the map
  std::size_t capacity_{};   ///< Maximum number of keys that can be inserted
  float max_load_factor_{};  ///< Max load factor before capacity growth

  std::vector<std::unique_ptr<cuco::legacy::static_map<key_type, mapped_type, Scope>>>
    submaps_;                                      ///< vector of pointers to each submap
  thrust::device_vector<view_type> submap_views_;  ///< vector of device views for each submap
  thrust::device_vector<mutable_view_type>
    submap_mutable_views_;         ///< vector of mutable device views for each submap
  std::size_t min_insert_size_{};  ///< min remaining capacity of submap for insert
  thrust::device_vector<atomic_ctr_type*>
    submap_num_successes_;  ///< Number of successfully erased keys for each submap
  Allocator alloc_{};       ///< Allocator passed to submaps to allocate their device storage
};
}  // namespace cuco

#include <cuco/detail/dynamic_map.inl>
#include <cuco/detail/dynamic_map/dynamic_map.inl>
