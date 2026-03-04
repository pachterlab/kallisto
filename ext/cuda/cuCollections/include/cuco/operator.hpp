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

namespace cuco {
// TODO point to container<->operator support matrix
/**
 * @brief Namespace containing device API operator tags
 *
 * This namespace defines tag types that are used to specify different operations
 * that can be performed on container reference objects in device code. These
 * tags serve as compile-time indicators for the type of operation being requested
 * in device code, and are not intended for use with host-side bulk APIs.
 *
 * Example Usage:
 * When applied to container reference object, e.g., `container_ref.rebind_operators(cuco::insert)`,
 * it enables using `container_ref.insert(...)`.
 * For a full example on how to use operators see `examples/static_set/device_ref_example.cu`.
 *
 * Available operators:
 * - `insert`: Inserts an element into the container
 * - `insert_and_find`: Inserts an element and returns an iterator to the stored location
 * - `insert_or_assign`: Inserts a new element or updates an existing element
 * - `insert_or_apply`: Inserts a new element or applies a user-defined function to an existing
 * element
 * - `erase`: Removes an element from the container
 * - `contains`: Checks element existence
 * - `count`: For a given key, returns the number of matching elements in the container
 * - `find`: Locates an element in the container
 * - `retrieve`: Retrieves all matching elements in the container for a given key
 * - `for_each`: Applies a user-defined function to each element in the container (or for a given
 * key)
 *
 */
inline namespace op {

/**
 * @brief Tag type for `insert` operator
 *
 * API Signature:
 * ```cpp
 * template <typename Value>
 * __device__ bool insert(Value const& value) noexcept
 *
 * template <typename Value, typename ParentCG>
 * __device__ bool insert(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
 *                        Value const& value) noexcept
 * ```
 *
 * Where:
 * @see @tparam Value Input type which is convertible to the container's `value_type`
 * @see @tparam ParentCG Type of parent Cooperative Group
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param value The element to insert
 *
 * @see @return `True` iff the given element is successfully inserted
 *
 */
struct insert_tag {
} inline constexpr insert;  ///< `cuco::insert` operator

/**
 * @brief Tag type for `insert_and_find` operator
 *
 * This API returns a pair consisting of an iterator to the inserted element (or to the
 * element that prevented the insertion) and a `bool` denoting whether the insertion took place or
 * not.
 *
 * API Signature:
 * ```cpp
 * template <typename Value>
 * __device__ cuda::std::pair<iterator, bool> insert_and_find(Value const& value) noexcept
 *
 * template <typename Value, typename ParentCG>
 * __device__ cuda::std::pair<iterator, bool> insert_and_find(
 *   cooperative_groups::thread_block_tile<cg_size, ParentCG> group, Value const& value) noexcept
 * ```
 *
 * Where:
 * @see @tparam Value Input type which is convertible to the container's `value_type`
 * @see @tparam ParentCG Type of parent Cooperative Group
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param value The element to insert
 *
 * @see @return a pair consisting of an iterator to the element and a bool indicating whether the
 * insertion is successful or not.
 *
 */
struct insert_and_find_tag {
} inline constexpr insert_and_find;  ///< `cuco::insert_and_find` operator

/**
 * @brief Tag type for `insert_or_assign` operator
 *
 * Inserts an element if it's not present in the container. Otherwise, assigns the input payload to
 * the existing element.
 *
 * API Signature:
 * ```cpp
 * template <typename Value>
 * __device__ void insert_or_assign(Value const& value) noexcept
 *
 * template <typename Value, typename ParentCG>
 * __device__ void insert_or_assign(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
 *                                  Value const& value) noexcept
 * ```
 *
 * Where:
 * @see @tparam Value Input type which is convertible to the container's `value_type`
 * @see @tparam ParentCG Type of parent Cooperative Group
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param value The element to insert
 *
 */
struct insert_or_assign_tag {
} inline constexpr insert_or_assign;  ///< `cuco::insert_or_assign` operator

/**
 * @brief Tag type for `insert_or apply` operator
 *
 * Inserts a new element or applies a user-defined function to an existing element.
 *
 * API Signature:
 * ```cpp
 * template <typename Value, typename Op>
 * __device__ bool insert_or_apply(Value const& value, Op op)
 *
 * template <typename Value,
 *           typename Init,
 *           typename Op>
 * __device__ bool insert_or_apply(Value const& value, Init init, Op op)
 *
 * template <typename Value, typename Op, typename ParentCG>
 * __device__ bool insert_or_apply(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
 *                                 Value const& value,
 *                                 Op op)
 *
 * template <typename Value, typename Init, typename Op, typename ParentCG>
 * __device__ bool insert_or_apply(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
 *                                 Value const& value,
 *                                 Init init,
 *                                 Op op)
 * ```
 *
 * Where:
 * @see @tparam Value Input type which is convertible to the container's `value_type`
 * @see @tparam Init Type of init value convertible to payload type
 * @see @tparam Op Callable type which is used as `apply` operation and can be
 *   called with arguments as `Op(cuda::atomic_ref<T, Scope>, T)`. `Op` strictly must
 *   have this signature to atomically apply the operation.
 * @see @tparam ParentCG Type of parent Cooperative Group
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param value The element to insert
 * @see @param init The init value of the op
 * @see @param op The callable object to perform binary operation between existing value at the slot
 *
 * @see @return Returns `true` if the given `value` is inserted successfully.
 *
 */
struct insert_or_apply_tag {
} inline constexpr insert_or_apply;  ///< `cuco::insert_or_apply` operator

/**
 * @brief Tag type for `erase` operator
 *
 * API Signature:
 * ```cpp
 * template <typename ProbeKey>
 * __device__ bool erase(ProbeKey const& key) noexcept
 *
 * template <typename ProbeKey, typename ParentCG>
 * __device__ bool erase(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
 *                       ProbeKey const& key) noexcept
 * ```
 *
 * Where:
 * @see @tparam ProbeKey Input key type which is convertible to the container's 'key_type'
 * @see @tparam ParentCG Type of parent Cooperative Group
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param key The key to search for
 *
 * @see @return 'True' if the given element is successfully erased
 *
 */
struct erase_tag {
} inline constexpr erase;  ///< `cuco::erase` operator

/**
 * @brief Tag type for `contains` operator
 *
 * API Signature:
 * ```cpp
 * template <typename ProbeKey>
 * __device__ bool contains(ProbeKey const& key) const noexcept
 *
 * template <typename ProbeKey, typename ParentCG>
 * __device__ bool contains(
 *   cooperative_groups::thread_block_tile<cg_size, ParentCG> group, ProbeKey const& key) const
 * noexcept
 * ```
 *
 * Where:
 * @see @tparam ProbeKey Input key type which is convertible to the containser's 'key_type'
 * @see @tparam ParentCG Type of parent Cooperative Group
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param key The key to search for
 *
 * @see @return A boolean indicating whether the probe key is present
 *
 */
struct contains_tag {
} inline constexpr contains;  ///< `cuco::contains` operator

/**
 * @brief Tag type for `count` operator
 *
 * API Signature:
 * ```cpp
 * template <typename ProbeKey>
 * __device__ size_type count(ProbeKey const& key) const noexcept
 *
 * template <typename ProbeKey, typename ParentCG>
 * __device__ size_type count(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
 *                            ProbeKey const& key) const noexcept
 * ```
 *
 * Where:
 * @see @tparam ProbeKey Input key type which is convertible to the containser's 'key_type'
 * @see @tparam ParentCG Type of parent Cooperative Group
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param key The key to search for
 *
 * @see @return Number of occurrences found by the current probing thread
 *
 */
struct count_tag {
} inline constexpr count;  ///< `cuco::contains` operator

/**
 * @brief Tag type for `find` operator
 *
 * API Signature:
 * ```cpp
 * template <typename ProbeKey>
 * __device__ const_iterator find(ProbeKey const& key) const noexcept
 *
 * template <typename ProbeKey, typename ParentCG>
 * __device__ const_iterator find(
 *   cooperative_groups::thread_block_tile<cg_size, ParentCG> group, ProbeKey const& key) const
 * noexcept
 * ```
 *
 * Where:
 * @see @tparam ProbeKey Input key type which is convertible to the containser's 'key_type'
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param key The key to search for
 *
 * @see @return An iterator to the position at which the equivalent element is stored
 *
 */
struct find_tag {
} inline constexpr find;  ///< `cuco::find` operator

/**
 * @brief Tag type for `retrieve` operator
 *
 * Retrieves all the matching elements corresponding to all keys in the range `[input_probe_begin,
 * input_probe_end)`.
 *
 * If key `k = *(first + i)` exists in the container, copies `k` to `output_probe` and associated
 * slot content to `output_match`, respectively. The output order is unspecified.
 *
 * Behavior is undefined if the size of the output range exceeds the number of retrieved slots.
 * Use `count()` to determine the size of the output range.
 *
 * API Signatures:
 * ```cpp
 * // Basic retrieve
 * template <int BlockSize,
 *           class InputProbeIt,
 *           class OutputProbeIt,
 *           class OutputMatchIt,
 *           class AtomicCounter>
 * __device__ void retrieve(cooperative_groups::thread_block  const& block,
 *                          InputProbeIt input_probe_begin,
 *                          InputProbeIt input_probe_end,
 *                          OutputProbeIt output_probe,
 *                          OutputMatchIt output_match,
 *                          AtomicCounter* atomic_counter) const
 *
 * // Conditional retrieve with predicate
 * template <int BlockSize,
 *           class InputProbeIt,
 *           class StencilIt,
 *           class Predicate,
 *           class OutputProbeIt,
 *           class OutputMatchIt,
 *           class AtomicCounter>
 * __device__ void retrieve_if(cooperative_groups::thread_block  const& block,
 *                             InputProbeIt input_probe_begin,
 *                             InputProbeIt input_probe_end,
 *                             StencilIt stencil,
 *                             Predicate pred,
 *                             OutputProbeIt output_probe,
 *                             OutputMatchIt output_match,
 *                             AtomicCounter* atomic_counter) const
 * ```
 *
 * Where:

 * @see @tparam BlockSize Size of the thread block this operation is executed in
 * @see @tparam InputProbeIt Device accessible input iterator whose `value_type` is
 * convertible to the container's `key_type`
 * @see @tparam StencilIt Device accessible random access iterator whose value_type is
 * convertible to Predicate's argument type (retrieve_if only)
 * @see @tparam Predicate Unary predicate callable whose return type must be convertible to `bool`
 * @see @tparam OutputProbeIt Device accessible input iterator whose `value_type` is
 * convertible to the container's `key_type`
 * @see @tparam OutputMatchIt Device accessible input iterator whose `value_type` is
 * convertible to the container's `value_type`
 * @see @tparam AtomicCounter Atomic counter type that follows the same semantics as
 * `cuda::atomic(_ref)`
 *
 * @see @param block Thread block this operation is executed in
 * @see @param input_probe_begin Beginning of the input sequence of keys
 * @see @param input_probe_end End of the input sequence of keys
 * @see @param stencil Beginning of the stencil sequence (retrieve_if only)
 * @see @param pred Predicate to test on every element in the range `[stencil, stencil + n)`
 (retrieve_if only)
 * @see @param output_probe Beginning of the sequence of keys corresponding to matching elements in
 * `output_match`
 * @see @param output_match Beginning of the sequence of matching elements
 * @see @param atomic_counter Counter that is used to determine the next free position in the output
 * sequences
 *
 */
struct retrieve_tag {
} inline constexpr retrieve;  ///< `cuco::retrieve` operator

/**
 * @brief Tag type for `for_each` operator
 *
 * Invokes a user-defined callback function on every element in the container with key equivalent to
 * the probe key and can additionally perform work that requires synchronizing the Cooperative Group
 * performing this operation.
 *
 * @note Passes an un-incrementable input iterator to the element whose key is equivalent to
 * `key` to the callback.
 *
 * @note Synchronizing `group` within `callback_op` is undefined behavior.
 *
 * @note The `sync_op` function can be used to perform work that requires synchronizing threads in
 * `group` inbetween probing steps, where the number of probing steps performed between
 * synchronization points is capped by `bucket_size * cg_size`. The functor will be called right
 * after the current probing bucket has been traversed.
 *
 * API Signature:
 * ```cpp
 * template <class ProbeKey, class CallbackOp>
 * __device__ void for_each(ProbeKey const& key, CallbackOp&& callback_op) const noexcept
 *
 * template <class ProbeKey, class CallbackOp, typename ParentCG>
 * __device__ void for_each(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
 *                          ProbeKey const& key,
 *                          CallbackOp&& callback_op) const noexcept
 *
 * template <class ProbeKey, class CallbackOp, class SyncOp, typename ParentCG>
 * __device__ void for_each(cooperative_groups::thread_block_tile<cg_size, ParentCG> group,
 *                          ProbeKey const& key,
 *                          CallbackOp&& callback_op,
 *                          SyncOp&& sync_op) const noexcept
 * ```
 *
 * Where:
 * @see @tparam ProbeKey Probe key type
 * @see @tparam CallbackOp Type of unary callback function object
 * @see @tparam SyncOp Functor or device lambda which accepts the current `group` object
 * @see @tparam ParentCG Type of parent Cooperative Group
 *
 * @see @param group The Cooperative Group used to perform this operation
 * @see @param key The key to search for
 * @see @param callback_op Function to call on every element found
 * @see @param sync_op Function that is allowed to synchronize `group` inbetween probing buckets
 *
 */
struct for_each_tag {
} inline constexpr for_each;  ///< `cuco::for_each` operator

}  // namespace op
}  // namespace cuco

#include <cuco/detail/operator.inl>
