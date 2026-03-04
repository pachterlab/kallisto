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

#include <cuco/detail/storage/storage_base.cuh>
#include <cuco/extent.cuh>
#include <cuco/utility/allocator.hpp>

#include <cuda/std/array>
#include <cuda/std/functional>
#include <cuda/stream_ref>

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>

namespace cuco {
/**
 * @brief Non-owning array of slots storage reference type.
 *
 * @tparam T Storage element type
 * @tparam BucketSize Number of slots in each bucket
 * @tparam Extent Type of extent denoting storage capacity
 */
template <typename T, int32_t BucketSize, typename Extent = cuco::extent<std::size_t>>
class bucket_storage_ref {
 public:
  static constexpr int32_t bucket_size = BucketSize;  ///< Number of elements processed per bucket
  static constexpr std::size_t alignment =
    cuda::std::min(sizeof(T) * bucket_size, std::size_t{16});  ///< Required alignment

  using extent_type = Extent;                            ///< Storage extent type
  using size_type   = typename extent_type::value_type;  ///< Storage size type
  using value_type  = T;                                 ///< Slot type
  using bucket_type = cuda::std::array<T, BucketSize>;   ///< Slot bucket type

  /**
   * @brief Constructor of slot storage ref.
   *
   * @param size Number of slots
   * @param slots Pointer to the slots array
   */
  __host__ __device__ explicit constexpr bucket_storage_ref(Extent size,
                                                            value_type* slots) noexcept;

  using iterator       = value_type*;        ///< Iterator type
  using const_iterator = value_type const*;  ///< Const forward iterator type

  /**
   * @brief Returns an iterator to one past the last slot.
   *
   * This is provided for convenience for those familiar with checking
   * an iterator returned from `find()` against the `end()` iterator.
   *
   * @return An iterator to one past the last slot
   */
  [[nodiscard]] __host__ __device__ constexpr iterator end() noexcept;

  /**
   * @brief Returns a const_iterator to one past the last slot.
   *
   * This is provided for convenience for those familiar with checking
   * an iterator returned from `find()` against the `end()` iterator.
   *
   * @return A const_iterator to one past the last slot
   */
  [[nodiscard]] __host__ __device__ constexpr iterator end() const noexcept;

  /**
   * @brief Gets slots array.
   *
   * @return Pointer to the first slot
   */
  [[nodiscard]] __host__ __device__ constexpr value_type* data() noexcept;

  /**
   * @brief Gets slots array.
   *
   * @return Pointer to the first slot
   */
  [[nodiscard]] __host__ __device__ constexpr value_type* data() const noexcept;

  /**
   * @brief Returns an array of slots (or a bucket) for a given index.
   *
   * @param index Index of the slot
   * @return An array of slots
   */
  [[nodiscard]] __device__ constexpr bucket_type operator[](size_type index) const noexcept;

  /**
   * @brief Gets the total number of slot buckets in the current storage.
   *
   * @return The total number of slot buckets
   */
  [[nodiscard]] __host__ __device__ constexpr size_type num_buckets() const noexcept;

  /**
   * @brief Gets the total number of slots in the current storage.
   *
   * @return The total number of slots
   */
  [[nodiscard]] __host__ __device__ constexpr size_type capacity() const noexcept;

  /**
   * @brief Gets the bucket extent of the current storage.
   *
   * @return The bucket extent.
   */
  [[nodiscard]] __host__ __device__ constexpr extent_type extent() const noexcept;

 private:
  extent_type extent_;  ///< Storage extent
  value_type* slots_;   ///< Pointer to the slots array
};

/**
 * @brief Array of slots open addressing storage class.
 *
 * @tparam T Slot type
 * @tparam BucketSize Number of slots in each bucket
 * @tparam Extent Type of extent denoting number of slots
 * @tparam Allocator Type of allocator used for device storage (de)allocation
 */
template <typename T,
          int32_t BucketSize,
          typename Extent    = cuco::extent<std::size_t>,
          typename Allocator = cuco::cuda_allocator<T>>
class bucket_storage {
 public:
  static constexpr int32_t bucket_size = BucketSize;  ///< Number of elements processed per bucket

  using extent_type = Extent;                            ///< Storage extent type
  using size_type   = typename extent_type::value_type;  ///< Storage size type
  using value_type  = T;                                 ///< Slot type
  using bucket_type = cuda::std::array<T, BucketSize>;   ///< Slot bucket type

  /// Type of the allocator to (de)allocate buckets
  using allocator_type =
    typename std::allocator_traits<Allocator>::template rebind_alloc<value_type>;
  using ref_type = bucket_storage_ref<value_type, bucket_size, extent_type>;  ///< Storage ref type

  /**
   * @brief Constructor of bucket slot storage.
   *
   * @note The input `size` should be exclusively determined by the return value of
   * `make_valid_extent` since it depends on the requested low-bound value, the probing scheme, and
   * the storage.
   *
   * @param size Number of slots to (de)allocate
   * @param allocator Allocator used for (de)allocating device storage
   * @param stream Stream to use for (de)allocating device storage
   */
  explicit constexpr bucket_storage(Extent size,
                                    Allocator const& allocator,
                                    cuda::stream_ref stream = cuda::stream_ref{
                                      cudaStream_t{nullptr}});

  bucket_storage(bucket_storage&&) = default;  ///< Move constructor
  /**
   * @brief Replaces the contents of the storage with another storage.
   *
   * @return Reference of the current storage object
   */
  bucket_storage& operator=(bucket_storage&&) = default;
  ~bucket_storage()                           = default;  ///< Destructor

  bucket_storage(bucket_storage const&)            = delete;
  bucket_storage& operator=(bucket_storage const&) = delete;

  /**
   * @brief Gets bucket slots array.
   *
   * @return Pointer to the first slot
   */
  [[nodiscard]] constexpr value_type* data() const noexcept;

  /**
   * @brief Gets the storage allocator.
   *
   * @return The storage allocator
   */
  [[nodiscard]] constexpr allocator_type allocator() const noexcept;

  /**
   * @brief Gets bucket storage reference.
   *
   * @return Reference of bucket storage
   */
  [[nodiscard]] constexpr ref_type ref() const noexcept;

  /**
   * @brief Initializes each slot in the bucket slot storage to contain `key`.
   *
   * @param key Key to which all keys in `slots` are initialized
   * @param stream Stream used for executing the kernel
   */
  void initialize(value_type key,
                  cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Asynchronously initializes each slot in the bucket storage to contain `key`.
   *
   * @param key Key to which all keys in `slots` are initialized
   * @param stream Stream used for executing the kernel
   */
  void initialize_async(value_type key,
                        cuda::stream_ref stream = cuda::stream_ref{cudaStream_t{nullptr}});

  /**
   * @brief Gets the total number of slot buckets in the current storage.
   *
   * @return The total number of slot buckets
   */
  [[nodiscard]] __host__ __device__ constexpr size_type num_buckets() const noexcept;

  /**
   * @brief Gets the total number of slots in the current storage.
   *
   * @return The total number of slots
   */
  [[nodiscard]] __host__ __device__ constexpr size_type capacity() const noexcept;

  /**
   * @brief Gets the bucket extent of the current storage.
   *
   * @return The bucket extent.
   */
  [[nodiscard]] __host__ __device__ constexpr extent_type extent() const noexcept;

 private:
  using slot_deleter_type =
    detail::custom_deleter<size_type, allocator_type>;  ///< Type of slot deleter

  extent_type extent_;        ///< Storage extent
  allocator_type allocator_;  ///< Allocator used to (de)allocate slots
  /// Pointer to the slot storage
  std::unique_ptr<value_type, slot_deleter_type> slots_;
};
}  // namespace cuco

#include <cuco/detail/storage/bucket_storage.inl>
