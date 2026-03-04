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

#include <cuco/detail/hash_functions/utils.cuh>
#include <cuco/extent.cuh>

#include <cuda/std/array>
#include <cuda/std/cstddef>
#include <cuda/std/type_traits>

#include <cstdint>

namespace cuco::detail {

/**
 * @brief The 32-bit integer finalizer function of `MurmurHash3`.
 *
 * This function implements the final mixing step of the `MurmurHash3` algorithm for 32-bit values.
 * It is designed to improve the avalanche behavior of the hash, ensuring that changes in input bits
 * have a more uniform effect on all output bits.
 *
 * @throw Key type must be 4 bytes in size
 *
 * @tparam Key The type of the value to finalize
 *
 * @param key The input value to finalize
 * @param seed Optional seed value
 * @return The finalized 32-bit hash value
 */
template <typename Key>
__host__ __device__ constexpr std::uint32_t fmix32(Key key, std::uint32_t seed = 0) noexcept
{
  static_assert(sizeof(Key) == 4, "Key type must be 4 bytes in size.");

  std::uint32_t h = static_cast<std::uint32_t>(key) ^ seed;
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

/**
 * @brief The 64-bit integer finalizer function of `MurmurHash3`.
 *
 * This function implements the final mixing step of the `MurmurHash3` algorithm for 64-bit values.
 * It is designed to improve the avalanche behavior of the hash, ensuring that changes in input bits
 * have a more uniform effect on all output bits.
 *
 * @throw Key type must be 8 bytes in size
 *
 * @tparam Key The type of the value to finalize
 *
 * @param key The input value to finalize
 * @param seed Optional seed value
 * @return The finalized 64-bit hash value
 */
template <typename Key>
__host__ __device__ constexpr std::uint64_t fmix64(Key key, std::uint64_t seed = 0) noexcept
{
  static_assert(sizeof(Key) == 8, "Key type must be 8 bytes in size.");

  std::uint64_t h = static_cast<std::uint64_t>(key) ^ seed;
  h ^= h >> 33;
  h *= 0xff51afd7ed558ccdULL;
  h ^= h >> 33;
  h *= 0xc4ceb9fe1a85ec53ULL;
  h ^= h >> 33;
  return h;
}

/**
 * @brief The 32bit integer finalizer hash function of `MurmurHash3`.
 *
 * @tparam Key The type of the values to hash
 */
template <typename Key>
struct MurmurHash3_fmix32 {
  using argument_type = Key;            ///< The type of the values taken as argument
  using result_type   = std::uint32_t;  ///< The type of the hash values produced

  /**
   * @brief Constructs a MurmurHash3_fmix32 hash function with the given `seed`.
   *
   * @param seed A custom number to randomize the resulting hash value
   */
  __host__ __device__ constexpr MurmurHash3_fmix32(std::uint32_t seed = 0) : seed_{seed} {}

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @param key The input argument to hash
   * @return A resulting hash value for `key`
   */
  constexpr result_type __host__ __device__ operator()(Key const& key) const noexcept
  {
    return fmix32(key, seed_);
  }

 private:
  std::uint32_t seed_;
};

/**
 * @brief The 64bit integer finalizer hash function of `MurmurHash3`.
 *
 * @tparam Key The type of the values to hash
 */
template <typename Key>
struct MurmurHash3_fmix64 {
  using argument_type = Key;            ///< The type of the values taken as argument
  using result_type   = std::uint64_t;  ///< The type of the hash values produced

  /**
   * @brief Constructs a MurmurHash3_fmix64 hash function with the given `seed`.
   *
   * @param seed A custom number to randomize the resulting hash value
   */
  __host__ __device__ constexpr MurmurHash3_fmix64(std::uint64_t seed = 0) : seed_{seed} {}

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @param key The input argument to hash
   * @return A resulting hash value for `key`
   */
  constexpr result_type __host__ __device__ operator()(Key const& key) const noexcept
  {
    return fmix64(key, seed_);
  }

 private:
  std::uint64_t seed_;
};

/**
 * @brief A `MurmurHash3_32` hash function to hash the given argument on host and device.
 *
 * MurmurHash3_32 implementation from
 * https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
 * -----------------------------------------------------------------------------
 * MurmurHash3 was written by Austin Appleby, and is placed in the public domain. The author
 * hereby disclaims copyright to this source code.
 *
 * Note - The x86 and x64 versions do _not_ produce the same results, as the algorithms are
 * optimized for their respective platforms. You can still compile and run any of them on any
 * platform, but your performance with the non-native version will be less than optimal.
 *
 * @tparam Key The type of the values to hash
 */
template <typename Key>
struct MurmurHash3_32 {
  using argument_type = Key;            ///< The type of the values taken as argument
  using result_type   = std::uint32_t;  ///< The type of the hash values produced

  /**
   * @brief Constructs a MurmurHash3_32 hash function with the given `seed`.
   *
   * @param seed A custom number to randomize the resulting hash value
   */
  __host__ __device__ constexpr MurmurHash3_32(std::uint32_t seed = 0) : seed_{seed} {}

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @param key The input argument to hash
   * @return The resulting hash value for `key`
   */
  constexpr result_type __host__ __device__ operator()(Key const& key) const noexcept
  {
    return compute_hash(reinterpret_cast<cuda::std::byte const*>(&key),
                        cuco::extent<std::size_t, sizeof(Key)>{});
  }

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @tparam Extent The extent type
   *
   * @param bytes The input argument to hash
   * @param size The extent of the data in bytes
   * @return The resulting hash value
   */
  template <typename Extent>
  constexpr result_type __host__ __device__ compute_hash(cuda::std::byte const* bytes,
                                                         Extent size) const noexcept
  {
    auto const nblocks = size / 4;

    std::uint32_t h1           = seed_;
    constexpr std::uint32_t c1 = 0xcc9e2d51;
    constexpr std::uint32_t c2 = 0x1b873593;
    //----------
    // body
    for (cuda::std::remove_const_t<decltype(nblocks)> i = 0; size >= 4 && i < nblocks; i++) {
      std::uint32_t k1 = load_chunk<std::uint32_t>(bytes, i);
      k1 *= c1;
      k1 = rotl32(k1, 15);
      k1 *= c2;
      h1 ^= k1;
      h1 = rotl32(h1, 13);
      h1 = h1 * 5 + 0xe6546b64;
    }
    //----------
    // tail
    std::uint32_t k1 = 0;
    switch (size & 3) {
      case 3:
        k1 ^= cuda::std::to_integer<std::uint32_t>(bytes[nblocks * 4 + 2]) << 16;
        [[fallthrough]];
      case 2:
        k1 ^= cuda::std::to_integer<std::uint32_t>(bytes[nblocks * 4 + 1]) << 8;
        [[fallthrough]];
      case 1:
        k1 ^= cuda::std::to_integer<std::uint32_t>(bytes[nblocks * 4 + 0]);
        k1 *= c1;
        k1 = rotl32(k1, 15);
        k1 *= c2;
        h1 ^= k1;
        [[fallthrough]];
      case 0: break;
    };
    //----------
    // finalization
    h1 ^= size;
    h1 = fmix32(h1);
    return h1;
  }

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @note This API is to ensure backward compatibility with existing use cases using `std::byte`.
   * Users are encouraged to use the appropriate `cuda::std::byte` overload whenever possible for
   * better support and performance on the device.
   *
   * @tparam Extent The extent type
   *
   * @param bytes The input argument to hash
   * @param size The extent of the data in bytes
   * @return The resulting hash value
   */
  template <typename Extent>
  constexpr result_type __host__ __device__ compute_hash(std::byte const* bytes,
                                                         Extent size) const noexcept
  {
    return this->compute_hash(reinterpret_cast<cuda::std::byte const*>(bytes), size);
  }

 private:
  std::uint32_t seed_;
};

/**
 * @brief A `MurmurHash3_x64_128` hash function to hash the given argument on host and device.
 *
 * MurmurHash3_x64_128 implementation from
 * https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
 * -----------------------------------------------------------------------------
 * MurmurHash3 was written by Austin Appleby, and is placed in the public domain. The author
 * hereby disclaims copyright to this source code.
 *
 * Note - The x86 and x64 versions do _not_ produce the same results, as the algorithms are
 * optimized for their respective platforms. You can still compile and run any of them on any
 * platform, but your performance with the non-native version will be less than optimal.
 *
 * @tparam Key The type of the values to hash
 */
template <typename Key>
struct MurmurHash3_x64_128 {
  using argument_type = Key;  ///< The type of the values taken as argument
  using result_type = cuda::std::array<std::uint64_t, 2>;  ///< The type of the hash values produced

  /**
   * @brief Constructs a MurmurHash3_x64_128 hash function with the given `seed`.
   *
   * @param seed A custom number to randomize the resulting hash value
   */
  __host__ __device__ constexpr MurmurHash3_x64_128(std::uint64_t seed = 0) : seed_{seed} {}

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @param key The input argument to hash
   * @return The resulting hash value for `key`
   */
  constexpr result_type __host__ __device__ operator()(Key const& key) const noexcept
  {
    return compute_hash(reinterpret_cast<cuda::std::byte const*>(&key),
                        cuco::extent<std::size_t, sizeof(Key)>{});
  }

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @tparam Extent The extent type
   *
   * @param bytes The input argument to hash
   * @param size The extent of the data in bytes
   * @return The resulting hash value
   */
  template <typename Extent>
  constexpr result_type __host__ __device__ compute_hash(cuda::std::byte const* bytes,
                                                         Extent size) const noexcept
  {
    constexpr std::uint32_t block_size = 16;
    auto const nblocks                 = size / block_size;

    std::uint64_t h1           = seed_;
    std::uint64_t h2           = seed_;
    constexpr std::uint64_t c1 = 0x87c37b91114253d5ull;
    constexpr std::uint64_t c2 = 0x4cf5ad432745937full;
    //----------
    // body
    for (cuda::std::remove_const_t<decltype(nblocks)> i = 0; size >= block_size && i < nblocks;
         i++) {
      std::uint64_t k1 = load_chunk<std::uint64_t>(bytes, 2 * i);
      std::uint64_t k2 = load_chunk<std::uint64_t>(bytes, 2 * i + 1);

      k1 *= c1;
      k1 = rotl64(k1, 31);
      k1 *= c2;

      h1 ^= k1;
      h1 = rotl64(h1, 27);
      h1 += h2;
      h1 = h1 * 5 + 0x52dce729;

      k2 *= c2;
      k2 = rotl64(k2, 33);
      k2 *= c1;

      h2 ^= k2;
      h2 = rotl64(h2, 31);
      h2 += h1;
      h2 = h2 * 5 + 0x38495ab5;
    }
    //----------
    // tail
    std::uint64_t k1 = 0;
    std::uint64_t k2 = 0;
    auto const tail  = reinterpret_cast<uint8_t const*>(bytes) + nblocks * block_size;
    switch (size & (block_size - 1)) {
      case 15: k2 ^= static_cast<std::uint64_t>(tail[14]) << 48; [[fallthrough]];
      case 14: k2 ^= static_cast<std::uint64_t>(tail[13]) << 40; [[fallthrough]];
      case 13: k2 ^= static_cast<std::uint64_t>(tail[12]) << 32; [[fallthrough]];
      case 12: k2 ^= static_cast<std::uint64_t>(tail[11]) << 24; [[fallthrough]];
      case 11: k2 ^= static_cast<std::uint64_t>(tail[10]) << 16; [[fallthrough]];
      case 10: k2 ^= static_cast<std::uint64_t>(tail[9]) << 8; [[fallthrough]];
      case 9:
        k2 ^= static_cast<std::uint64_t>(tail[8]) << 0;
        k2 *= c2;
        k2 = rotl64(k2, 33);
        k2 *= c1;
        h2 ^= k2;
        [[fallthrough]];

      case 8: k1 ^= static_cast<std::uint64_t>(tail[7]) << 56; [[fallthrough]];
      case 7: k1 ^= static_cast<std::uint64_t>(tail[6]) << 48; [[fallthrough]];
      case 6: k1 ^= static_cast<std::uint64_t>(tail[5]) << 40; [[fallthrough]];
      case 5: k1 ^= static_cast<std::uint64_t>(tail[4]) << 32; [[fallthrough]];
      case 4: k1 ^= static_cast<std::uint64_t>(tail[3]) << 24; [[fallthrough]];
      case 3: k1 ^= static_cast<std::uint64_t>(tail[2]) << 16; [[fallthrough]];
      case 2: k1 ^= static_cast<std::uint64_t>(tail[1]) << 8; [[fallthrough]];
      case 1:
        k1 ^= static_cast<std::uint64_t>(tail[0]) << 0;
        k1 *= c1;
        k1 = rotl64(k1, 31);
        k1 *= c2;
        h1 ^= k1;
        [[fallthrough]];
      case 0: break;
    };
    //----------
    // finalization
    h1 ^= size;
    h2 ^= size;

    h1 += h2;
    h2 += h1;

    h1 = fmix64(h1);
    h2 = fmix64(h2);

    h1 += h2;
    h2 += h1;

    return {h1, h2};
  }

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @note This API is to ensure backward compatibility with existing use cases using `std::byte`.
   * Users are encouraged to use the appropriate `cuda::std::byte` overload whenever possible for
   * better support and performance on the device.
   *
   * @tparam Extent The extent type
   *
   * @param bytes The input argument to hash
   * @param size The extent of the data in bytes
   * @return The resulting hash value
   */
  template <typename Extent>
  constexpr result_type __host__ __device__ compute_hash(std::byte const* bytes,
                                                         Extent size) const noexcept
  {
    return this->compute_hash(reinterpret_cast<cuda::std::byte const*>(bytes), size);
  }

 private:
  std::uint64_t seed_;
};

/**
 * @brief A `MurmurHash3_x86_128` hash function to hash the given argument on host and device.
 *
 * MurmurHash3_x86_128 implementation from
 * https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
 * -----------------------------------------------------------------------------
 * MurmurHash3 was written by Austin Appleby, and is placed in the public domain. The author
 * hereby disclaims copyright to this source code.
 *
 * Note - The x86 and x64 versions do _not_ produce the same results, as the algorithms are
 * optimized for their respective platforms. You can still compile and run any of them on any
 * platform, but your performance with the non-native version will be less than optimal.
 *
 * @tparam Key The type of the values to hash
 */
template <typename Key>
struct MurmurHash3_x86_128 {
  using argument_type = Key;  ///< The type of the values taken as argument
  using result_type = cuda::std::array<std::uint32_t, 4>;  ///< The type of the hash values produced

  /**
   * @brief Constructs a MurmurHash3_x86_128 hash function with the given `seed`.
   *
   * @param seed A custom number to randomize the resulting hash value
   */
  __host__ __device__ constexpr MurmurHash3_x86_128(std::uint32_t seed = 0) : seed_{seed} {}

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @param key The input argument to hash
   * @return The resulting hash value for `key`
   */
  constexpr result_type __host__ __device__ operator()(Key const& key) const noexcept
  {
    return compute_hash(reinterpret_cast<cuda::std::byte const*>(&key),
                        cuco::extent<std::size_t, sizeof(Key)>{});
  }

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @tparam Extent The extent type
   *
   * @param bytes The input argument to hash
   * @param size The extent of the data in bytes
   * @return The resulting hash value
   */
  template <typename Extent>
  constexpr result_type __host__ __device__ compute_hash(cuda::std::byte const* bytes,
                                                         Extent size) const noexcept
  {
    constexpr std::uint32_t block_size = 16;
    auto const nblocks                 = size / block_size;

    std::uint32_t h1 = seed_;
    std::uint32_t h2 = seed_;
    std::uint32_t h3 = seed_;
    std::uint32_t h4 = seed_;

    constexpr std::uint32_t c1 = 0x239b961b;
    constexpr std::uint32_t c2 = 0xab0e9789;
    constexpr std::uint32_t c3 = 0x38b34ae5;
    constexpr std::uint32_t c4 = 0xa1e38b93;
    //----------
    // body
    for (cuda::std::remove_const_t<decltype(nblocks)> i = 0; size >= block_size && i < nblocks;
         i++) {
      std::uint32_t k1 = load_chunk<std::uint32_t>(bytes, 4 * i);
      std::uint32_t k2 = load_chunk<std::uint32_t>(bytes, 4 * i + 1);
      std::uint32_t k3 = load_chunk<std::uint32_t>(bytes, 4 * i + 2);
      std::uint32_t k4 = load_chunk<std::uint32_t>(bytes, 4 * i + 3);

      k1 *= c1;
      k1 = rotl32(k1, 15);
      k1 *= c2;
      h1 ^= k1;

      h1 = rotl32(h1, 19);
      h1 += h2;
      h1 = h1 * 5 + 0x561ccd1b;

      k2 *= c2;
      k2 = rotl32(k2, 16);
      k2 *= c3;
      h2 ^= k2;

      h2 = rotl32(h2, 17);
      h2 += h3;
      h2 = h2 * 5 + 0x0bcaa747;

      k3 *= c3;
      k3 = rotl32(k3, 17);
      k3 *= c4;
      h3 ^= k3;

      h3 = rotl32(h3, 15);
      h3 += h4;
      h3 = h3 * 5 + 0x96cd1c35;

      k4 *= c4;
      k4 = rotl32(k4, 18);
      k4 *= c1;
      h4 ^= k4;

      h4 = rotl32(h4, 13);
      h4 += h1;
      h4 = h4 * 5 + 0x32ac3b17;
    }
    //----------
    // tail
    std::uint32_t k1 = 0;
    std::uint32_t k2 = 0;
    std::uint32_t k3 = 0;
    std::uint32_t k4 = 0;

    auto const tail = reinterpret_cast<uint8_t const*>(bytes) + nblocks * block_size;
    switch (size & (block_size - 1)) {
      case 15: k4 ^= static_cast<std::uint32_t>(tail[14]) << 16; [[fallthrough]];
      case 14: k4 ^= static_cast<std::uint32_t>(tail[13]) << 8; [[fallthrough]];
      case 13:
        k4 ^= static_cast<std::uint32_t>(tail[12]) << 0;
        k4 *= c4;
        k4 = rotl32(k4, 18);
        k4 *= c1;
        h4 ^= k4;
        [[fallthrough]];

      case 12: k3 ^= static_cast<std::uint32_t>(tail[11]) << 24; [[fallthrough]];
      case 11: k3 ^= static_cast<std::uint32_t>(tail[10]) << 16; [[fallthrough]];
      case 10: k3 ^= static_cast<std::uint32_t>(tail[9]) << 8; [[fallthrough]];
      case 9:
        k3 ^= static_cast<std::uint32_t>(tail[8]) << 0;
        k3 *= c3;
        k3 = rotl32(k3, 17);
        k3 *= c4;
        h3 ^= k3;
        [[fallthrough]];

      case 8: k2 ^= static_cast<std::uint32_t>(tail[7]) << 24; [[fallthrough]];
      case 7: k2 ^= static_cast<std::uint32_t>(tail[6]) << 16; [[fallthrough]];
      case 6: k2 ^= static_cast<std::uint32_t>(tail[5]) << 8; [[fallthrough]];
      case 5:
        k2 ^= static_cast<std::uint32_t>(tail[4]) << 0;
        k2 *= c2;
        k2 = rotl32(k2, 16);
        k2 *= c3;
        h2 ^= k2;
        [[fallthrough]];

      case 4: k1 ^= static_cast<std::uint32_t>(tail[3]) << 24; [[fallthrough]];
      case 3: k1 ^= static_cast<std::uint32_t>(tail[2]) << 16; [[fallthrough]];
      case 2: k1 ^= static_cast<std::uint32_t>(tail[1]) << 8; [[fallthrough]];
      case 1:
        k1 ^= static_cast<std::uint32_t>(tail[0]) << 0;
        k1 *= c1;
        k1 = rotl32(k1, 15);
        k1 *= c2;
        h1 ^= k1;
        [[fallthrough]];
      case 0: break;
    };
    //----------
    // finalization
    h1 ^= size;
    h2 ^= size;
    h3 ^= size;
    h4 ^= size;

    h1 += h2;
    h1 += h3;
    h1 += h4;
    h2 += h1;
    h3 += h1;
    h4 += h1;

    h1 = fmix32(h1);
    h2 = fmix32(h2);
    h3 = fmix32(h3);
    h4 = fmix32(h4);

    h1 += h2;
    h1 += h3;
    h1 += h4;
    h2 += h1;
    h3 += h1;
    h4 += h1;

    return {h1, h2, h3, h4};
  }

  /**
   * @brief Returns a hash value for its argument, as a value of type `result_type`.
   *
   * @note This API is to ensure backward compatibility with existing use cases using `std::byte`.
   * Users are encouraged to use the appropriate `cuda::std::byte` overload whenever possible for
   * better support and performance on the device.
   *
   * @tparam Extent The extent type
   *
   * @param bytes The input argument to hash
   * @param size The extent of the data in bytes
   * @return The resulting hash value
   */
  template <typename Extent>
  constexpr result_type __host__ __device__ compute_hash(std::byte const* bytes,
                                                         Extent size) const noexcept
  {
    return this->compute_hash(reinterpret_cast<cuda::std::byte const*>(bytes), size);
  }

 private:
  std::uint32_t seed_;
};

}  //  namespace cuco::detail
