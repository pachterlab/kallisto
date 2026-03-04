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

#include <test_utils.hpp>

#include <cuco/detail/__config>
#include <cuco/hash_functions.cuh>

#include <cuda/std/cstddef>
#include <cuda/std/functional>
#include <cuda/std/limits>
#include <thrust/device_vector.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

template <int32_t Words>
struct large_key {
  constexpr __host__ __device__ large_key(int32_t value) noexcept
  {
    for (int32_t i = 0; i < Words; ++i) {
      data_[i] = value;
    }
  }

 private:
  int32_t data_[Words];
};

template <typename Hash, typename... HashConstructorArgs>
static __host__ __device__ bool check_hash_result(
  typename Hash::argument_type const& key,
  typename Hash::result_type expected,
  HashConstructorArgs&&... hash_constructor_args) noexcept
{
  Hash h(cuda::std::forward<HashConstructorArgs>(hash_constructor_args)...);
  return (h(key) == expected);
}

template <typename OutputIter>
__global__ void check_identity_hash_result_kernel(OutputIter result)
{
  int i = 0;

  result[i++] = check_hash_result<cuco::identity_hash<signed char>>(0, 0);
  result[i++] = check_hash_result<cuco::identity_hash<signed char>>(
    cuda::std::numeric_limits<signed char>::max(), cuda::std::numeric_limits<signed char>::max());

  result[i++] = check_hash_result<cuco::identity_hash<int32_t>>(0, 0);
  result[i++] = check_hash_result<cuco::identity_hash<int32_t>>(
    cuda::std::numeric_limits<int32_t>::max(), cuda::std::numeric_limits<int32_t>::max());

  result[i++] = check_hash_result<cuco::identity_hash<int64_t>>(0, 0);
  result[i++] = check_hash_result<cuco::identity_hash<int64_t>>(
    cuda::std::numeric_limits<int64_t>::max(), cuda::std::numeric_limits<int64_t>::max());
}

TEST_CASE("utility cuco::identity_hash test", "")
{
  SECTION("Check if host-generated hash values match the identity function.")
  {
    CHECK(check_hash_result<cuco::identity_hash<signed char>>(0, 0));
    CHECK(check_hash_result<cuco::identity_hash<signed char>>(
      cuda::std::numeric_limits<signed char>::max(),
      cuda::std::numeric_limits<signed char>::max()));

    CHECK(check_hash_result<cuco::identity_hash<int32_t>>(0, 0));
    CHECK(check_hash_result<cuco::identity_hash<int32_t>>(
      cuda::std::numeric_limits<int32_t>::max(), cuda::std::numeric_limits<int32_t>::max()));

    CHECK(check_hash_result<cuco::identity_hash<int64_t>>(0, 0));
    CHECK(check_hash_result<cuco::identity_hash<int64_t>>(
      cuda::std::numeric_limits<int64_t>::max(), cuda::std::numeric_limits<int64_t>::max()));
  }
  SECTION("Check if device-generated hash values match the identity function.")
  {
    thrust::device_vector<bool> result(7, true);

    check_identity_hash_result_kernel<<<1, 1>>>(result.begin());

    CHECK(cuco::test::all_of(result.begin(), result.end(), cuda::std::identity{}));
  }
}

template <typename OutputIter>
__global__ void check_hash_result_kernel_64(OutputIter result)
{
  int i = 0;

  result[i++] = check_hash_result<cuco::xxhash_64<char>>(0, 16804241149081757544, 0);
  result[i++] = check_hash_result<cuco::xxhash_64<char>>(42, 765293966243412708, 0);
  result[i++] = check_hash_result<cuco::xxhash_64<char>>(0, 9486749600008296231, 42);

  result[i++] = check_hash_result<cuco::xxhash_64<int32_t>>(0, 4246796580750024372, 0);
  result[i++] = check_hash_result<cuco::xxhash_64<int32_t>>(0, 3614696996920510707, 42);
  result[i++] = check_hash_result<cuco::xxhash_64<int32_t>>(42, 15516826743637085169, 0);
  result[i++] = check_hash_result<cuco::xxhash_64<int32_t>>(123456789, 9462334144942111946, 0);

  result[i++] = check_hash_result<cuco::xxhash_64<int64_t>>(0, 3803688792395291579, 0);
  result[i++] = check_hash_result<cuco::xxhash_64<int64_t>>(0, 13194218611613725804, 42);
  result[i++] = check_hash_result<cuco::xxhash_64<int64_t>>(42, 13066772586158965587, 0);
  result[i++] = check_hash_result<cuco::xxhash_64<int64_t>>(123456789, 14662639848940634189, 0);

#if defined(CUCO_HAS_INT128)
  result[i++] = check_hash_result<cuco::xxhash_64<__int128>>(123456789, 7986913354431084250, 0);
#endif

  result[i++] =
    check_hash_result<cuco::xxhash_64<large_key<32>>>(123456789, 2031761887105658523, 0);
}

TEST_CASE("utility cuco::xxhash_64 test", "")
{
  // Reference hash values were computed using https://github.com/Cyan4973/xxHash
  SECTION("Check if host-generated hash values match the reference implementation.")
  {
    CHECK(check_hash_result<cuco::xxhash_64<char>>(0, 16804241149081757544, 0));
    CHECK(check_hash_result<cuco::xxhash_64<char>>(42, 765293966243412708, 0));
    CHECK(check_hash_result<cuco::xxhash_64<char>>(0, 9486749600008296231, 42));

    CHECK(check_hash_result<cuco::xxhash_64<int32_t>>(0, 4246796580750024372, 0));
    CHECK(check_hash_result<cuco::xxhash_64<int32_t>>(0, 3614696996920510707, 42));
    CHECK(check_hash_result<cuco::xxhash_64<int32_t>>(42, 15516826743637085169, 0));
    CHECK(check_hash_result<cuco::xxhash_64<int32_t>>(123456789, 9462334144942111946, 0));

    CHECK(check_hash_result<cuco::xxhash_64<int64_t>>(0, 3803688792395291579, 0));
    CHECK(check_hash_result<cuco::xxhash_64<int64_t>>(0, 13194218611613725804, 42));
    CHECK(check_hash_result<cuco::xxhash_64<int64_t>>(42, 13066772586158965587, 0));
    CHECK(check_hash_result<cuco::xxhash_64<int64_t>>(123456789, 14662639848940634189, 0));

#if defined(CUCO_HAS_INT128)
    CHECK(check_hash_result<cuco::xxhash_64<__int128>>(123456789, 7986913354431084250, 0));
#endif

    // 32*4=128-byte key to test the pipelined outermost hashing loop
    CHECK(check_hash_result<cuco::xxhash_64<large_key<32>>>(123456789, 2031761887105658523, 0));
  }

  SECTION("Check if device-generated hash values match the reference implementation.")
  {
    thrust::device_vector<bool> result(10);

    check_hash_result_kernel_64<<<1, 1>>>(result.begin());

    CHECK(cuco::test::all_of(result.begin(), result.end(), cuda::std::identity{}));
  }
}

template <typename OutputIter>
__global__ void check_hash_result_kernel_32(OutputIter result)
{
  int i = 0;

  result[i++] = check_hash_result<cuco::xxhash_32<char>>(0, 3479547966, 0);
  result[i++] = check_hash_result<cuco::xxhash_32<char>>(42, 3774771295, 0);
  result[i++] = check_hash_result<cuco::xxhash_32<char>>(0, 2099223482, 42);

  result[i++] = check_hash_result<cuco::xxhash_32<int32_t>>(0, 148298089, 0);
  result[i++] = check_hash_result<cuco::xxhash_32<int32_t>>(0, 2132181312, 42);
  result[i++] = check_hash_result<cuco::xxhash_32<int32_t>>(42, 1161967057, 0);
  result[i++] = check_hash_result<cuco::xxhash_32<int32_t>>(123456789, 2987034094, 0);

  result[i++] = check_hash_result<cuco::xxhash_32<int64_t>>(0, 3736311059, 0);
  result[i++] = check_hash_result<cuco::xxhash_32<int64_t>>(0, 1076387279, 42);
  result[i++] = check_hash_result<cuco::xxhash_32<int64_t>>(42, 2332451213, 0);
  result[i++] = check_hash_result<cuco::xxhash_32<int64_t>>(123456789, 1561711919, 0);

#if defined(CUCO_HAS_INT128)
  result[i++] = check_hash_result<cuco::xxhash_32<__int128>>(123456789, 1846633701, 0);
#endif

  result[i++] = check_hash_result<cuco::xxhash_32<large_key<32>>>(123456789, 3715432378, 0);
}

TEST_CASE("utility cuco::xxhash_32 test", "")
{
  // Reference hash values were computed using https://github.com/Cyan4973/xxHash
  SECTION("Check if host-generated hash values match the reference implementation.")
  {
    CHECK(check_hash_result<cuco::xxhash_32<char>>(0, 3479547966, 0));
    CHECK(check_hash_result<cuco::xxhash_32<char>>(42, 3774771295, 0));
    CHECK(check_hash_result<cuco::xxhash_32<char>>(0, 2099223482, 42));

    CHECK(check_hash_result<cuco::xxhash_32<int32_t>>(0, 148298089, 0));
    CHECK(check_hash_result<cuco::xxhash_32<int32_t>>(0, 2132181312, 42));
    CHECK(check_hash_result<cuco::xxhash_32<int32_t>>(42, 1161967057, 0));
    CHECK(check_hash_result<cuco::xxhash_32<int32_t>>(123456789, 2987034094, 0));

    CHECK(check_hash_result<cuco::xxhash_32<int64_t>>(0, 3736311059, 0));
    CHECK(check_hash_result<cuco::xxhash_32<int64_t>>(0, 1076387279, 42));
    CHECK(check_hash_result<cuco::xxhash_32<int64_t>>(42, 2332451213, 0));
    CHECK(check_hash_result<cuco::xxhash_32<int64_t>>(123456789, 1561711919, 0));

#if defined(CUCO_HAS_INT128)
    CHECK(check_hash_result<cuco::xxhash_32<__int128>>(123456789, 1846633701, 0));
#endif

    // 32*4=128-byte key to test the pipelined outermost hashing loop
    CHECK(check_hash_result<cuco::xxhash_32<large_key<32>>>(123456789, 3715432378, 0));
  }

  SECTION("Check if device-generated hash values match the reference implementation.")
  {
    thrust::device_vector<bool> result(20, true);

    check_hash_result_kernel_32<<<1, 1>>>(result.begin());

    CHECK(cuco::test::all_of(result.begin(), result.end(), cuda::std::identity{}));
  }
}

TEMPLATE_TEST_CASE_SIG("utility hasher compute_hash tests",
                       "",
                       ((typename Hash), Hash),
                       (cuco::murmurhash3_32<char>),
                       (cuco::murmurhash3_32<int32_t>),
                       (cuco::xxhash_32<char>),
                       (cuco::xxhash_32<int32_t>),
                       (cuco::xxhash_64<char>),
                       (cuco::xxhash_64<int32_t>))
{
  using key_type = typename Hash::argument_type;

  Hash hash;
  key_type key = 42;

  SECTION("Identical keys with static and dynamic key size should have the same hash value.")
  {
    CHECK(hash(key) ==
          hash.compute_hash(reinterpret_cast<cuda::std::byte const*>(&key), sizeof(key_type)));
  }
}

template <typename OutputIter>
__global__ void check_murmurhash3_128_result_kernel(OutputIter result)
{
  int i = 0;

  result[i++] = check_hash_result<cuco::murmurhash3_x64_128<int32_t>, uint64_t>(
    0, {14961230494313510588ull, 6383328099726337777ull}, 0);
  result[i++] = check_hash_result<cuco::murmurhash3_x64_128<int32_t>, uint64_t>(
    9, {1779292183511753683ull, 16298496441448380334ull}, 0);
  result[i++] = check_hash_result<cuco::murmurhash3_x64_128<int32_t>, uint64_t>(
    42, {2913627637088662735ull, 16344193523890567190ull}, 0);
  result[i++] = check_hash_result<cuco::murmurhash3_x64_128<int32_t>, uint64_t>(
    42, {2248879576374326886ull, 18006515275339376488ull}, 42);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int32_t, 2>>, uint64_t>(
      {2, 2}, {12221386834995143465ull, 6690950894782946573ull}, 0);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int32_t, 3>>, uint64_t>(
      {1, 4, 9}, {299140022350411792ull, 9891903873182035274ull}, 42);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int32_t, 4>>, uint64_t>(
      {42, 64, 108, 1024}, {4333511168876981289ull, 4659486988434316416ull}, 63);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int32_t, 16>>, uint64_t>(
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {3302412811061286680ull, 7070355726356610672ull},
      1024);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int64_t, 2>>, uint64_t>(
      {2, 2}, {8554944597931919519ull, 14938998000509429729ull}, 0);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int64_t, 3>>, uint64_t>(
      {1, 4, 9}, {13442629947720186435ull, 7061727494178573325ull}, 42);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int64_t, 4>>, uint64_t>(
      {42, 64, 108, 1024}, {8786399719555989948ull, 14954183901757012458ull}, 63);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int64_t, 16>>, uint64_t>(
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {15409921801541329777ull, 10546487400963404004ull},
      1024);

  result[i++] = check_hash_result<cuco::murmurhash3_x86_128<int32_t>, uint32_t>(
    0, {3422973727, 2656139328, 2656139328, 2656139328}, 0);
  result[i++] = check_hash_result<cuco::murmurhash3_x86_128<int32_t>, uint32_t>(
    9, {2808089785, 314604614, 314604614, 314604614}, 0);
  result[i++] = check_hash_result<cuco::murmurhash3_x86_128<int32_t>, uint32_t>(
    42, {3611919118, 1962256489, 1962256489, 1962256489}, 0);
  result[i++] = check_hash_result<cuco::murmurhash3_x86_128<int32_t>, uint32_t>(
    42, {3399017053, 732469929, 732469929, 732469929}, 42);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int32_t, 2>>, uint32_t>(
      {2, 2}, {1234494082, 1431451587, 431049201, 431049201}, 0);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int32_t, 3>>, uint32_t>(
      {1, 4, 9}, {2516796247, 2757675829, 778406919, 2453259553}, 42);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int32_t, 4>>, uint32_t>(
      {42, 64, 108, 1024}, {2686265656, 591236665, 3797082165, 2731908938}, 63);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int32_t, 16>>, uint32_t>(
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {3918256832, 4205523739, 1707810111, 1625952473},
      1024);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int64_t, 2>>, uint32_t>(
      {2, 2}, {3811075945, 727160712, 3510740342, 235225510}, 0);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int64_t, 3>>, uint32_t>(
      {1, 4, 9}, {2817194959, 206796677, 3391242768, 248681098}, 42);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int64_t, 4>>, uint32_t>(
      {42, 64, 108, 1024}, {2335912146, 1566515912, 760710030, 452077451}, 63);
  result[i++] =
    check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int64_t, 16>>, uint32_t>(
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {1101169764, 1758958147, 2406511780, 2903571412},
      1024);
}

TEST_CASE("utility cuco::murmurhash3_x64_128 test", "")
{
  // Reference hash values were computed using
  // https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp

  SECTION("Check if host-generated hash values match the reference implementation.")
  {
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<int32_t>, uint64_t>(
      0, {14961230494313510588ull, 6383328099726337777ull}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<int32_t>, uint64_t>(
      9, {1779292183511753683ull, 16298496441448380334ull}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<int32_t>, uint64_t>(
      42, {2913627637088662735ull, 16344193523890567190ull}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<int32_t>, uint64_t>(
      42, {2248879576374326886ull, 18006515275339376488ull}, 42));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int32_t, 2>>, uint64_t>(
      {2, 2}, {12221386834995143465ull, 6690950894782946573ull}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int32_t, 3>>, uint64_t>(
      {1, 4, 9}, {299140022350411792ull, 9891903873182035274ull}, 42));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int32_t, 4>>, uint64_t>(
      {42, 64, 108, 1024}, {4333511168876981289ull, 4659486988434316416ull}, 63));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int32_t, 16>>, uint64_t>(
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {3302412811061286680ull, 7070355726356610672ull},
      1024));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int64_t, 2>>, uint64_t>(
      {2, 2}, {8554944597931919519ull, 14938998000509429729ull}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int64_t, 3>>, uint64_t>(
      {1, 4, 9}, {13442629947720186435ull, 7061727494178573325ull}, 42));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int64_t, 4>>, uint64_t>(
      {42, 64, 108, 1024}, {8786399719555989948ull, 14954183901757012458ull}, 63));
    CHECK(check_hash_result<cuco::murmurhash3_x64_128<cuda::std::array<int64_t, 16>>, uint64_t>(
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {15409921801541329777ull, 10546487400963404004ull},
      1024));

    CHECK(check_hash_result<cuco::murmurhash3_x86_128<int32_t>, uint32_t>(
      0, {3422973727, 2656139328, 2656139328, 2656139328}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<int32_t>, uint32_t>(
      9, {2808089785, 314604614, 314604614, 314604614}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<int32_t>, uint32_t>(
      42, {3611919118, 1962256489, 1962256489, 1962256489}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<int32_t>, uint32_t>(
      42, {3399017053, 732469929, 732469929, 732469929}, 42));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int32_t, 2>>, uint32_t>(
      {2, 2}, {1234494082, 1431451587, 431049201, 431049201}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int32_t, 3>>, uint32_t>(
      {1, 4, 9}, {2516796247, 2757675829, 778406919, 2453259553}, 42));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int32_t, 4>>, uint32_t>(
      {42, 64, 108, 1024}, {2686265656, 591236665, 3797082165, 2731908938}, 63));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int32_t, 16>>, uint32_t>(
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {3918256832, 4205523739, 1707810111, 1625952473},
      1024));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int64_t, 2>>, uint32_t>(
      {2, 2}, {3811075945, 727160712, 3510740342, 235225510}, 0));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int64_t, 3>>, uint32_t>(
      {1, 4, 9}, {2817194959, 206796677, 3391242768, 248681098}, 42));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int64_t, 4>>, uint32_t>(
      {42, 64, 108, 1024}, {2335912146, 1566515912, 760710030, 452077451}, 63));
    CHECK(check_hash_result<cuco::murmurhash3_x86_128<cuda::std::array<int64_t, 16>>, uint32_t>(
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {1101169764, 1758958147, 2406511780, 2903571412},
      1024));
  }

  SECTION("Check if device-generated hash values match the reference implementation.")
  {
    thrust::device_vector<bool> result(24, true);

    check_murmurhash3_128_result_kernel<<<1, 1>>>(result.begin());

    CHECK(cuco::test::all_of(result.begin(), result.end(), cuda::std::identity{}));
  }
}
