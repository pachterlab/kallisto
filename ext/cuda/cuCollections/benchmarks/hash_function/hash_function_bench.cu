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

#include <benchmark_defaults.hpp>

#include <cuco/hash_functions.cuh>
#include <cuco/utility/key_generator.cuh>

#include <nvbench/nvbench.cuh>

#include <cuda/std/cstddef>
#include <thrust/device_vector.h>

#include <cstdint>
#include <type_traits>

// repeat hash computation n times
static constexpr auto n_repeats = 100;

template <int32_t Words>
struct large_key {
  constexpr __host__ __device__ large_key(int32_t seed) noexcept
  {
    for (int32_t i = 0; i < Words; ++i) {
      data_[i] = seed;
    }
  }

 private:
  int32_t data_[Words];
};

template <typename T>
constexpr __host__ __device__ void hash_result_aggregate(T& agg, T hash_val)
{
  agg += hash_val;
}

template <>
constexpr __host__ __device__ void hash_result_aggregate(cuda::std::array<uint64_t, 2>& agg,
                                                         cuda::std::array<uint64_t, 2> hash_val)
{
  agg[0] += hash_val[0];
  agg[1] += hash_val[1];
}

template <>
constexpr __host__ __device__ void hash_result_aggregate(cuda::std::array<uint32_t, 4>& agg,
                                                         cuda::std::array<uint32_t, 4> hash_val)
{
  agg[0] += hash_val[0];
  agg[1] += hash_val[1];
  agg[2] += hash_val[2];
  agg[3] += hash_val[3];
}

template <int32_t BlockSize, typename Hasher, typename OutputIt>
__global__ void hash_bench_kernel(Hasher hash,
                                  cuco::detail::index_type n,
                                  OutputIt out,
                                  bool materialize_result)
{
  cuco::detail::index_type const gid         = BlockSize * blockIdx.x + threadIdx.x;
  cuco::detail::index_type const loop_stride = gridDim.x * BlockSize;
  cuco::detail::index_type idx               = gid;
  typename Hasher::result_type agg           = {};

  while (idx < n) {
    typename Hasher::argument_type key(idx);
    for (int32_t i = 0; i < n_repeats; ++i) {  // execute hash func n times
      hash_result_aggregate(agg, hash(key));
    }
    idx += loop_stride;
  }

  if (materialize_result) { out[gid] = agg; }
}

/**
 * @brief A benchmark evaluating performance of various hash functions
 */
template <typename Hash>
void hash_eval(nvbench::state& state, nvbench::type_list<Hash>)
{
  bool const materialize_result = false;
  constexpr auto block_size     = 128;
  auto const num_keys           = state.get_int64("NumInputs");
  auto const grid_size          = (num_keys + block_size * 16 - 1) / block_size * 16;

  thrust::device_vector<typename Hash::result_type> hash_values((materialize_result) ? num_keys
                                                                                     : 1);

  state.add_element_count(num_keys);

  state.exec([&](nvbench::launch& launch) {
    hash_bench_kernel<block_size><<<grid_size, block_size, 0, launch.get_stream()>>>(
      Hash{}, num_keys, hash_values.begin(), materialize_result);
  });
}

template <int32_t BlockSize, typename Hasher, typename InputIt, typename OutputIt>
__global__ void string_hash_bench_kernel(
  Hasher hash, InputIt in, cuco::detail::index_type n, OutputIt out, bool materialize_result)
{
  cuco::detail::index_type const gid         = BlockSize * blockIdx.x + threadIdx.x;
  cuco::detail::index_type const loop_stride = gridDim.x * BlockSize;
  cuco::detail::index_type idx               = gid;
  typename Hasher::result_type agg           = {};

  while (idx < n) {
    auto const key = thrust::raw_reference_cast(*(in + idx));
    for (int32_t i = 0; i < n_repeats; ++i) {  // execute hash func n times
      hash_result_aggregate(agg, hash.compute_hash(key.data(), key.size()));
    }
    idx += loop_stride;
  }

  if (materialize_result) { out[gid] = agg; }
}

/**
 * @brief A benchmark evaluating performance of various hash functions on random strings with
 * variable length and alignment
 */
template <typename Hash>
void string_hash_eval(nvbench::state& state, nvbench::type_list<Hash>)
{
  static_assert(std::is_same_v<typename Hash::argument_type, cuda::std::byte>,
                "Argument type must be cuda::std::byte");

  bool const materialize_result = false;
  constexpr auto block_size     = 128;
  auto const num_keys           = state.get_int64("NumInputs");
  auto const min_length         = state.get_int64("MinLength");
  auto const max_length         = state.get_int64("MaxLength");
  auto const grid_size          = (num_keys + block_size * 16 - 1) / block_size * 16;

  if (min_length > max_length) {
    state.skip("MinLength > MaxLength");
    return;
  }

  // auto const [keys, storage] = ... (can't capture structured bindings into lambdas in C++17)
  auto const sequences =
    cuco::utility::generate_random_byte_sequences(num_keys, min_length, max_length);
  auto const& keys = sequences.first;
  // auto const& storage = sequences.second;

  thrust::device_vector<typename Hash::result_type> hash_values((materialize_result) ? num_keys
                                                                                     : 1);

  state.add_element_count(num_keys);
  // state.add_global_memory_reads<cuda::std::byte>(storage.size() * n_repeats);

  state.exec([&](nvbench::launch& launch) {
    string_hash_bench_kernel<block_size><<<grid_size, block_size, 0, launch.get_stream()>>>(
      Hash{}, keys.begin(), num_keys, hash_values.begin(), materialize_result);
  });
}

NVBENCH_BENCH_TYPES(
  hash_eval,
  NVBENCH_TYPE_AXES(nvbench::type_list<cuco::murmurhash3_32<nvbench::int32_t>,
                                       cuco::murmurhash3_32<nvbench::int64_t>,
                                       cuco::murmurhash3_32<large_key<32>>,  // 32*4bytes
                                       cuco::xxhash_32<nvbench::int32_t>,
                                       cuco::xxhash_32<nvbench::int64_t>,
                                       cuco::xxhash_32<large_key<32>>,
                                       cuco::xxhash_64<nvbench::int32_t>,
                                       cuco::xxhash_64<nvbench::int64_t>,
                                       cuco::xxhash_64<large_key<32>>,
                                       cuco::murmurhash3_fmix_32<nvbench::int32_t>,
                                       cuco::murmurhash3_fmix_64<nvbench::int64_t>,
                                       cuco::murmurhash3_x86_128<nvbench::int32_t>,
                                       cuco::murmurhash3_x86_128<nvbench::int64_t>,
                                       cuco::murmurhash3_x86_128<large_key<32>>,
                                       cuco::murmurhash3_x64_128<nvbench::int32_t>,
                                       cuco::murmurhash3_x64_128<nvbench::int64_t>,
                                       cuco::murmurhash3_x64_128<large_key<32>>>))
  .set_name("hash_function_eval")
  .set_type_axes_names({"Hash"})
  .add_int64_axis("NumInputs", {cuco::benchmark::defaults::N * 10});

NVBENCH_BENCH_TYPES(
  string_hash_eval,
  NVBENCH_TYPE_AXES(nvbench::type_list<cuco::murmurhash3_32<cuda::std::byte>,
                                       cuco::xxhash_32<cuda::std::byte>,
                                       cuco::xxhash_64<cuda::std::byte>,
                                       cuco::murmurhash3_x86_128<cuda::std::byte>,
                                       cuco::murmurhash3_x64_128<cuda::std::byte>>))
  .set_name("string_hash_function_eval")
  .set_type_axes_names({"Hash"})
  .add_int64_axis("NumInputs", {cuco::benchmark::defaults::N / 4})
  .add_int64_axis("MinLength", {1, 4})
  .add_int64_axis("MaxLength", {4, 32, 64});
