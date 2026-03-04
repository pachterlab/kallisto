/*
 * Copyright (c) 2025, NVIDIA CORPORATION.
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

#include <cuco/static_map.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sequence.h>

#include <catch2/catch_template_test_macros.hpp>

using size_type = std::size_t;

template <class Container>
__global__ void test_retrieve_if_kernel(
  Container container_ref,
  typename Container::key_type* keys_begin,
  std::size_t num_keys,
  typename Container::key_type* stencil_begin,
  typename Container::key_type* output_probe,
  typename Container::value_type* output_match,
  cuda::atomic<int, cuda::thread_scope_device>* atomic_counter)
{
  using key_type = typename Container::key_type;
  namespace cg   = cooperative_groups;

  auto const block = cg::this_thread_block();
  auto const pred  = [] __device__(key_type k) { return k % 2 == 0; };

  container_ref.retrieve_if<128>(block,
                                 keys_begin,
                                 keys_begin + num_keys,
                                 stencil_begin,
                                 pred,
                                 output_probe,
                                 output_match,
                                 *atomic_counter);
}

template <class Container>
__global__ void test_retrieve_if_all_false_kernel(
  Container container_ref,
  typename Container::key_type* keys_begin,
  std::size_t num_keys,
  typename Container::key_type* stencil_begin,
  typename Container::key_type* output_probe,
  typename Container::value_type* output_match,
  cuda::atomic<int, cuda::thread_scope_device>* atomic_counter)
{
  using key_type = typename Container::key_type;
  namespace cg   = cooperative_groups;

  auto const block        = cg::this_thread_block();
  auto const always_false = [] __device__(key_type) { return false; };

  container_ref.retrieve_if<128>(block,
                                 keys_begin,
                                 keys_begin + num_keys,
                                 stencil_begin,
                                 always_false,
                                 output_probe,
                                 output_match,
                                 *atomic_counter);
}

template <class Container>
__global__ void test_retrieve_if_all_true_kernel(
  Container container_ref,
  typename Container::key_type* keys_begin,
  std::size_t num_keys,
  typename Container::key_type* stencil_begin,
  typename Container::key_type* output_probe,
  typename Container::value_type* output_match,
  cuda::atomic<int, cuda::thread_scope_device>* atomic_counter)
{
  using key_type = typename Container::key_type;
  namespace cg   = cooperative_groups;

  auto const block       = cg::this_thread_block();
  auto const always_true = [] __device__(key_type) { return true; };

  container_ref.retrieve_if<128>(block,
                                 keys_begin,
                                 keys_begin + num_keys,
                                 stencil_begin,
                                 always_true,
                                 output_probe,
                                 output_match,
                                 *atomic_counter);
}

TEMPLATE_TEST_CASE_SIG("static_map retrieve_if",
                       "",
                       ((typename Key, typename T), Key, T),
                       (int32_t, int32_t),
                       (int64_t, int64_t))
{
  constexpr size_type num_keys{400};

  using container_type = cuco::static_map<Key, T>;
  using value_type     = typename container_type::value_type;

  container_type container{num_keys * 2, cuco::empty_key<Key>{-1}, cuco::empty_value<T>{-1}};

  auto keys_begin  = thrust::counting_iterator<Key>(1);
  auto vals_begin  = thrust::counting_iterator<T>(1);
  auto pairs_begin = thrust::make_zip_iterator(thrust::make_tuple(keys_begin, vals_begin));

  container.insert(pairs_begin, pairs_begin + num_keys);

  SECTION("Testing retrieve_if with even predicate")
  {
    thrust::device_vector<Key> input_keys(keys_begin, keys_begin + num_keys);
    thrust::device_vector<Key> stencil_values(keys_begin, keys_begin + num_keys);
    thrust::device_vector<Key> probed_keys(num_keys);
    thrust::device_vector<value_type> matched_pairs(num_keys);

    cuda::atomic<int, cuda::thread_scope_device>* d_atomic_counter;
    CUCO_CUDA_TRY(
      cudaMalloc(&d_atomic_counter, sizeof(cuda::atomic<int, cuda::thread_scope_device>)));
    CUCO_CUDA_TRY(
      cudaMemset(d_atomic_counter, 0, sizeof(cuda::atomic<int, cuda::thread_scope_device>)));

    auto const container_ref = container.ref(cuco::op::retrieve);

    test_retrieve_if_kernel<<<1, 128>>>(container_ref,
                                        thrust::raw_pointer_cast(input_keys.data()),
                                        num_keys,
                                        thrust::raw_pointer_cast(stencil_values.data()),
                                        thrust::raw_pointer_cast(probed_keys.data()),
                                        thrust::raw_pointer_cast(matched_pairs.data()),
                                        d_atomic_counter);
    CUCO_CUDA_TRY(cudaDeviceSynchronize());

    int h_counter;
    CUCO_CUDA_TRY(cudaMemcpy(&h_counter, d_atomic_counter, sizeof(int), cudaMemcpyDeviceToHost));

    // Should retrieve even numbers only
    REQUIRE(h_counter > 0);
    REQUIRE(h_counter <= static_cast<int>(num_keys));

    CUCO_CUDA_TRY(cudaFree(d_atomic_counter));
  }

  SECTION("Testing retrieve_if with always false predicate")
  {
    thrust::device_vector<Key> input_keys(keys_begin, keys_begin + num_keys);
    thrust::device_vector<Key> stencil_values(keys_begin, keys_begin + num_keys);
    thrust::device_vector<Key> probed_keys(num_keys);
    thrust::device_vector<value_type> matched_pairs(num_keys);

    cuda::atomic<int, cuda::thread_scope_device>* d_atomic_counter;
    CUCO_CUDA_TRY(
      cudaMalloc(&d_atomic_counter, sizeof(cuda::atomic<int, cuda::thread_scope_device>)));
    CUCO_CUDA_TRY(
      cudaMemset(d_atomic_counter, 0, sizeof(cuda::atomic<int, cuda::thread_scope_device>)));

    auto const container_ref = container.ref(cuco::op::retrieve);

    test_retrieve_if_all_false_kernel<<<1, 128>>>(container_ref,
                                                  thrust::raw_pointer_cast(input_keys.data()),
                                                  num_keys,
                                                  thrust::raw_pointer_cast(stencil_values.data()),
                                                  thrust::raw_pointer_cast(probed_keys.data()),
                                                  thrust::raw_pointer_cast(matched_pairs.data()),
                                                  d_atomic_counter);
    CUCO_CUDA_TRY(cudaDeviceSynchronize());

    int h_counter;
    CUCO_CUDA_TRY(cudaMemcpy(&h_counter, d_atomic_counter, sizeof(int), cudaMemcpyDeviceToHost));

    // Should retrieve nothing
    REQUIRE(h_counter == 0);

    CUCO_CUDA_TRY(cudaFree(d_atomic_counter));
  }

  SECTION("Testing retrieve_if with always true predicate")
  {
    thrust::device_vector<Key> input_keys(keys_begin, keys_begin + num_keys);
    thrust::device_vector<Key> stencil_values(keys_begin, keys_begin + num_keys);
    thrust::device_vector<Key> probed_keys(num_keys);
    thrust::device_vector<value_type> matched_pairs(num_keys);

    cuda::atomic<int, cuda::thread_scope_device>* d_atomic_counter;
    CUCO_CUDA_TRY(
      cudaMalloc(&d_atomic_counter, sizeof(cuda::atomic<int, cuda::thread_scope_device>)));
    CUCO_CUDA_TRY(
      cudaMemset(d_atomic_counter, 0, sizeof(cuda::atomic<int, cuda::thread_scope_device>)));

    auto const container_ref = container.ref(cuco::op::retrieve);

    test_retrieve_if_all_true_kernel<<<1, 128>>>(container_ref,
                                                 thrust::raw_pointer_cast(input_keys.data()),
                                                 num_keys,
                                                 thrust::raw_pointer_cast(stencil_values.data()),
                                                 thrust::raw_pointer_cast(probed_keys.data()),
                                                 thrust::raw_pointer_cast(matched_pairs.data()),
                                                 d_atomic_counter);
    CUCO_CUDA_TRY(cudaDeviceSynchronize());

    int h_counter;
    CUCO_CUDA_TRY(cudaMemcpy(&h_counter, d_atomic_counter, sizeof(int), cudaMemcpyDeviceToHost));

    // Should retrieve all keys that exist in the container
    REQUIRE(h_counter == static_cast<int>(num_keys));

    CUCO_CUDA_TRY(cudaFree(d_atomic_counter));
  }
}