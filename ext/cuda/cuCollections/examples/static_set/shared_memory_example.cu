/*
 * Copyright (c) 2024-2025, NVIDIA CORPORATION.
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

#include <cuco/static_set.cuh>

#include <cuda/std/functional>

#include <cooperative_groups.h>

/**
 * @file shared_memory_example.cu
 * @brief Demonstrates usage of the static_set in shared memory.
 */

template <class SetRef>
__global__ void shmem_set_kernel(typename SetRef::extent_type valid_extent,
                                 cuco::empty_key<typename SetRef::key_type> empty_key_sentinel)
{
  // We first allocate the shared memory storage for the `set`.
  // The storage is comprised of contiguous buckets of slots,
  // which allow for vectorized loads.
  __shared__ typename SetRef::value_type slots[valid_extent.value()];

  // Next, we construct the actual storage object from the raw array.
  auto storage = SetRef::storage_ref_type(valid_extent, slots);
  // Now we can instantiate the set from the storage.
  auto set = SetRef(empty_key_sentinel, {}, {}, {}, storage);

  // The current thread blcok / CTA.
  auto const block = cooperative_groups::this_thread_block();

  // Initialize the raw storage using all threads in the block.
  set.initialize(block);
  // Synchronize the block to make sure initialization has finished.
  block.sync();

  // The `set` object does not come with any functionality. We first have to transform it into an
  // object that supports the function we need (in this case `insert`).
  auto insert_ref = set.rebind_operators(cuco::insert);

  // Each thread inserts its thread id into the set.
  typename SetRef::key_type const key = block.thread_rank();
  // Note that if you want to use a cg_size other then one, you have to use the cooperative
  // overload of this function, i.e., insert(cg, key);
  insert_ref.insert(key);

  // Synchronize the cta to make sure all insert operations have finished.
  block.sync();
  // Next, we want to check if the keys can be found again using the `contains` function. We create
  // a new non-owning object based on the `insert_ref` that supports `contains` but no longer
  // supports `insert`.
  // CAVEAT: concurrent use of `insert_ref` and `contains_ref` is undefined behavior.
  auto const contains_ref = insert_ref.rebind_operators(cuco::contains);

  // Check if all keys can be found
  if (not contains_ref.contains(key)) { printf("ERROR: Key %d not found\n", key); }
}

int main(void)
{
  using Key = int;

  // The "empty" sentinel is a reserved value required by our implementation.
  // Inserting or retrieving this value is UB.
  cuco::empty_key<Key> constexpr empty_key_sentinel{-1};
  // Width of vectorized loads during probing.
  auto constexpr bucket_size = 1;
  // Cooperative group size
  auto constexpr cg_size = 1;

  // Minimum number of slots in the set.
  // Static shared memory requires the size to be a compile time variable.
  // cuco::extent<int, 1000> is equivalent to `constexpr int capacity = 1000;`
  using extent_type = cuco::extent<int, 1000>;
  // Define the probing scheme.
  using probing_scheme_type = cuco::linear_probing<cg_size, cuco::default_hash_function<Key>>;

  // We define the set type given the parameters above.
  using set_type = cuco::static_set<Key,
                                    extent_type,
                                    cuda::thread_scope_block,
                                    cuda::std::equal_to<Key>,
                                    probing_scheme_type,
                                    cuco::cuda_allocator<Key>,
                                    cuco::storage<bucket_size>>;
  // Next, we can derive the non-owning reference type from the set type.
  // This is the type we use in the kernel to wrap a raw shared memory array as a `static_set`.
  using set_ref_type = typename set_type::ref_type<>;

  // Cuco imposes a number of non-trivial contraints on the capacity value.
  // This function will take the requested capacity (1000) and return the next larger
  // valid extent.
  auto constexpr valid_extent =
    cuco::make_valid_extent<probing_scheme_type, cuco::storage<bucket_size>>(extent_type{});

  // Launch the kernel with a single thread block.
  shmem_set_kernel<set_ref_type><<<1, 128>>>(valid_extent, empty_key_sentinel);
  cudaDeviceSynchronize();
}
