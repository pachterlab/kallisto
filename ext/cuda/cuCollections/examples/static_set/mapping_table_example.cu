/*
 * Copyright (c) 2024, NVIDIA CORPORATION.
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

#include <cuda/std/array>
#include <cuda/std/tuple>
#include <thrust/device_vector.h>

#include <iostream>

/**
 * @file mapping_table_example.cu
 *
 * @brief Demonstrates how to use hash set as a lookup table of the original data
 *
 * `cuco` hash tables such as `cuco::static_set` or `cuco::static_map` support only 4/8 byte keys.
 * This limitation arises because `cuco` hash tables rely on atomic Compare-And-Swap (CAS)
 * operations for key insertions (or queries), and the hardware natively supports only 4-byte and
 * 8-byte CAS. To enable support for larger keys, one approach is to implement atomic lock tables at
 * the software level. However, this approach would lead to a notable performance decrease due to
 * the high runtime cost of atomic lock tables.
 *
 * Additionally, `cuco` hash tables use open addressing as the hash collision resolution method.
 * This approach requires users to provide a sentinel that indicates unused slots in the data
 * structure. The sentinel value is a reserved value that must be never present in the problem. Note
 * that inserting or querying a sentinel value is undefined behavior. This can be problematic,
 * especially when the input data type is complex and determining a valid sentinel value is not
 * straightforward.
 *
 * This sample code demonstrates a solution to address these issues by using hash set as an
 * indirection mapping table to the original data:
 *  - The keys inserted in the hash table are indices of the original data array.
 *  - Using `-1` as a sentinel value is safe because accessing `data[-1]` is invalid.
 *  - Custom hashers and key equality comparators are required to hash and compare original keys
 *    based on indices.
 *
 * @note This example is for demonstration purposes only. It is not intended to show the most
 * performant way to do the example algorithm.
 */

/**
 * @brief User-defined key equal to compare two keys
 *
 * @tparam T Original key type
 */
template <typename T>
struct my_equal {
  my_equal(T const* data) : _data{data} {}
  /**
   * @brief Checks if two keys are identical based on their indices in the
   * original data array
   *
   * @param lhs The left hand side index
   * @param rhs The right hand side index
   * @return 'true' if two tuples are identical
   */
  __device__ constexpr bool operator()(int32_t lhs, int32_t rhs) const
  {
    // Check all 4 elements of a tuple to determine if two tuples are identical
    return cuda::std::get<0>(_data[lhs]) == cuda::std::get<0>(_data[rhs]) and
           cuda::std::get<1>(_data[lhs]) == cuda::std::get<1>(_data[rhs]) and
           cuda::std::get<2>(_data[lhs]) == cuda::std::get<2>(_data[rhs]) and
           cuda::std::get<3>(_data[lhs]) == cuda::std::get<3>(_data[rhs]);
  }
  T const* _data;
};

/**
 * @brief User-defined hash function to hash the original data based on its index
 *
 * @tparam T Original key type
 */
template <typename T>
struct my_hasher {
  my_hasher(T const* data) : _data{data} {}
  __device__ auto operator()(int32_t index) const
  {
    // Only hashes the first element of a tuple for demonstration purposes
    return cuda::std::get<0>(_data[index]);
  }
  T const* _data;
};

/**
 * @brief Utility to print the content of a given `tuple`
 *
 * @tparam T Type of the tuple
 */
template <typename T>
void print(T const& tuple)
{
  std::cout << "[" << cuda::std::get<0>(tuple) << ", " << cuda::std::get<1>(tuple) << ", "
            << cuda::std::get<2>(tuple) << ", "
            << "[" << cuda::std::get<3>(tuple)[0] << ", " << cuda::std::get<3>(tuple)[1] << ", "
            << cuda::std::get<3>(tuple)[2] << ", " << cuda::std::get<3>(tuple)[3] << "]]\n";
}

int main(void)
{
  // The original key type is larger than 8-byte and complex to spell the full type name
  using Key = cuda::std::tuple<uint32_t, char, bool, cuda::std::array<double, 4UL>>;
  // Imagine the array size is huge or the key type is more complex, it becomes impossible to
  // determine a valid sentinel value without introspecting the data
  auto const h_data =
    std::vector<Key>{cuda::std::tuple{11u, 'a', true, cuda::std::array{1., 2., 3., 4.}},
                     cuda::std::tuple{11u, 'a', true, cuda::std::array{1., 2., 3., 4.}},
                     cuda::std::tuple{22u, 'b', true, cuda::std::array{5., 6., 7., 8.}},
                     cuda::std::tuple{11u, 'a', true, cuda::std::array{5., 6., 7., 8.}},
                     cuda::std::tuple{11u, 'a', false, cuda::std::array{1., 2., 3., 4.}}};
  auto const size = h_data.size();
  thrust::device_vector<Key> d_data{h_data};

  // The actual key type is an index type, `int32_t` is large enough to cover the whole input range
  // and 4-byte atomic CAS is more efficient than the 8-byte one.
  using ActualKey = int32_t;
  // `-1` is a valid sentinel value since one will never access `data[-1]`
  ActualKey constexpr empty_key_sentinel = -1;

  auto const data_ptr = d_data.data().get();
  auto set = cuco::static_set{cuco::extent<std::size_t>{size * 2},  // about 50% load factor
                              cuco::empty_key{empty_key_sentinel},
                              my_equal{data_ptr},
                              cuco::linear_probing<1, my_hasher<Key>>{my_hasher<Key>{data_ptr}}};

  // The actual keys are indices of 5 elements
  auto const actual_keys = thrust::device_vector<ActualKey>{0, 1, 2, 3, 4};
  set.insert(actual_keys.begin(), actual_keys.end());

  auto unique_keys           = thrust::device_vector<ActualKey>(size);
  auto const unique_keys_end = set.retrieve_all(unique_keys.begin());
  auto const num             = std::distance(unique_keys.begin(), unique_keys_end);

  std::cout << "There are " << num << " distinct input elements:\n";
  for (auto i = 0; i < num; ++i) {
    // Retrieve query output based on indices
    print(h_data[unique_keys[i]]);
  }

  return 0;
}
