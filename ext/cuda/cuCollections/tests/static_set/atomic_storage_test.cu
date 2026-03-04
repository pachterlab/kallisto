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

#include <cuco/static_set.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <catch2/catch_test_macros.hpp>

using T    = int32_t;
using Hash = uint32_t;
using Key  = cuco::pair<Hash, T>;

struct hasher {
  __device__ Hash operator()(Key const& k) const { return k.first; }
};

struct always_not_equal {
  __device__ constexpr bool operator()(Key const&, Key const&) const noexcept
  {
    // All build table keys are distinct thus `false` no matter what
    return false;
  }
};

class build_fn {
 public:
  __device__ __forceinline__ auto operator()(T i) const noexcept { return cuco::pair{_hash(i), i}; }

 private:
  cuco::default_hash_function<T> _hash{};
};

// This test exercise is designed to replicate a Spark runtime failure scenario
// https://github.com/NVIDIA/spark-rapids/issues/12586 and
// https://github.com/rapidsai/cudf/issues/18587
// that is not addressed by the current test suite. It will result in a runtime
// crash if the CCCL atomic storage is not managed correctly.
TEST_CASE("atomic_storage_test", "")
{
  using probe = cuco::linear_probing<1, hasher>;

  auto const num_keys = 100'000;

  auto set = cuco::static_set{cuco::extent<int>{num_keys},
                              0.5,
                              cuco::empty_key<Key>{Key{std::numeric_limits<Hash>::max(), -1}},
                              always_not_equal{},
                              probe{},
                              {},
                              cuco::storage<1>{}};

  auto keys_begin = thrust::make_transform_iterator(thrust::counting_iterator{0}, build_fn{});

  set.insert_async(keys_begin, keys_begin + num_keys);
  auto const count = set.size();

  REQUIRE(count == num_keys);
}
