/*
 * Copyright (c) 2023-2024, NVIDIA CORPORATION.
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

#include <cuco/static_map.cuh>

#include <cuda/functional>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("static_map rehash test", "")
{
  using key_type    = int;
  using mapped_type = long;

  constexpr std::size_t num_keys{400};
  constexpr std::size_t num_erased_keys{100};

  cuco::static_map map{num_keys,
                       cuco::empty_key<key_type>{-1},
                       cuco::empty_value<mapped_type>{-1},
                       cuco::erased_key<key_type>{-2}};

  auto keys_begin = thrust::counting_iterator<key_type>(1);

  auto pairs_begin = thrust::make_transform_iterator(
    keys_begin,
    cuda::proclaim_return_type<cuco::pair<key_type, mapped_type>>([] __device__(key_type const& x) {
      return cuco::pair<key_type, mapped_type>(x, static_cast<mapped_type>(x));
    }));

  map.insert(pairs_begin, pairs_begin + num_keys);

  map.rehash();
  REQUIRE(map.size() == num_keys);

  map.rehash(num_keys * 2);
  REQUIRE(map.size() == num_keys);

  map.erase(keys_begin, keys_begin + num_erased_keys);
  map.rehash();
  REQUIRE(map.size() == num_keys - num_erased_keys);
}
