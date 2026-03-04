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

#include <cuco/static_map.cuh>

#include <cuda/std/functional>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("static_map capacity test", "")
{
  using Key        = int32_t;
  using T          = int32_t;
  using ProbeT     = cuco::double_hashing<1, cuco::default_hash_function<Key>>;
  using Equal      = cuda::std::equal_to<Key>;
  using AllocatorT = cuco::cuda_allocator<cuda::std::byte>;
  using StorageT   = cuco::storage<2>;

  SECTION("zero capacity is allowed.")
  {
    auto constexpr gold_capacity = 4;

    using extent_type = cuco::extent<std::size_t, 0>;
    cuco::static_map<Key,
                     T,
                     extent_type,
                     cuda::thread_scope_device,
                     Equal,
                     ProbeT,
                     AllocatorT,
                     StorageT>
      map{extent_type{}, cuco::empty_key<Key>{-1}, cuco::empty_value<T>{-1}};
    auto const capacity = map.capacity();
    REQUIRE(capacity == gold_capacity);

    auto ref                = map.ref(cuco::insert);
    auto const ref_capacity = ref.capacity();
    REQUIRE(ref_capacity == gold_capacity);
  }

  SECTION("negative capacity (ikr -_-||) is also allowed.")
  {
    auto constexpr gold_capacity = 4;

    using extent_type = cuco::extent<int32_t>;
    cuco::static_map<Key,
                     T,
                     extent_type,
                     cuda::thread_scope_device,
                     Equal,
                     ProbeT,
                     AllocatorT,
                     StorageT>
      map{extent_type{-10}, cuco::empty_key<Key>{-1}, cuco::empty_value<T>{-1}};
    auto const capacity = map.capacity();
    REQUIRE(capacity == gold_capacity);

    auto ref                = map.ref(cuco::insert);
    auto const ref_capacity = ref.capacity();
    REQUIRE(ref_capacity == gold_capacity);
  }

  constexpr std::size_t num_keys{400};

  SECTION("Static bucket extent can be evaluated at build time.")
  {
    std::size_t constexpr gold_extent = 422;  // 211 x 2

    using extent_type = cuco::extent<std::size_t, num_keys>;
    cuco::static_map<Key,
                     T,
                     extent_type,
                     cuda::thread_scope_device,
                     Equal,
                     ProbeT,
                     AllocatorT,
                     StorageT>
      map{extent_type{}, cuco::empty_key<Key>{-1}, cuco::empty_value<T>{-1}};

    auto ref               = map.ref(cuco::insert);
    auto const num_buckets = ref.bucket_extent();
    STATIC_REQUIRE(static_cast<std::size_t>(num_buckets) == gold_extent);
  }

  SECTION("Dynamic extent is evaluated at run time.")
  {
    auto constexpr gold_capacity = 422;  // 211 x 2

    using extent_type = cuco::extent<std::size_t>;
    cuco::static_map<Key,
                     T,
                     extent_type,
                     cuda::thread_scope_device,
                     Equal,
                     ProbeT,
                     AllocatorT,
                     StorageT>
      map{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<T>{-1}};
    auto const capacity = map.capacity();
    REQUIRE(capacity == gold_capacity);

    auto ref                = map.ref(cuco::insert);
    auto const ref_capacity = ref.capacity();
    REQUIRE(ref_capacity == gold_capacity);
  }

  SECTION("map can be constructed from plain integer.")
  {
    auto constexpr gold_capacity = 422;  // 211 x 2

    cuco::static_map<Key,
                     T,
                     std::size_t,
                     cuda::thread_scope_device,
                     Equal,
                     ProbeT,
                     AllocatorT,
                     StorageT>
      map{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<T>{-1}};
    auto const capacity = map.capacity();
    REQUIRE(capacity == gold_capacity);

    auto ref                = map.ref(cuco::insert);
    auto const ref_capacity = ref.capacity();
    REQUIRE(ref_capacity == gold_capacity);
  }

  SECTION("map can be constructed from plain integer and load factor.")
  {
    auto constexpr gold_capacity = 502;  // 251 x 2

    cuco::static_map<Key,
                     T,
                     std::size_t,
                     cuda::thread_scope_device,
                     Equal,
                     ProbeT,
                     AllocatorT,
                     StorageT>
      map{num_keys, 0.8, cuco::empty_key<Key>{-1}, cuco::empty_value<T>{-1}};
    auto const capacity = map.capacity();
    REQUIRE(capacity == gold_capacity);

    auto ref                = map.ref(cuco::insert);
    auto const ref_capacity = ref.capacity();
    REQUIRE(ref_capacity == gold_capacity);
  }

  SECTION("Dynamic extent of linear probing is evaluated at run time.")
  {
    auto constexpr gold_capacity = 400;

    using probe = cuco::linear_probing<2, cuco::default_hash_function<Key>>;
    auto map =
      cuco::static_map<Key,
                       T,
                       cuco::extent<std::size_t>,
                       cuda::thread_scope_device,
                       Equal,
                       probe,
                       AllocatorT,
                       StorageT>{num_keys, cuco::empty_key<Key>{-1}, cuco::empty_value<T>{-1}};

    auto const capacity = map.capacity();
    REQUIRE(capacity == gold_capacity);

    auto ref                = map.ref(cuco::insert);
    auto const ref_capacity = ref.capacity();
    REQUIRE(ref_capacity == gold_capacity);
  }
}
