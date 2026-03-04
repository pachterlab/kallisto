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

#include <cuco/static_multiset.cuh>

#include <catch2/catch_template_test_macros.hpp>

using size_type = int32_t;

TEMPLATE_TEST_CASE_SIG(
  "static_multiset load factor tests",
  "",
  ((typename Key, cuco::test::probe_sequence Probe, int CGSize), Key, Probe, CGSize),
  (int32_t, cuco::test::probe_sequence::double_hashing, 1),
  (int32_t, cuco::test::probe_sequence::double_hashing, 2),
  (int64_t, cuco::test::probe_sequence::double_hashing, 1),
  (int64_t, cuco::test::probe_sequence::double_hashing, 2),
  (int32_t, cuco::test::probe_sequence::linear_probing, 1),
  (int32_t, cuco::test::probe_sequence::linear_probing, 2),
  (int64_t, cuco::test::probe_sequence::linear_probing, 1),
  (int64_t, cuco::test::probe_sequence::linear_probing, 2))
{
  constexpr size_type num_keys{10};

  using probe = std::conditional_t<Probe == cuco::test::probe_sequence::linear_probing,
                                   cuco::linear_probing<CGSize, cuco::default_hash_function<Key>>,
                                   cuco::double_hashing<CGSize, cuco::default_hash_function<Key>>>;

  SECTION("Negative load factor will throw exception")
  {
    REQUIRE_THROWS(cuco::static_multiset{
      num_keys, -0.1, cuco::empty_key<Key>{-1}, {}, probe{}, {}, cuco::storage<2>{}});
  }

  SECTION("Zero load factor will throw exception")
  {
    REQUIRE_THROWS(cuco::static_multiset{
      num_keys, 0.0, cuco::empty_key<Key>{-1}, {}, probe{}, {}, cuco::storage<2>{}});
  }

  SECTION("Load factor larger than one will throw exception")
  {
    REQUIRE_THROWS(cuco::static_multiset{
      num_keys, 1.1, cuco::empty_key<Key>{-1}, {}, probe{}, {}, cuco::storage<2>{}});
  }
}
