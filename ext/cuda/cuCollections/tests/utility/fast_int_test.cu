/*
 * Copyright (c) 2023, NVIDIA CORPORATION.
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

#include <cuco/utility/fast_int.cuh>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <cstdint>
#include <type_traits>

TEMPLATE_TEST_CASE(
  "utility::fast_int tests", "", std::int32_t, std::uint32_t, std::int64_t, std::uint64_t)
{
  TestType value           = GENERATE(1, 2, 9, 32, 4123, 8192, 4312456);
  TestType lhs             = GENERATE(1, 2, 9, 32, 4123, 8192, 4312456);
  constexpr auto max_value = std::numeric_limits<TestType>::max();

  cuco::utility::fast_int fast_value{value};

  SECTION("Should be explicitly convertible to the underlying integer type.")
  {
    REQUIRE(static_cast<TestType>(fast_value) == value);
  }

  SECTION("Fast div/mod should produce correct result.")
  {
    INFO(lhs << " /% " << value);
    REQUIRE(lhs / fast_value == lhs / value);
    REQUIRE(lhs % fast_value == lhs % value);
  }

  SECTION("Fast div/mod with maximum rhs value should produce correct result.")
  {
    INFO(lhs << " /% " << max_value);
    cuco::utility::fast_int fast_max{max_value};
    REQUIRE(lhs / fast_max == lhs / max_value);
    REQUIRE(lhs % fast_max == lhs % max_value);
  }

  SECTION("Fast div/mod with maximum lhs value should produce correct result.")
  {
    INFO(max_value << " /% " << value);
    REQUIRE(max_value / fast_value == max_value / value);
    REQUIRE(max_value % fast_value == max_value % value);
  }
}
