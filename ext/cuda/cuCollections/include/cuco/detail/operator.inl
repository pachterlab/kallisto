/*
 * Copyright (c) 2022-2025, NVIDIA CORPORATION.
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

#pragma once

#include <cuco/utility/traits.hpp>

#include <type_traits>

namespace cuco {
namespace detail {

/**
 * @brief CRTP mixin which augments a given `Reference` with an `Operator`.
 *
 * @throw If the operator is not defined in `include/cuco/operator.hpp`
 *
 * @tparam Operator Operator type, i.e., `cuco::op::*_tag`
 * @tparam Reference The reference type.
 *
 * @note This primary template should never be instantiated.
 */
template <typename Operator, typename Reference>
class operator_impl {
  static_assert(cuco::dependent_false<Operator, Reference>,
                "Operator type is not supported by reference type.");
};

/**
 * @brief Checks if the given `Operator` is contained in a list of `Operators`.
 *
 * @tparam Operator Operator type, i.e., `cuco::op::*_tag`
 * @tparam Operators List of operators to search in
 *
 * @return `true` if `Operator` is contained in `Operators`, `false` otherwise.
 */
template <typename Operator, typename... Operators>
__host__ __device__ static constexpr bool has_operator()
{
  return ((std::is_same_v<Operators, Operator>) || ...);
}

}  // namespace detail
}  // namespace cuco
