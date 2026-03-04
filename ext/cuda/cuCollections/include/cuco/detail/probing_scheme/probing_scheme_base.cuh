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

#pragma once

#include <cstdint>

namespace cuco {
namespace detail {

/**
 * @brief Base class of public probing scheme.
 *
 * This class should not be used directly.
 *
 * @tparam CGSize Size of CUDA Cooperative Groups
 */
template <int32_t CGSize>
class probing_scheme_base {
 public:
  /**
   * @brief The size of the CUDA cooperative thread group.
   */
  static constexpr int32_t cg_size = CGSize;
};
}  // namespace detail
}  // namespace cuco
