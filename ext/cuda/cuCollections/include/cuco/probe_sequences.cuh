/*
 * Copyright (c) 2021-2025, NVIDIA CORPORATION.
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

#include <cuco/detail/probe_sequence_impl.cuh>

namespace cuco::legacy {

/**
 * @brief Public linear probing scheme class.
 *
 * Linear probing is efficient when few collisions are present. Performance hints:
 * - Use linear probing when collisions are rare. e.g. low occupancy or low multiplicity.
 * - `CGSize` = 1 or 2 when hash map is small (10'000'000 or less), 4 or 8 otherwise.
 *
 * `Hash` should be callable object type.
 *
 * @tparam CGSize Size of CUDA Cooperative Groups
 * @tparam Hash Unary callable type
 */
template <int CGSize, typename Hash>
class linear_probing : public detail::probe_sequence_base<CGSize> {
 public:
  using probe_sequence_base_type =
    detail::probe_sequence_base<CGSize>;  ///< The base probe scheme type
  using probe_sequence_base_type::cg_size;
  using probe_sequence_base_type::vector_width;

  /// Type of implementation details
  template <typename Key, typename Value, cuda::thread_scope Scope>
  using impl = detail::linear_probing_impl<Key, Value, Scope, vector_width(), CGSize, Hash>;
};

/**
 *
 * @brief Public double hashing scheme class.
 *
 * Default probe sequence for `cuco::static_multimap`. Double hashing shows superior
 * performance when dealing with high multiplicty and/or high occupancy use cases. Performance
 * hints:
 * - `CGSize` = 1 or 2 when hash map is small (10'000'000 or less), 4 or 8 otherwise.
 *
 * `Hash1` and `Hash2` should be callable object type.
 *
 * @tparam CGSize Size of CUDA Cooperative Groups
 * @tparam Hash1 Unary callable type
 * @tparam Hash2 Unary callable type
 */
template <int CGSize, typename Hash1, typename Hash2 = Hash1>
class double_hashing : public detail::probe_sequence_base<CGSize> {
 public:
  using probe_sequence_base_type =
    detail::probe_sequence_base<CGSize>;  ///< The base probe scheme type
  using probe_sequence_base_type::cg_size;
  using probe_sequence_base_type::vector_width;

  /// Type of implementation details
  template <typename Key, typename Value, cuda::thread_scope Scope>
  using impl = detail::double_hashing_impl<Key, Value, Scope, vector_width(), CGSize, Hash1, Hash2>;
};

}  // namespace cuco::legacy
