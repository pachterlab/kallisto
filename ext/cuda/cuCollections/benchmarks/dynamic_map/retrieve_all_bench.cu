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

#include <benchmark_defaults.hpp>
#include <benchmark_utils.hpp>

#include <cuco/dynamic_map.cuh>
#include <cuco/utility/key_generator.cuh>

#include <nvbench/nvbench.cuh>

#include <thrust/device_vector.h>
#include <thrust/transform.h>

using namespace cuco::benchmark;  // defaults, dist_from_state
using namespace cuco::utility;    // key_generator, distribution

/**
 * @brief A benchmark evaluating `cuco::dynamic_map::retrieve_all` performance
 */
template <typename Key, typename Value, typename Dist>
std::enable_if_t<(sizeof(Key) == sizeof(Value)), void> dynamic_map_retrieve_all(
  nvbench::state& state, nvbench::type_list<Key, Value, Dist>)
{
  using pair_type = cuco::pair<Key, Value>;

  auto const num_keys     = state.get_int64("NumInputs");
  auto const initial_size = state.get_int64("InitSize");

  thrust::device_vector<Key> keys(num_keys);

  key_generator gen{};
  gen.generate(dist_from_state<Dist>(state), keys.begin(), keys.end());

  thrust::device_vector<pair_type> pairs(num_keys);
  thrust::transform(keys.begin(), keys.end(), pairs.begin(), [] __device__(Key const& key) {
    return pair_type(key, {});
  });

  cuco::dynamic_map<Key, Value> map{
    static_cast<size_t>(initial_size), cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{-1}};
  map.insert(pairs.begin(), pairs.end());
  // Prepare output buffers
  thrust::device_vector<Key> retrieved_keys(map.get_size());
  thrust::device_vector<Value> retrieved_values(map.get_size());

  state.add_element_count(map.get_size());

  state.exec(nvbench::exec_tag::sync, [&](nvbench::launch& launch) {
    map.retrieve_all(retrieved_keys.begin(), retrieved_values.begin(), launch.get_stream());
  });
}

template <typename Key, typename Value, typename Dist>
std::enable_if_t<(sizeof(Key) != sizeof(Value)), void> dynamic_map_retrieve_all(
  nvbench::state& state, nvbench::type_list<Key, Value, Dist>)
{
  state.skip("Key should be the same type as Value.");
}

NVBENCH_BENCH_TYPES(dynamic_map_retrieve_all,
                    NVBENCH_TYPE_AXES(defaults::KEY_TYPE_RANGE,
                                      defaults::VALUE_TYPE_RANGE,
                                      nvbench::type_list<distribution::unique>))
  .set_name("dynamic_map_retrieve_all_unique_capacity")
  .set_type_axes_names({"Key", "Value", "Distribution"})
  .add_int64_axis("NumInputs", defaults::N_RANGE)
  .add_int64_axis("InitSize", {defaults::INITIAL_SIZE});

NVBENCH_BENCH_TYPES(dynamic_map_retrieve_all,
                    NVBENCH_TYPE_AXES(defaults::KEY_TYPE_RANGE,
                                      defaults::VALUE_TYPE_RANGE,
                                      nvbench::type_list<distribution::unique>))
  .set_name("dynamic_map_retrieve_all_fixed_capacity")
  .set_type_axes_names({"Key", "Value", "Distribution"})
  .add_int64_axis("NumInputs", {defaults::N})
  .add_int64_axis("InitSize", {defaults::INITIAL_SIZE});
