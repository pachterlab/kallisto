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

#include <benchmark_defaults.hpp>
#include <benchmark_utils.hpp>

#include <cuco/static_multimap.cuh>
#include <cuco/utility/key_generator.cuh>

#include <nvbench/nvbench.cuh>

#include <thrust/device_vector.h>
#include <thrust/transform.h>

using namespace cuco::benchmark;  // defaults, dist_from_state
using namespace cuco::utility;    // key_generator, distribution

/**
 * @brief A benchmark evaluating `cuco::static_multimap::count` performance
 */
template <typename Key, typename Value, typename Dist>
std::enable_if_t<(sizeof(Key) == sizeof(Value)), void> static_multimap_count(
  nvbench::state& state, nvbench::type_list<Key, Value, Dist>)
{
  using pair_type = cuco::pair<Key, Value>;

  auto const num_keys      = state.get_int64("NumInputs");
  auto const occupancy     = state.get_float64("Occupancy");
  auto const matching_rate = state.get_float64("MatchingRate");

  std::size_t const size = num_keys / occupancy;

  thrust::device_vector<Key> keys(num_keys);

  key_generator gen{};
  gen.generate(dist_from_state<Dist>(state), keys.begin(), keys.end());

  thrust::device_vector<pair_type> pairs(num_keys);
  thrust::transform(keys.begin(), keys.end(), pairs.begin(), [] __device__(Key const& key) {
    return pair_type{key, {}};
  });

  gen.dropout(keys.begin(), keys.end(), matching_rate);

  state.add_element_count(num_keys);

  cuco::experimental::static_multimap<Key, Value> map{
    size, cuco::empty_key<Key>{-1}, cuco::empty_value<Value>{-1}};
  map.insert(pairs.begin(), pairs.end());

  state.exec(nvbench::exec_tag::sync, [&](nvbench::launch& launch) {
    auto count = map.count(keys.begin(), keys.end(), {launch.get_stream()});
  });
}

template <typename Key, typename Value, typename Dist>
std::enable_if_t<(sizeof(Key) != sizeof(Value)), void> static_multimap_count(
  nvbench::state& state, nvbench::type_list<Key, Value, Dist>)
{
  state.skip("Key should be the same type as Value.");
}

NVBENCH_BENCH_TYPES(static_multimap_count,
                    NVBENCH_TYPE_AXES(defaults::KEY_TYPE_RANGE,
                                      defaults::VALUE_TYPE_RANGE,
                                      nvbench::type_list<distribution::uniform>))
  .set_name("static_multimap_count_uniform_capacity")
  .set_type_axes_names({"Key", "Value", "Distribution"})
  .add_int64_axis("NumInputs", defaults::N_RANGE_CACHE)
  .add_float64_axis("Occupancy", {defaults::OCCUPANCY})
  .add_int64_axis("Multiplicity", {defaults::MULTIPLICITY})
  .add_float64_axis("MatchingRate", {defaults::MATCHING_RATE});

NVBENCH_BENCH_TYPES(static_multimap_count,
                    NVBENCH_TYPE_AXES(defaults::KEY_TYPE_RANGE,
                                      defaults::VALUE_TYPE_RANGE,
                                      nvbench::type_list<distribution::uniform>))
  .set_name("static_multimap_count_uniform_occupancy")
  .set_type_axes_names({"Key", "Value", "Distribution"})
  .add_int64_axis("NumInputs", {defaults::N})
  .add_float64_axis("Occupancy", defaults::OCCUPANCY_RANGE)
  .add_int64_axis("Multiplicity", {defaults::MULTIPLICITY})
  .add_float64_axis("MatchingRate", {defaults::MATCHING_RATE});

NVBENCH_BENCH_TYPES(static_multimap_count,
                    NVBENCH_TYPE_AXES(defaults::KEY_TYPE_RANGE,
                                      defaults::VALUE_TYPE_RANGE,
                                      nvbench::type_list<distribution::uniform>))
  .set_name("static_multimap_count_uniform_multiplicity")
  .set_type_axes_names({"Key", "Value", "Distribution"})
  .add_int64_axis("NumInputs", {defaults::N})
  .add_float64_axis("Occupancy", {defaults::OCCUPANCY})
  .add_int64_axis("Multiplicity", defaults::MULTIPLICITY_RANGE)
  .add_float64_axis("MatchingRate", {defaults::MATCHING_RATE});

NVBENCH_BENCH_TYPES(static_multimap_count,
                    NVBENCH_TYPE_AXES(defaults::KEY_TYPE_RANGE,
                                      defaults::VALUE_TYPE_RANGE,
                                      nvbench::type_list<distribution::uniform>))
  .set_name("static_multimap_count_uniform_matching_rate")
  .set_type_axes_names({"Key", "Value", "Distribution"})
  .add_int64_axis("NumInputs", {defaults::N})
  .add_float64_axis("Occupancy", {defaults::OCCUPANCY})
  .add_int64_axis("Multiplicity", {defaults::MULTIPLICITY})
  .add_float64_axis("MatchingRate", defaults::MATCHING_RATE_RANGE);
