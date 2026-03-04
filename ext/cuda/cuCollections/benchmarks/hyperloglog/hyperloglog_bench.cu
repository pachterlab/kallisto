/*
 * Copyright (c) 2024-2025, NVIDIA CORPORATION.
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

#include <cuco/hyperloglog.cuh>
#include <cuco/static_set.cuh>
#include <cuco/utility/key_generator.cuh>

#include <nvbench/nvbench.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/iterator/transform_iterator.h>

#include <cmath>
#include <cstddef>

using namespace cuco::benchmark;  // defaults, dist_from_state
using namespace cuco::utility;    // key_generator, distribution

template <typename InputIt>
[[nodiscard]] std::size_t exact_distinct_count(InputIt first, std::size_t n)
{
  // TODO static_set currently only supports types up-to 8-bytes in size.
  // Casting is valid since the keys generated are representable in int64_t.
  using T = std::int64_t;

  auto cast_iter = thrust::make_transform_iterator(
    first, cuda::proclaim_return_type<T>([] __device__(auto i) { return static_cast<T>(i); }));

  auto set = cuco::static_set{n, 0.8, cuco::empty_key<T>{-1}};
  set.insert(cast_iter, cast_iter + n);
  return set.size();
}

template <class Estimator, class Dist>
[[nodiscard]] double relative_error(nvbench::state& state, std::size_t num_samples)
{
  using T = typename Estimator::value_type;

  auto const num_items      = state.get_int64("NumInputs");
  auto const sketch_size_kb = state.get_int64("SketchSizeKB");

  thrust::device_vector<T> items(num_items);

  key_generator gen{};
  Estimator estimator{cuco::sketch_size_kb(sketch_size_kb)};
  double error_sum = 0;
  for (std::size_t i = 0; i < num_samples; ++i) {
    gen.generate(dist_from_state<Dist>(state), items.begin(), items.end());
    estimator.add(items.begin(), items.end());
    double estimated_cardinality = estimator.estimate();
    double true_cardinality      = exact_distinct_count(items.begin(), num_items);
    error_sum += std::abs(estimated_cardinality / true_cardinality - 1.0);
    estimator.clear();
  }

  return error_sum / num_samples;
}

/**
 * @brief A benchmark evaluating `cuco::hyperloglog` end-to-end performance
 */
template <typename T, typename Dist>
void hyperloglog_e2e(nvbench::state& state, nvbench::type_list<T, Dist>)
{
  using estimator_type = cuco::hyperloglog<T>;

  auto const num_items      = state.get_int64("NumInputs");
  auto const sketch_size_kb = state.get_int64("SketchSizeKB");

  state.add_element_count(num_items);
  state.add_global_memory_reads<T>(num_items, "InputSize");

  auto const err_samples = (cuda::std::is_same_v<Dist, distribution::unique>) ? 1 : 5;
  auto const err         = relative_error<estimator_type, Dist>(state, err_samples);
  auto& summ             = state.add_summary("MeanRelativeError");
  summ.set_string("hint", "MRelErr");
  summ.set_string("short_name", "MeanRelativeError");
  summ.set_string("description", "Mean relatve approximation error.");
  summ.set_float64("value", err);

  thrust::device_vector<T> items(num_items);

  key_generator gen{};
  gen.generate(dist_from_state<Dist>(state), items.begin(), items.end());

  estimator_type estimator{cuco::sketch_size_kb(sketch_size_kb)};
  std::size_t estimated_cardinality = 0;
  state.exec(nvbench::exec_tag::sync | nvbench::exec_tag::timer,
             [&](nvbench::launch& launch, auto& timer) {
               timer.start();
               estimator.add_async(items.begin(), items.end(), {launch.get_stream()});
               estimated_cardinality = estimator.estimate({launch.get_stream()});
               timer.stop();

               estimator.clear_async({launch.get_stream()});
             });
}

/**
 * @brief A benchmark evaluating `cuco::hyperloglog::add_async` performance
 */
template <typename T, typename Dist>
void hyperloglog_add(nvbench::state& state, nvbench::type_list<T, Dist>)
{
  using estimator_type = cuco::hyperloglog<T>;

  auto const num_items      = state.get_int64("NumInputs");
  auto const sketch_size_kb = state.get_int64("SketchSizeKB");

  thrust::device_vector<T> items(num_items);

  key_generator gen{};
  gen.generate(dist_from_state<Dist>(state), items.begin(), items.end());

  state.add_element_count(num_items);
  state.add_global_memory_reads<T>(num_items, "InputSize");

  estimator_type estimator{cuco::sketch_size_kb(sketch_size_kb)};
  state.exec(nvbench::exec_tag::timer, [&](nvbench::launch& launch, auto& timer) {
    timer.start();
    estimator.add_async(items.begin(), items.end(), {launch.get_stream()});
    timer.stop();

    estimator.clear_async({launch.get_stream()});
  });
}

using TYPE_RANGE = nvbench::type_list<nvbench::int32_t, nvbench::int64_t, __int128_t>;

NVBENCH_BENCH_TYPES(hyperloglog_e2e,
                    NVBENCH_TYPE_AXES(TYPE_RANGE, nvbench::type_list<distribution::uniform>))
  .set_name("hyperloglog_e2e_uniform")
  .set_type_axes_names({"T", "Distribution"})
  .add_int64_power_of_two_axis("NumInputs", {28, 29, 30})
  .add_int64_axis("SketchSizeKB", {8, 16, 32, 64, 128, 256})  // 256KB uses gmem fallback kernel
  .add_int64_axis("Multiplicity", {1});

NVBENCH_BENCH_TYPES(hyperloglog_add,
                    NVBENCH_TYPE_AXES(TYPE_RANGE, nvbench::type_list<distribution::uniform>))
  .set_name("hyperloglog_add_uniform")
  .set_type_axes_names({"T", "Distribution"})
  .add_int64_power_of_two_axis("NumInputs", {28, 29, 30})
  .add_int64_axis("SketchSizeKB", {8, 16, 32, 64, 128, 256})
  .add_int64_axis("Multiplicity", {1});
