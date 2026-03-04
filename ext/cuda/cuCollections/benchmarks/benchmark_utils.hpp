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

#pragma once

#include <cuco/detail/error.hpp>
#include <cuco/utility/key_generator.cuh>

#include <nvbench/nvbench.cuh>

#include <thrust/iterator/iterator_traits.h>
#include <thrust/iterator/tabulate_output_iterator.h>

#include <nv/target>

namespace cuco::benchmark {

template <typename Dist>
auto dist_from_state(nvbench::state const& state)
{
  if constexpr (std::is_same_v<Dist, cuco::utility::distribution::unique>) {
    return Dist{};
  } else if constexpr (std::is_same_v<Dist, cuco::utility::distribution::uniform>) {
    auto const multiplicity = state.get_int64("Multiplicity");
    return Dist{multiplicity};
  } else if constexpr (std::is_same_v<Dist, cuco::utility::distribution::gaussian>) {
    auto const skew = state.get_float64("Skew");
    return Dist{skew};
  } else {
    CUCO_FAIL("Unexpected distribution type");
  }
}

template <typename T, typename NewType>
struct rebind_hasher;

template <template <typename> class Template, typename OldType, typename NewType>
struct rebind_hasher<Template<OldType>, NewType> {
  using type = Template<NewType>;
};

template <typename T, typename NewType>
using rebind_hasher_t = typename rebind_hasher<T, NewType>::type;

template <class OutputIt>
struct lazy_discard {
  OutputIt it;

  using index_type = typename thrust::iterator_traits<OutputIt>::difference_type;
  using value_type = typename thrust::iterator_traits<OutputIt>::value_type;

  __device__ void device_dispatch(index_type index, value_type const& value) const
  {
    // pick some predicate that is always false, but depends on the runtime value
    if (threadIdx.x > 2025 + *reinterpret_cast<char const*>(&value)) { *(it + index) = value; }
  }
  __host__ __device__ void operator()(index_type index, value_type const& value) const
  {
    NV_IF_TARGET(NV_IS_DEVICE,
                 this->device_dispatch(index, value);)  // we don't care about the host path for now
  }
};

/**
 * @brief An output iterator similar to `thrust::discard_iterator` but prevents the write from being
 * optimized out by the compiler.
 */
template <class OutputIt>
auto make_lazy_discard_iterator(OutputIt it)
{
  return thrust::tabulate_output_iterator(lazy_discard<OutputIt>{it});
}

}  // namespace cuco::benchmark

NVBENCH_DECLARE_TYPE_STRINGS(cuco::utility::distribution::unique, "UNIQUE", "distribution::unique");
NVBENCH_DECLARE_TYPE_STRINGS(cuco::utility::distribution::uniform,
                             "UNIFORM",
                             "distribution::uniform");
NVBENCH_DECLARE_TYPE_STRINGS(cuco::utility::distribution::gaussian,
                             "GAUSSIAN",
                             "distribution::gaussian");
