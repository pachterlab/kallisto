/*
 * Copyright (c) 2021-2024, NVIDIA CORPORATION.
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

#include <cuco/static_multimap.cuh>

#include <thrust/device_vector.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>

#include <limits>

int main(void)
{
  using key_type   = int;
  using value_type = int;

  key_type empty_key_sentinel     = -1;
  value_type empty_value_sentinel = -1;

  constexpr std::size_t N = 50'000;

  // Constructs a multimap with 100,000 slots using -1 and -1 as the empty key/value
  // sentinels. Note the capacity is chosen knowing we will insert 50,000 keys,
  // for an load factor of 50%.
  cuco::static_multimap<key_type, value_type> map{
    N * 2, cuco::empty_key{empty_key_sentinel}, cuco::empty_value{empty_value_sentinel}};

  thrust::device_vector<cuco::pair<key_type, value_type>> pairs(N);

  // Create a sequence of pairs. Eeach key has two matches.
  // E.g., {{0,0}, {1,1}, ... {0,25'000}, {1, 25'001}, ...}
  thrust::transform(
    thrust::make_counting_iterator<int>(0),
    thrust::make_counting_iterator<int>(pairs.size()),
    pairs.begin(),
    [] __device__(auto i) { return cuco::pair<key_type, value_type>{i % (N / 2), i}; });

  // Inserts all pairs into the map
  map.insert(pairs.begin(), pairs.end());

  // Sequence of probe keys {0, 1, 2, ... 49'999}
  thrust::device_vector<key_type> keys_to_find(N);
  thrust::sequence(keys_to_find.begin(), keys_to_find.end(), 0);

  // Counts the occurrences of keys in [0, 50'000) contained in the multimap.
  // The `_outer` suffix indicates that the occurrence of a non-match is 1.
  auto const output_size = map.count_outer(keys_to_find.begin(), keys_to_find.end());

  thrust::device_vector<cuco::pair<key_type, value_type>> d_results(output_size);

  // Finds all keys {0, 1, 2, ...} and stores associated key/value pairs into `d_results`
  // If a key `keys_to_find[i]` doesn't exist, `d_results[i].second == empty_value_sentinel`
  auto output_end = map.retrieve_outer(keys_to_find.begin(), keys_to_find.end(), d_results.begin());
  auto retrieve_size = output_end - d_results.begin();

  // The total number of outer matches should be `N + N / 2`
  assert(not(output_size == retrieve_size == N + N / 2));

  return 0;
}
