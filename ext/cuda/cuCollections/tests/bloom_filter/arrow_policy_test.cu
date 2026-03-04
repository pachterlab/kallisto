/*
 * Copyright (c) 2024, NVIDIA CORPORATION.
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

#include <cuco/bloom_filter.cuh>

#include <cuda/functional>
#include <thrust/device_vector.h>
#include <thrust/functional.h>

#include <catch2/catch_template_test_macros.hpp>

#include <random>
#include <type_traits>

namespace {

template <typename Key>
thrust::device_vector<uint32_t> get_arrow_filter_reference_bitset()
{
  static std::vector<thrust::device_vector<uint32_t>> const reference_bitsets{
    {
      3017764846,
      4219371383,
      4160077310,
      3214786543,
      4020088765,
      4294437885,
      2013200345,
      2550116063,
      855631359,
      4290436829,
      2884632042,
      1592483646,
      4281695998,
      2080111551,
      3220060030,
      4021279731,
    },  // type = int32, blocks = 2, num_keys = 100
    {
      860053560,  186397876,  1518788617, 2013987426, 545522943,  79856155,   103371656,
      20265733,   2168586373, 1210138712, 2437452036, 1342183988, 1107366672, 3560981000,
      2184221186, 1661010032, 2317009736, 1442875878, 1116227467, 3458613792, 114398528,
      679658134,  206734656,  340863450,  2220104352, 141846788,  948331524,  2344943952,
      4030989912, 3239203139, 2941256193, 4035057968,
    },  // type = int64, blocks = 4, num_keys = 50
    {
      3807057303, 3207519405, 2508188120, 1491024175, 2073585514, 2094743110, 2533287591,
      691662424,  1498889215, 2069126314, 2270481639, 796401059,  1961968732, 3512881027,
      3162306144, 2277085974, 3477648628, 1090385857, 4035761415, 1165385841, 4047856262,
      2297893848, 902599838,  418175153,  1437192944, 3673877288, 1536198910, 98677451,
      3620189521, 3794688342, 3625373537, 3550967313, 2119503598, 1805574667, 4076413870,
      2999897588, 3050286944, 4146882307, 3459690182, 167235913,  2078961096, 1863964920,
      1408130860, 4190644775, 532451008,  1563872186, 2529714129, 465761275,  3161649891,
      4204002248, 3931628891, 3251515903, 1421507581, 3849056446, 1748476671, 4223388125,
      1627644727, 2717076288, 2992639576, 3864567831, 190096788,  1885360347, 724608293,
      2768994330,
    }  // type = float, blocks = 8, num_keys = 200
  };

  if constexpr (std::is_same_v<Key, int32_t>) {
    return reference_bitsets[0];  // int32
  } else if constexpr (std::is_same_v<Key, int64_t>) {
    return reference_bitsets[1];  // int64
  } else if constexpr (std::is_same_v<Key, float>) {
    return reference_bitsets[2];  // float
  } else {
    throw std::invalid_argument("Reference bitsets available for int32, int64, float only.\n\n");
  }
}

template <typename Key>
std::pair<size_t, size_t> get_arrow_filter_test_settings()
{
  static std::vector<std::pair<size_t, size_t>> const test_settings = {
    {2, 100},  // type = int32, blocks = 2, num_keys = 100
    {4, 50},   // type = int64, blocks = 4, num_keys = 50
    {8, 200}   // type = float, blocks = 8, num_keys = 200
  };

  if constexpr (std::is_same_v<Key, int32_t>) {
    return test_settings[0];  // int32
  } else if constexpr (std::is_same_v<Key, int64_t>) {
    return test_settings[1];  // int64
  } else if constexpr (std::is_same_v<Key, float>) {
    return test_settings[2];  // float
  } else {
    throw std::invalid_argument("Test settings available for int32, int64, float only.\n\n");
  }
}

template <typename Key>
std::vector<Key> sequence_values(size_t size)
{
  std::vector<Key> values(size);
  std::iota(values.begin(), values.end(), Key{1});
  return values;
}

}  // namespace

template <typename Filter>
void test_filter_bitset(Filter& filter, size_t num_keys)
{
  using key_type  = typename Filter::key_type;
  using word_type = typename Filter::word_type;

  // Generate keys
  auto const h_keys = sequence_values<key_type>(num_keys);
  thrust::device_vector<key_type> d_keys(h_keys.begin(), h_keys.end());

  // Insert to the bloom filter
  filter.add(d_keys.begin(), d_keys.begin() + num_keys);

  // Get reference words device_vector
  auto const reference_words = get_arrow_filter_reference_bitset<key_type>();

  // Number of words in the filter
  auto const num_words = filter.block_extent() * filter.words_per_block;

  // Get the bitset
  thrust::device_vector<word_type> filter_words(filter.data(), filter.data() + num_words);

  REQUIRE(cuco::test::equal(
    filter_words.begin(),
    filter_words.end(),
    reference_words.begin(),
    cuda::proclaim_return_type<bool>([] __device__(auto const& filter_word, auto const& ref_word) {
      return filter_word == ref_word;
    })));
}

TEMPLATE_TEST_CASE_SIG("bloom_filter arrow filter policy bitset validation",
                       "",
                       (class Key),
                       (int32_t),
                       (int64_t),
                       (float))
{
  // Get test settings
  auto const [sub_filters, num_keys] = get_arrow_filter_test_settings<Key>();

  using policy_type = cuco::arrow_filter_policy<Key>;
  cuco::bloom_filter<Key, cuco::extent<size_t>, cuda::thread_scope_device, policy_type> filter{
    sub_filters};

  test_filter_bitset(filter, num_keys);
}
