#=============================================================================
# Copyright (c) 2018-2025, NVIDIA CORPORATION.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=============================================================================
cmake_minimum_required(VERSION 3.23.1 FATAL_ERROR)

include(CTest)

###################################################################################################
# - Fetch Catch   ---------------------------------------------------------------------------------

CPMAddPackage(
  NAME Catch2
  GITHUB_REPOSITORY catchorg/Catch2
  VERSION 3.3.0
)

# Header for catch_discover_tests
if(Catch2_ADDED)
    include(${Catch2_SOURCE_DIR}/extras/Catch.cmake)
endif()

###################################################################################################
function(ConfigureTest TEST_NAME)
    add_executable(${TEST_NAME} ${ARGN})
    target_link_libraries(${TEST_NAME} PRIVATE Catch2::Catch2WithMain cuco CUDA::cudart)
    target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    set_target_properties(${TEST_NAME} PROPERTIES
                                       RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
    target_compile_options(${TEST_NAME} PRIVATE --compiler-options=-Wall --compiler-options=-Wextra
      --compiler-options=-Werror -Wno-deprecated-gpu-targets --expt-extended-lambda)
    # Add GCC-specific warning suppression only for GCC
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${TEST_NAME} PRIVATE -Xcompiler -Wno-subobject-linkage)
    endif()
    catch_discover_tests(${TEST_NAME} EXTRA_ARGS --allow-running-no-tests)
endfunction(ConfigureTest)

###################################################################################################
### test sources ##################################################################################
###################################################################################################

###################################################################################################
# - utility tests ---------------------------------------------------------------------------------
ConfigureTest(UTILITY_TEST
    utility/extent_test.cu
    utility/storage_test.cu
    utility/fast_int_test.cu
    utility/hash_test.cu
    utility/probing_scheme_test.cu)

###################################################################################################
# - static_set tests ------------------------------------------------------------------------------
ConfigureTest(STATIC_SET_TEST
    static_set/atomic_storage_test.cu
    static_set/capacity_test.cu
    static_set/for_each_test.cu
    static_set/heterogeneous_lookup_test.cu
    static_set/insert_and_find_test.cu
    static_set/large_input_test.cu
    static_set/retrieve_test.cu
    static_set/retrieve_all_test.cu
    static_set/retrieve_if_test.cu
    static_set/rehash_test.cu
    static_set/size_test.cu
    static_set/shared_memory_test.cu
    static_set/stream_test.cu
    static_set/unique_sequence_test.cu)

###################################################################################################
# - static_map tests ------------------------------------------------------------------------------
ConfigureTest(STATIC_MAP_TEST
    static_map/capacity_test.cu
    static_map/contains_test.cu
    static_map/custom_type_test.cu
    static_map/duplicate_keys_test.cu
    static_map/erase_test.cu
    static_map/find_test.cu
    static_map/for_each_test.cu
    static_map/hash_test.cu
    static_map/heterogeneous_lookup_test.cu
    static_map/insert_and_find_test.cu
    static_map/insert_or_assign_test.cu
    static_map/insert_or_apply_test.cu
    static_map/key_sentinel_test.cu
    static_map/shared_memory_test.cu
    static_map/stream_test.cu
    static_map/rehash_test.cu
    static_map/retrieve_test.cu
    static_map/retrieve_if_test.cu)

###################################################################################################
# - dynamic_map tests -----------------------------------------------------------------------------
ConfigureTest(DYNAMIC_MAP_TEST
    dynamic_map/unique_sequence_test.cu
    dynamic_map/unique_sequence_test_experimental.cu
    dynamic_map/erase_test.cu
    dynamic_map/find_test.cu
    dynamic_map/retrieve_all_test.cu)

###################################################################################################
# - static_multiset tests -------------------------------------------------------------------------
ConfigureTest(STATIC_MULTISET_TEST
    static_multiset/contains_test.cu
    static_multiset/count_test.cu
    static_multiset/custom_count_test.cu
    static_multiset/find_test.cu
    static_multiset/insert_test.cu
    static_multiset/for_each_test.cu
    static_multiset/retrieve_test.cu
    static_multiset/retrieve_if_test.cu
    static_multiset/large_input_test.cu
    static_multiset/load_factor_test.cu
    static_multiset/stream_test.cu)

###################################################################################################
# - static_multimap tests -------------------------------------------------------------------------
ConfigureTest(STATIC_MULTIMAP_TEST
    static_multimap/count_test.cu
    static_multimap/find_test.cu
    static_multimap/heterogeneous_lookup_test.cu
    static_multimap/insert_contains_test.cu
    static_multimap/insert_if_test.cu
    static_multimap/multiplicity_test.cu
    static_multimap/for_each_test.cu
    static_multimap/retrieve_test.cu
    static_multimap/retrieve_if_test.cu)

###################################################################################################
# - dynamic_bitset tests --------------------------------------------------------------------------
ConfigureTest(DYNAMIC_BITSET_TEST
    dynamic_bitset/find_next_test.cu
    dynamic_bitset/get_test.cu
    dynamic_bitset/rank_test.cu
    dynamic_bitset/select_test.cu
    dynamic_bitset/size_test.cu)

###################################################################################################
# - hyperloglog ----------------------------------------------------------------------
ConfigureTest(HYPERLOGLOG_TEST
    hyperloglog/unique_sequence_test.cu
    hyperloglog/spark_parity_test.cu
    hyperloglog/device_ref_test.cu)

###################################################################################################
# - bloom_filter ----------------------------------------------------------------------------------
ConfigureTest(BLOOM_FILTER_TEST
    bloom_filter/unique_sequence_test.cu
    bloom_filter/arrow_policy_test.cu
    bloom_filter/variable_cg_test.cu)

###################################################################################################
# - roaring_bitmap ---------------------------------------------------------------------------------
ConfigureTest(ROARING_BITMAP_TEST
    roaring_bitmap/contains_test.cu)
