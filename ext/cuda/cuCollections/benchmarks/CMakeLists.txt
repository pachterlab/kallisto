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

###################################################################################################
### compiler function #############################################################################
###################################################################################################

###################################################################################################
function(ConfigureBench BENCH_NAME)
    add_executable(${BENCH_NAME} ${ARGN})
    set_target_properties(${BENCH_NAME} PROPERTIES
                                        POSITION_INDEPENDENT_CODE ON
                                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/benchmarks")
    target_include_directories(${BENCH_NAME} PRIVATE
                                             "${CMAKE_CURRENT_SOURCE_DIR}")
    target_compile_options(${BENCH_NAME} PRIVATE --compiler-options=-Wall --compiler-options=-Wextra
      --compiler-options=-Werror -Wno-deprecated-gpu-targets --expt-extended-lambda -lineinfo)
    # Add GCC-specific warning suppression only for GCC
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${BENCH_NAME} PRIVATE -Xcompiler -Wno-subobject-linkage)
    endif()
    target_link_libraries(${BENCH_NAME} PRIVATE
                                        nvbench::main
                                        pthread
                                        cuco)
endfunction(ConfigureBench)

###################################################################################################
### benchmark sources #############################################################################
###################################################################################################

###################################################################################################
# - static_set benchmarks -------------------------------------------------------------------------
ConfigureBench(STATIC_SET_BENCH
  static_set/contains_bench.cu
  static_set/find_bench.cu
  static_set/insert_bench.cu
  static_set/retrieve_bench.cu
  static_set/retrieve_all_bench.cu
  static_set/size_bench.cu
  static_set/rehash_bench.cu)

###################################################################################################
# - static_map benchmarks -------------------------------------------------------------------------
ConfigureBench(STATIC_MAP_BENCH
  static_map/insert_bench.cu
  static_map/find_bench.cu
  static_map/contains_bench.cu
  static_map/erase_bench.cu
  static_map/insert_or_apply_bench.cu)

###################################################################################################
# - static_multiset benchmarks --------------------------------------------------------------------
ConfigureBench(STATIC_MULTISET_BENCH
  static_multiset/contains_bench.cu
  static_multiset/retrieve_bench.cu
  static_multiset/count_bench.cu
  static_multiset/find_bench.cu
  static_multiset/insert_bench.cu)

###################################################################################################
# - static_multimap benchmarks --------------------------------------------------------------------
ConfigureBench(STATIC_MULTIMAP_BENCH
  static_multimap/insert_bench.cu
  static_multimap/retrieve_bench.cu
  static_multimap/query_bench.cu
  static_multimap/count_bench.cu)

###################################################################################################
# - dynamic_map benchmarks ------------------------------------------------------------------------
ConfigureBench(DYNAMIC_MAP_BENCH
  dynamic_map/insert_bench.cu
  dynamic_map/find_bench.cu
  dynamic_map/contains_bench.cu
  dynamic_map/erase_bench.cu
  dynamic_map/retrieve_all_bench.cu)

###################################################################################################
# - hash function benchmarks ----------------------------------------------------------------------
ConfigureBench(HASH_FUNCTION_BENCH
  hash_function/hash_function_bench.cu)

###################################################################################################
# - hyperloglog benchmarks -----------------------------------------------------------
ConfigureBench(HYPERLOGLOG_BENCH
  hyperloglog/hyperloglog_bench.cu)

###################################################################################################
# - bloom_filter benchmarks -----------------------------------------------------------------------
ConfigureBench(BLOOM_FILTER_BENCH
  bloom_filter/add_bench.cu
  bloom_filter/contains_bench.cu)

###################################################################################################
# - roaring_bitmap benchmarks ---------------------------------------------------------------------
ConfigureBench(ROARING_BITMAP_BENCH
  roaring_bitmap/contains_bench.cu)