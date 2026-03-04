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
# - compiler function -----------------------------------------------------------------------------

function(ConfigureExample EXAMPLE_NAME EXAMPLE_SRC)
    add_executable(${EXAMPLE_NAME} "${EXAMPLE_SRC}")
    set_target_properties(${EXAMPLE_NAME} PROPERTIES
                                          RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/examples")
    target_include_directories(${EXAMPLE_NAME} PRIVATE
                                             "${CMAKE_CURRENT_SOURCE_DIR}")
    target_compile_options(${EXAMPLE_NAME} PRIVATE --compiler-options=-Wall --compiler-options=-Wextra
      --compiler-options=-Werror -Wno-deprecated-gpu-targets --expt-extended-lambda)
    # Add GCC-specific warning suppression only for GCC
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${EXAMPLE_NAME} PRIVATE -Xcompiler -Wno-subobject-linkage)
    endif()
    target_link_libraries(${EXAMPLE_NAME} PRIVATE cuco CUDA::cudart)
endfunction(ConfigureExample)

###################################################################################################
### Example sources ###############################################################################
###################################################################################################

ConfigureExample(STATIC_SET_HOST_BULK_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_set/host_bulk_example.cu")
ConfigureExample(STATIC_SET_DEVICE_REF_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_set/device_ref_example.cu")
ConfigureExample(STATIC_SET_DEVICE_SUBSETS_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_set/device_subsets_example.cu")
ConfigureExample(STATIC_SET_SHARED_MEMORY_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_set/shared_memory_example.cu")
ConfigureExample(STATIC_SET_MAPPING_TABLE_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_set/mapping_table_example.cu")
ConfigureExample(STATIC_MULTISET_HOST_BULK_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_multiset/host_bulk_example.cu")
ConfigureExample(STATIC_MAP_HOST_BULK_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_map/host_bulk_example.cu")
ConfigureExample(STATIC_MAP_DEVICE_SIDE_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_map/device_ref_example.cu")
ConfigureExample(STATIC_MAP_CUSTOM_TYPE_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_map/custom_type_example.cu")
ConfigureExample(STATIC_MAP_COUNT_BY_KEY_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_map/count_by_key_example.cu")
ConfigureExample(STATIC_MAP_PREALLOCATED_MEMORY_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_map/preallocated_memory_example.cu")
ConfigureExample(STATIC_MAP_HETEROGENEOUS_LOOKUP_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_map/heterogeneous_lookup_example.cu")
ConfigureExample(STATIC_MULTIMAP_HOST_BULK_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/static_multimap/host_bulk_example.cu")
ConfigureExample(HYPERLOGLOG_HOST_BULK_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/hyperloglog/host_bulk_example.cu")
ConfigureExample(HYPERLOGLOG_DEVICE_REF_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/hyperloglog/device_ref_example.cu")
ConfigureExample(BLOOM_FILTER_HOST_BULK_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/bloom_filter/host_bulk_example.cu")
ConfigureExample(ROARING_BITMAP_HOST_BULK_EXAMPLE "${CMAKE_CURRENT_SOURCE_DIR}/roaring_bitmap/host_bulk_example.cu")
