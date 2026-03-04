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

set(rapids-cmake-version 25.12)
if(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/CUCO_RAPIDS.cmake)
    file(DOWNLOAD
      https://raw.githubusercontent.com/rapidsai/rapids-cmake/branch-${rapids-cmake-version}/RAPIDS.cmake
         ${CMAKE_CURRENT_BINARY_DIR}/CUCO_RAPIDS.cmake)
endif()
include(${CMAKE_CURRENT_BINARY_DIR}/CUCO_RAPIDS.cmake)

include(rapids-cmake)
include(rapids-cpm)
include(rapids-cuda)
include(rapids-export)
include(rapids-find)

# * Determine GPU architectures
# * Enable the CMake CUDA language
rapids_cuda_init_architectures(CUCO)

project(CUCO VERSION 0.0.1 LANGUAGES CXX CUDA)

###################################################################################################
# - build options ---------------------------------------------------------------------------------

set(default_build_option_state OFF)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_LIST_DIR}")
  set(default_build_option_state ON)
endif()

option(BUILD_TESTS "Configure CMake to build tests" ${default_build_option_state})
option(BUILD_BENCHMARKS "Configure CMake to build (google) benchmarks" ${default_build_option_state})
option(BUILD_EXAMPLES "Configure CMake to build examples" ${default_build_option_state})
option(INSTALL_CUCO "Enable CMake install rules for cuco" ${default_build_option_state})

# Write the version header
rapids_cmake_write_version_file(include/cuco/version_config.hpp)

##############################################################################
# - build type ---------------------------------------------------------------

# Set a default build type if none was specified
rapids_cmake_build_type(Release)

# needed for clangd and clang-tidy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

##############################################################################
# - compiler options ---------------------------------------------------------

# * Find CUDAToolkit package
# * Offers support for CMAKE_CUDA_ARCHITECTURES=NATIVE
rapids_find_package(
    CUDAToolkit REQUIRED
    BUILD_EXPORT_SET cuco-exports
    INSTALL_EXPORT_SET cuco-exports
)

###################################################################################################
# - find packages we depend on --------------------------------------------------------------------

rapids_cpm_init()

include(cmake/thirdparty/get_cccl.cmake)

###################################################################################################
# - cuco target   ---------------------------------------------------------------------------------
add_library(cuco INTERFACE)
add_library(cuco::cuco ALIAS cuco)
target_include_directories(cuco INTERFACE
    INTERFACE   $<BUILD_INTERFACE:${CUCO_SOURCE_DIR}/include>
                $<INSTALL_INTERFACE:include>)
target_link_libraries(cuco INTERFACE CCCL::CCCL CUDA::toolkit)
target_compile_features(cuco INTERFACE cxx_std_17 cuda_std_17)

###################################################################################################
# - Optionally download RoaringFormatSpec test data -----------------------------------------------

option(CUCO_DOWNLOAD_ROARING_TESTDATA "Download RoaringFormatSpec test data" ON)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/roaring_testdata.cmake)

###################################################################################################
# - optionally build tests ------------------------------------------------------------------------

if(BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
    # Include header tests as part of the regular test suite
    include(cmake/header_testing.cmake)
endif(BUILD_TESTS)

###################################################################################################
# - Optionally build nvbench benchmarks -----------------------------------------------------------

if(BUILD_BENCHMARKS)
    include(${rapids-cmake-dir}/cpm/nvbench.cmake)
    include(${rapids-cmake-dir}/cpm/package_override.cmake)
    rapids_cpm_nvbench(BUILD_STATIC)

    add_subdirectory(benchmarks)
endif(BUILD_BENCHMARKS)

###################################################################################################
# - Optionally build examples ---------------------------------------------------------------------

if(BUILD_EXAMPLES)
    add_subdirectory(examples)
endif(BUILD_EXAMPLES)


###################################################################################################
# - Install targets -------------------------------------------------------------------------------

install(TARGETS cuco EXPORT cuco-exports)

set(doc_string
    [=[
Provide targets for cuCollections.

cuCollections (cuco) is an open-source, header-only library of GPU-accelerated,
concurrent data structures.

Similar to how Thrust and CUB provide STL-like, GPU accelerated algorithms and
primitives, cuCollections provides STL-like concurrent data structures.
cuCollections is not a one-to-one, drop-in replacement for STL data structures
like std::unordered_map. Instead, it provides functionally similar data
structures tailored for efficient use with GPUs.

]=])

# build directory cuco-config generation
rapids_export(
    BUILD cuco
    EXPORT_SET cuco-exports
    GLOBAL_TARGETS cuco
    NAMESPACE cuco::
    DOCUMENTATION doc_string)

if(INSTALL_CUCO)
    install(DIRECTORY include/cuco/ DESTINATION include/cuco)
    install(FILES ${CUCO_BINARY_DIR}/include/cuco/version_config.hpp DESTINATION include/cuco)
    # install directory cuco-config generation
    rapids_export(
        INSTALL cuco
        EXPORT_SET cuco-exports
        GLOBAL_TARGETS cuco
        NAMESPACE cuco::
        DOCUMENTATION doc_string)
endif()
