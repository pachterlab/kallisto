#=============================================================================
# Copyright (c) 2025, NVIDIA CORPORATION.
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

# For every public header, build a translation unit containing `#include <header>`
# to let the compiler try to figure out warnings in that header if it is not otherwise
# included in tests, and also to verify if the headers are modular enough.
# .inl files are not globbed for, because they are not supposed to be used as public
# entrypoints.

function(cuco_add_header_tests)
  file(GLOB_RECURSE headers
    RELATIVE "${CUCO_SOURCE_DIR}/include"
    CONFIGURE_DEPENDS
    "${CUCO_SOURCE_DIR}/include/cuco/*.cuh"
    "${CUCO_SOURCE_DIR}/include/cuco/*.hpp"
  )
  
  list(LENGTH headers headers_count)
  message(STATUS "Found ${headers_count} headers for testing")

  # List of headers that have known issues or are not meant to be included directly
  set(excluded_headers
    # Add any headers that should be excluded from testing here
    # Example: cuco/internal_header.cuh
  )
  
  # Remove excluded headers
  if(excluded_headers)
    list(REMOVE_ITEM headers ${excluded_headers})
  endif()

  foreach (header IN LISTS headers)
    # Create a safe target name by replacing path separators and dots
    string(REPLACE "/" "_" header_target_name "${header}")
    string(REPLACE "." "_" header_target_name "${header_target_name}")
    # Use a hash to ensure uniqueness in case of similar names
    string(MD5 header_hash "${header}")
    string(SUBSTRING "${header_hash}" 0 8 header_hash_short)
    set(headertest_target "cuco_header_${header_target_name}_${header_hash_short}")
    
    set(header_src "${CMAKE_CURRENT_BINARY_DIR}/headers/${headertest_target}/${header}.cu")
    
    # Create the directory if it doesn't exist
    get_filename_component(header_dir "${header_src}" DIRECTORY)
    file(MAKE_DIRECTORY "${header_dir}")
    
    # Write simple test file that includes the header
    file(WRITE "${header_src}" "#include <${header}>\nint main() { return 0; }\n")

    # Create executable test for this specific header
    add_executable(${headertest_target} ${header_src})
    target_link_libraries(${headertest_target} PRIVATE cuco::cuco CUDA::cudart)

    target_compile_options(${headertest_target} PRIVATE
      $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>
      --compiler-options=-Wall --compiler-options=-Wextra
      --compiler-options=-Werror -Wno-deprecated-gpu-targets
    )

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      target_compile_options(${headertest_target} PRIVATE -Xcompiler -Wno-subobject-linkage)
    endif()

    set_target_properties(${headertest_target} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/headers"
    )

    # Add as a CTest test
    add_test(NAME ${headertest_target} COMMAND ${headertest_target})
  endforeach()
endfunction()

cuco_add_header_tests()
