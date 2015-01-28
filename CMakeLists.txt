cmake_minimum_required(VERSION 2.8.8)

project(kallisto)

set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)

# Includes Catch in the project:
add_subdirectory(${EXT_PROJECTS_DIR}/catch)
include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
enable_testing(true)  # Enables unit-testing.

# add_compile_options(-Wall)
add_compile_options(-std=c++11)
add_compile_options(-Wdeprecated-register)

add_subdirectory(src)

add_subdirectory(test)
# enable_testing()
# add_test(MainTest test/tests)
