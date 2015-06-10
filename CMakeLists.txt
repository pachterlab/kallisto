cmake_minimum_required(VERSION 2.8.12)

project(kallisto)

set(EXT_PROJECTS_DIR ${PROJECT_SOURCE_DIR}/ext)

#add_compile_options(-Wall -Wno-unused-function)
add_compile_options(-std=c++11)
IF(CMAKE_BUILD_TYPE MATCHES Debug)
    message("debug mode")
		add_compile_options(-O0 -g)
ELSE(CMAKE_BUILD_TYPE MATCHES Debug)
		add_compile_options(-O3)
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

if(CMAKE_BUILD_TYPE MATCHES Profile)
	add_compile_options(-g)
endif(CMAKE_BUILD_TYPE MATCHES Profile)

if(LINK MATCHES static)
    message("static build")
ELSE(LINK MATCHES shared)
    message("shared build")
ENDIF(LINK MATCHES static)

# add_compile_options(-Wdeprecated-register)

add_subdirectory(src)
include_directories(${EXT_PROJECTS_DIR})

if (BUILD_TESTING)
    add_subdirectory(${EXT_PROJECTS_DIR}/catch)

    # Includes Catch in the project:
    include_directories(${CATCH_INCLUDE_DIR} ${COMMON_INCLUDES})
    enable_testing(true)  # Enables unit-testing.

    add_subdirectory(unit_tests)
endif(BUILD_TESTING)

# enable_testing()
# add_test(MainTest test/tests)
