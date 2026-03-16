file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)
# Note: gpumain.cu is not in sources since it's a .cu file, not .cpp

include_directories(../ext/htslib)
include_directories(../ext/cuda/cuCollections/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

if(NOT MAX_KMER_SIZE)
set(MAX_KMER_SIZE "32")
endif()

add_compile_options(-Wno-subobject-linkage) # Suppress bifrost warning

add_library(kallisto_core ${sources} ${headers})
target_include_directories(kallisto_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# add_executable(kallisto main.cpp)
add_executable(kallisto main.cpp ${sources} ${headers})
find_package( Threads REQUIRED )
ExternalProject_Get_Property(bifrost install_dir)

target_link_libraries(kallisto kallisto_core pthread ${install_dir}/build/src/libbifrost.a)

if (USE_BAM)
target_link_libraries(kallisto ${CMAKE_CURRENT_SOURCE_DIR}/../ext/htslib/libhts.a)
endif(USE_BAM)

add_dependencies(kallisto bifrost)

if(CUDA_AVAILABLE)
    add_executable(gpukallisto gpumain.cu GPUProcessReads.cu GPUKernels.cu GPUIndex.cu)

    # Set CUDA properties for gpukallisto
    set_target_properties(gpukallisto PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
        CUDA_STANDARD 17
    )

    # Add required CUDA flags for cuCollections and per-thread default streams
    target_compile_options(gpukallisto PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-extended-lambda>
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
        $<$<COMPILE_LANGUAGE:CUDA>:--default-stream per-thread>
    )

    find_package(nvcomp REQUIRED)
    target_link_libraries(gpukallisto kallisto_core pthread
        ${CMAKE_CURRENT_SOURCE_DIR}/../ext/htslib/libhts.a
        cuco CUDA::cudart nvcomp::nvcomp)

    # Link gpukallisto against bifrost
    target_link_libraries(gpukallisto ${install_dir}/build/src/libbifrost.a)

    add_dependencies(gpukallisto bifrost)
endif(CUDA_AVAILABLE)

# ============================================================================
# metalkallisto — Apple Metal port (macOS / Apple Silicon only)
# ============================================================================
if(APPLE)
    # Compile kallisto.metal -> .air -> .metallib at build time
    set(METAL_SRC  ${CMAKE_CURRENT_SOURCE_DIR}/kernels/kallisto.metal)
    set(METAL_AIR  ${CMAKE_CURRENT_BINARY_DIR}/kallisto.air)
    set(METAL_LIB  ${CMAKE_CURRENT_BINARY_DIR}/kallisto.metallib)

    execute_process(
        COMMAND xcrun -f metal
        OUTPUT_VARIABLE METAL_EXECUTABLE
        OUTPUT_STRIP_TRAILING_WHITESPACE
        RESULT_VARIABLE METAL_FIND_RESULT
    )
    if(NOT METAL_FIND_RESULT EQUAL 0 OR METAL_EXECUTABLE STREQUAL "")
        message(FATAL_ERROR "Unable to locate the Metal compiler via xcrun -f metal")
    endif()
    get_filename_component(METAL_BIN_DIR "${METAL_EXECUTABLE}" DIRECTORY)
    set(METALLIB_EXECUTABLE "${METAL_BIN_DIR}/metallib")

    add_custom_command(
        OUTPUT  ${METAL_AIR}
        COMMAND ${METAL_EXECUTABLE} -c ${METAL_SRC} -o ${METAL_AIR}
        DEPENDS ${METAL_SRC}
        COMMENT "Compiling Metal shader to AIR"
    )
    add_custom_command(
        OUTPUT  ${METAL_LIB}
        COMMAND ${METALLIB_EXECUTABLE} ${METAL_AIR} -o ${METAL_LIB}
        DEPENDS ${METAL_AIR}
        COMMENT "Linking Metal AIR to metallib"
    )
    add_custom_target(kallisto_metallib DEPENDS ${METAL_LIB})

    add_executable(metalkallisto
        metalmain.mm
        MetalUtils.mm
        MetalIndex.mm
        MetalPipeline.mm
        MetalEM.mm
        MetalProcessReads.mm
    )

    add_dependencies(metalkallisto kallisto_metallib bifrost)

    set_target_properties(metalkallisto PROPERTIES
        OBJCXX_STANDARD 17
        OBJCXX_STANDARD_REQUIRED ON
    )

    target_compile_definitions(metalkallisto PRIVATE
        METALLIB_PATH="${METAL_LIB}"
    )

    target_include_directories(metalkallisto PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${EXT_PROJECTS_DIR}/bifrost/build/src
        ${EXT_PROJECTS_DIR}/bifrost/src
    )

    # install_dir already set by ExternalProject_Get_Property(bifrost install_dir) above
    target_link_libraries(metalkallisto
        kallisto_core
        pthread
        ${install_dir}/build/src/libbifrost.a
        ${METAL_FRAMEWORK}
        ${MPS_FRAMEWORK}
        ${FOUNDATION_FRAMEWORK}
    )
    # zlib is linked below, after find_package(ZLIB) is called
endif(APPLE)

if(LINK MATCHES static)
    set(BUILD_SHARED_LIBS OFF)
    set(HDF5_USE_STATIC_LIBRARIES 1)

    if (UNIX AND NOT APPLE)
        #set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libstdc++")
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
    endif(UNIX AND NOT APPLE)

    SET_TARGET_PROPERTIES(kallisto  PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINK MATCHES static)


if(USE_HDF5)
    find_package( HDF5 COMPONENTS CXX REQUIRED )
endif(USE_HDF5)

if (NOT ZLIBNG)
find_package( ZLIB REQUIRED )

if ( ZLIB_FOUND )
    if (DEFINED ZLIB_INCLUDE_DIRS)
        include_directories( ${ZLIB_INCLUDE_DIRS} )
    elseif (DEFINED ZLIB_INCLUDE_DIR)
        include_directories( ${ZLIB_INCLUDE_DIR} )
    else()
        message(FATAL_ERROR "zlib found but no include directories are set.")
    endif()
    target_link_libraries(kallisto kallisto_core ${ZLIB_LIBRARIES})
    if(CUDA_AVAILABLE AND TARGET gpukallisto)
        target_link_libraries(gpukallisto ${ZLIB_LIBRARIES})
    endif()
    if(APPLE AND TARGET metalkallisto)
        target_link_libraries(metalkallisto ${ZLIB_LIBRARIES})
    endif()
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )
endif(NOT ZLIBNG)

if (ZLIBNG)
    add_definitions(-DKALLISTO_USE_ZLIB_NG)
    include_directories(../ext/zlib-ng/zlib-ng)
    include_directories(../ext/zlib-ng)
    if(WIN32)
    target_link_libraries(kallisto ${CMAKE_CURRENT_SOURCE_DIR}/../ext/zlib-ng/zlib-ng/libz.lib)
    else()
    target_link_libraries(kallisto ${CMAKE_CURRENT_SOURCE_DIR}/../ext/zlib-ng/zlib-ng/libz.a)
    endif(WIN32)
endif(ZLIBNG)

if(USE_HDF5)
    #if(HDF5_FOUND)
        include_directories( ${HDF5_INCLUDE_DIRS} )
        #target_link_libraries( kallisto_core ${HDF5_LIBRARIES} )
        target_link_libraries( kallisto ${HDF5_LIBRARIES} )
    #else()
    #    message(FATAL_ERROR "HDF5 not found. Required to output files")
    #endif()
endif(USE_HDF5)

if(LINK MATCHES static)
    if (UNIX AND NOT APPLE)
        target_link_libraries(kallisto librt.a)
    endif()
else()
    if (UNIX AND NOT APPLE)
        target_link_libraries(kallisto rt)
    endif()
endif(LINK MATCHES static)


install(TARGETS kallisto DESTINATION "${CMAKE_INSTALL_BINDIR}")
