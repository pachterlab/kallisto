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

find_package( Threads REQUIRED )
ExternalProject_Get_Property(bifrost install_dir)

target_link_libraries(kallisto kallisto_core pthread ${CMAKE_CURRENT_SOURCE_DIR}/../ext/htslib/libhts.a)
find_package(nvcomp REQUIRED)
target_link_libraries(gpukallisto kallisto_core pthread ${CMAKE_CURRENT_SOURCE_DIR}/../ext/htslib/libhts.a cuco CUDA::cudart nvcomp::nvcomp)

if (USE_BAM)
target_link_libraries(kallisto kallisto_core pthread ${CMAKE_CURRENT_SOURCE_DIR}/../ext/htslib/libhts.a ${install_dir}/build/src/libbifrost.a  -DMAX_KMER_SIZE=${MAX_KMER_SIZE})
else()
target_link_libraries(kallisto kallisto_core pthread ${install_dir}/build/src/libbifrost.a -DMAX_KMER_SIZE=${MAX_KMER_SIZE})
endif(USE_BAM)

# Link gpukallisto against bifrost (provides roaring bitmap, Kmer, UnitigMap, etc.)
target_link_libraries(gpukallisto ${install_dir}/build/src/libbifrost.a)

add_dependencies(kallisto bifrost)
add_dependencies(gpukallisto bifrost)

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
    target_link_libraries(gpukallisto ${ZLIB_LIBRARIES})
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
