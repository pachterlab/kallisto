file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)

if (USE_BAM)
include_directories(../ext/htslib)
endif(USE_BAM)

add_compile_options(-Wno-subobject-linkage) # Suppress bifrost warning

#add_library(kallisto_core ${sources} ${headers})
#target_include_directories( PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB_RECURSE sources *.cpp)
file(GLOB_RECURSE headers *.h *.hpp)
list(REMOVE_ITEM sources ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_executable(kallisto ${sources} ${headers} main.cpp)

if (DEFINED ENV{RELEASE_OS} AND "$ENV{RELEASE_OS}" STREQUAL "windows")
    set_target_properties(kallisto PROPERTIES LINK_FLAGS "-Wl,--allow-multiple-definition")
endif()
#add_executable(kallisto main.cpp)

find_package( Threads REQUIRED )
ExternalProject_Get_Property(bifrost install_dir)
if (USE_BAM)
target_link_libraries(kallisto  pthread ${CMAKE_CURRENT_SOURCE_DIR}/../ext/htslib/libhts.a ${install_dir}/build/src/libbifrost.a)
else()
target_link_libraries(kallisto  pthread ${install_dir}/build/src/libbifrost.a)
endif(USE_BAM)

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
    find_package( HDF5 REQUIRED )
endif(USE_HDF5)

if (NOT ZLIBNG)
find_package( ZLIB REQUIRED )

if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries(kallisto  ${ZLIB_LIBRARIES})
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
    if(HDF5_FOUND)
        include_directories( ${HDF5_INCLUDE_DIRS} )
        #target_link_libraries(  ${HDF5_LIBRARIES} )
        target_link_libraries( kallisto ${HDF5_LIBRARIES} )
    else()
        message(FATAL_ERROR "HDF5 not found. Required to output files")
    endif()
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
