file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)

include_directories(../ext/)

add_library(kallisto_core ${sources} ${headers})
target_include_directories(kallisto_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(kallisto main.cpp)

if(LINK MATCHES static)
    set(BUILD_SHARED_LIBS OFF)
    set(HDF5_USE_STATIC_LIBRARIES 1)

    if (UNIX AND NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
    endif(UNIX AND NOT APPLE)

endif(LINK MATCHES static)

find_package( HDF5 REQUIRED )
find_package( ZLIB REQUIRED )

if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )

if(HDF5_FOUND)
    include_directories( ${HDF5_INCLUDE_DIR} )
    target_link_libraries( kallisto_core ${HDF5_LIBRARIES} )
    target_link_libraries( kallisto ${HDF5_LIBRARIES} )
else()
    message(FATAL_ERROR "HDF5 not found. Required to output files")
endif()

target_link_libraries(kallisto kallisto_core pthread)

if (UNIX AND NOT APPLE)
    target_link_libraries(kallisto rt)
endif()

install(TARGETS kallisto DESTINATION bin)
