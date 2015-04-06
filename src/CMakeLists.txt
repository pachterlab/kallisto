file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp h5kal.cpp)

include_directories(../ext/)
add_compile_options(-DSEQAN_HAS_ZLIB)
add_library(kallisto_core ${sources} ${headers})
target_include_directories(kallisto_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(kallisto main.cpp)

# add_executable(h5kal h5kal.cpp)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( kallisto ${ZLIB_LIBRARIES} )
else()
    message(FATA_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )

if(HDF5_FOUND)
    include_directories( ${HDF5_INCLUDE_DIR} )
    # target_link_libraries( h5kal ${HDF5_LIBRARIES} )
    target_link_libraries( kallisto ${HDF5_LIBRARIES} )
else()
    message(FATAL_ERROR "HDF5 not found. Required to output files")
endif()

target_link_libraries(kallisto kallisto_core pthread)

if (UNIX AND NOT APPLE)
    target_link_libraries(kallisto rt)
endif()

install(TARGETS kallisto DESTINATION bin)
