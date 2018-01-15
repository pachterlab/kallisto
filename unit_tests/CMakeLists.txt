project(Tests)

# include_directories(../Catch/include)
# include_directories(../src)

file(GLOB sources *.cpp)

add_executable(tests ${sources})
add_test(unittest tests)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( tests ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

target_link_libraries(tests kallisto_core)

find_package( HDF5 REQUIRED )
if(HDF5_FOUND)
    include_directories( ${HDF5_INCLUDE_DIR} )
    target_link_libraries( tests ${HDF5_LIBRARIES} )
else()
    message(FATAL_ERROR "HDF5 not found. Required to output files")
endif()
