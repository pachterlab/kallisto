project(Tests)

include_directories(../Catch/include)
include_directories(../src)

file(GLOB sources *.cpp)

add_executable(tests ${sources})

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( tests ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

target_link_libraries(tests kallisto_core)
