file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)

add_library(kallisto_core ${sources} ${headers})
target_include_directories(kallisto_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(kallisto main.cpp)


find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries( kallisto ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

target_link_libraries(kallisto kallisto_core)

install(TARGETS kallisto DESTINATION bin)
