file(GLOB sources *.cpp)
file(GLOB headers *.h)

list(REMOVE_ITEM sources main.cpp)

add_library(kallisto_core ${sources} ${headers})
target_include_directories(kallisto_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(kallisto main.cpp)

target_link_libraries(kallisto kallisto_core)
