project(Tests)

include_directories(../Catch/include)
include_directories(../src)

file(GLOB sources *.cpp)

add_executable(tests ${sources})

target_link_libraries(tests kallisto_core)
