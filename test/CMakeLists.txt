project(Tests)

include_directories(../Catch/include)
include_directories(../src)

add_executable(tests tests.cpp)

target_link_libraries(tests kallisto_core)
