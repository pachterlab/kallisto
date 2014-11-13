project(Tests)

include_directories(../Catch/include)
include_directories(../src)

add_executable(exTests example_tests.cpp)

target_link_libraries(exTests xwcore)
