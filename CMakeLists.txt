cmake_minimum_required(VERSION 3.0)

project(kallisto)

# add_compile_options(-Wall)
add_compile_options(-std=c++11)
add_compile_options(-Wdeprecated-register)

add_subdirectory(src)

add_subdirectory(test)
enable_testing()
add_test(ExampleTest test/exTests)
