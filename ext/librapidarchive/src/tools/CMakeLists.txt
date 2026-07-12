if(LIBRAPIDARCHIVE_BUILD_TOOLS)
add_executable(ibzip2)
target_sources(ibzip2 PRIVATE ${CMAKE_CURRENT_LIST_DIR}/ibzip2.cpp)
target_link_libraries(ibzip2 PRIVATE cxxopts indexed_bzip2_parallel)

add_executable(bitpatternstats)
target_sources(bitpatternstats PRIVATE ${CMAKE_CURRENT_LIST_DIR}/bitpatternstats.cpp)
target_link_libraries(bitpatternstats PRIVATE cxxopts indexed_bzip2_parallel)

add_executable(generateLUTs)
target_sources(generateLUTs PRIVATE ${CMAKE_CURRENT_LIST_DIR}/generateLUTs.cpp)
target_link_libraries(generateLUTs PRIVATE cxxopts librapidgzip)

add_executable(rapidgzip)
target_sources(rapidgzip PRIVATE ${CMAKE_CURRENT_LIST_DIR}/rapidgzip.cpp)
target_link_libraries(rapidgzip PRIVATE cxxopts librapidgzip)
if(LIBRAPIDARCHIVE_WITH_RPMALLOC)
    target_link_libraries(rapidgzip PRIVATE rpmalloc)
    target_compile_definitions(rapidgzip PRIVATE LIBRAPIDARCHIVE_WITH_RPMALLOC)
endif()

if(LIBRAPIDARCHIVE_RUN_CLANG_TIDY)
    set_target_properties(ibzip2 PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
    set_target_properties(rapidgzip PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_COMMAND}")
endif()

endif()
