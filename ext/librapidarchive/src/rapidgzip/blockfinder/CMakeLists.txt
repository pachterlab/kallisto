add_library(libblockfinder INTERFACE)
add_library(rapidgzip::libblockfinder ALIAS libblockfinder)
target_include_directories(libblockfinder INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../..)
target_sources(libblockfinder INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/Bgzf.hpp
    ${CMAKE_CURRENT_LIST_DIR}/DynamicHuffman.hpp
    ${CMAKE_CURRENT_LIST_DIR}/Interface.hpp
    ${CMAKE_CURRENT_LIST_DIR}/PigzStringView.hpp
    ${CMAKE_CURRENT_LIST_DIR}/precodecheck/CountAllocatedLeaves.hpp
)
target_link_libraries(libblockfinder INTERFACE core)
target_compile_features(libblockfinder INTERFACE cxx_std_17)
