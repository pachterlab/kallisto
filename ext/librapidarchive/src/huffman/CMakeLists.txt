add_library(libhuffman INTERFACE)
target_include_directories(libhuffman INTERFACE ${CMAKE_CURRENT_LIST_DIR}/..)
target_sources(libhuffman INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/HuffmanCodingBase.hpp
    ${CMAKE_CURRENT_LIST_DIR}/HuffmanCodingLinearSearch.hpp
    ${CMAKE_CURRENT_LIST_DIR}/HuffmanCodingShortBitsCached.hpp
    ${CMAKE_CURRENT_LIST_DIR}/HuffmanCodingSymbolsPerLength.hpp
)
target_link_libraries(libhuffman INTERFACE core)
target_compile_features(libhuffman INTERFACE cxx_std_17)
