add_library(librapidgzipgzip INTERFACE)
add_library(rapidgzip::gzip ALIAS librapidgzipgzip)
target_include_directories(librapidgzipgzip INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../..)
target_sources(librapidgzipgzip INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/crc32.hpp
    ${CMAKE_CURRENT_LIST_DIR}/definitions.hpp
    ${CMAKE_CURRENT_LIST_DIR}/deflate.hpp
    ${CMAKE_CURRENT_LIST_DIR}/format.hpp
    ${CMAKE_CURRENT_LIST_DIR}/GzipAnalyzer.hpp
    ${CMAKE_CURRENT_LIST_DIR}/gzip.hpp
    ${CMAKE_CURRENT_LIST_DIR}/GzipReader.hpp
    ${CMAKE_CURRENT_LIST_DIR}/InflateWrapper.hpp
    ${CMAKE_CURRENT_LIST_DIR}/isal.hpp
    ${CMAKE_CURRENT_LIST_DIR}/precode.hpp
    ${CMAKE_CURRENT_LIST_DIR}/RFCTables.hpp
    ${CMAKE_CURRENT_LIST_DIR}/zlib.hpp
)
target_link_libraries(librapidgzipgzip INTERFACE core rapidgzip::gzip::huffman zlibstatic)
if(LIBRAPIDARCHIVE_WITH_ISAL)
    target_link_libraries(librapidgzipgzip INTERFACE
        core
        rapidgzip::gzip::huffman
        rapidgzip::libblockfinder
        ISAL::isal
        zlibstatic)
endif()
target_compile_features(librapidgzipgzip INTERFACE cxx_std_17)
