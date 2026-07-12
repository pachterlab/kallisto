add_library(libchunkdecoding INTERFACE)
add_library(rapidgzip::chunkdecoding ALIAS libchunkdecoding)
target_include_directories(libchunkdecoding INTERFACE ${CMAKE_CURRENT_LIST_DIR}/../..)
target_sources(libchunkdecoding INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/Bzip2Chunk.hpp
    ${CMAKE_CURRENT_LIST_DIR}/GzipChunk.hpp
    ${CMAKE_CURRENT_LIST_DIR}/DecompressionError.hpp
)
target_link_libraries(libchunkdecoding INTERFACE core rapidgzip::gzip rapidgzip::libblockfinder indexed_bzip2)
target_compile_features(libchunkdecoding INTERFACE cxx_std_17)
