add_subdirectory(huffman)
add_subdirectory(gzip)
add_subdirectory(blockfinder)
add_subdirectory(chunkdecoding)


add_library(librapidgzip INTERFACE)
add_library(rapidgzip::rapidgzip ALIAS librapidgzip)
target_include_directories(librapidgzip INTERFACE ${CMAKE_CURRENT_LIST_DIR}/..)
target_sources(librapidgzip INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/ChunkData.hpp
    ${CMAKE_CURRENT_LIST_DIR}/CompressedVector.hpp
    ${CMAKE_CURRENT_LIST_DIR}/DecodedData.hpp
    ${CMAKE_CURRENT_LIST_DIR}/DecodedDataView.hpp
    ${CMAKE_CURRENT_LIST_DIR}/GzipBlockFinder.hpp
    ${CMAKE_CURRENT_LIST_DIR}/GzipChunkFetcher.hpp
    ${CMAKE_CURRENT_LIST_DIR}/IndexFileFormat.hpp
    ${CMAKE_CURRENT_LIST_DIR}/MarkerReplacement.hpp
    ${CMAKE_CURRENT_LIST_DIR}/ParallelGzipReader.hpp
    ${CMAKE_CURRENT_LIST_DIR}/rapidgzip.hpp
    ${CMAKE_CURRENT_LIST_DIR}/WindowMap.hpp
)
target_link_libraries(librapidgzip INTERFACE
    core
    rapidgzip::chunkdecoding
)
target_compile_features(librapidgzip INTERFACE cxx_std_17)
if(LIBRAPIDARCHIVE_WITH_ISAL)
    target_compile_definitions(librapidgzip INTERFACE LIBRAPIDARCHIVE_WITH_ISAL)
endif()
