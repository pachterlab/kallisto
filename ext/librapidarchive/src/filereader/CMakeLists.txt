add_library(filereader INTERFACE)
target_include_directories(filereader INTERFACE ${CMAKE_CURRENT_LIST_DIR}/..)
target_sources(filereader INTERFACE
    ${CMAKE_CURRENT_LIST_DIR}/BitReader.hpp
    ${CMAKE_CURRENT_LIST_DIR}/Buffered.hpp
    ${CMAKE_CURRENT_LIST_DIR}/BufferView.hpp
    ${CMAKE_CURRENT_LIST_DIR}/FileReader.hpp
    ${CMAKE_CURRENT_LIST_DIR}/Python.hpp
    ${CMAKE_CURRENT_LIST_DIR}/Shared.hpp
    ${CMAKE_CURRENT_LIST_DIR}/SinglePass.hpp
    ${CMAKE_CURRENT_LIST_DIR}/Standard.hpp
    ${CMAKE_CURRENT_LIST_DIR}/StreamAdapter.hpp
)
target_compile_features(filereader INTERFACE cxx_std_17)
