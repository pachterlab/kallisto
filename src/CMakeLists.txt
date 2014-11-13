file(GLOB sources *.cpp)
file(GLOB headers *.h)

list(REMOVE_ITEM sources main.cpp)

add_library(xwcore ${sources} ${headers})
target_link_libraries(xwcore gff)
target_include_directories(xwcore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${gff_SOURCE_DIR})

add_executable(xprs_wiggle main.cpp)

target_link_libraries(xprs_wiggle xwcore)
