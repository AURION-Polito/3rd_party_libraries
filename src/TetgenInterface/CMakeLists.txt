list(APPEND 3rd_party_libraries_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/TetgenInterface.cpp)

list(APPEND 3rd_party_libraries_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/TetgenInterface.hpp)

list(APPEND 3rd_party_libraries_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(3rd_party_libraries_src_sources ${3rd_party_libraries_src_sources} PARENT_SCOPE)
set(3rd_party_libraries_src_headers ${3rd_party_libraries_src_headers} PARENT_SCOPE)
set(3rd_party_libraries_src_includes ${3rd_party_libraries_src_includes} PARENT_SCOPE)
