# INCLUDE SUB DIRECTORIES
add_subdirectory(Common)
add_subdirectory(VoroInterface)
add_subdirectory(TriangleInterface)
add_subdirectory(TetgenInterface)
add_subdirectory(VtkInterface)
add_subdirectory(MetisInterface)
add_subdirectory(PETScInterface)
add_subdirectory(SuiteSparseInterface)
add_subdirectory(MKLInterface)

list(APPEND 3rd_party_libraries_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(3rd_party_libraries_src_sources ${3rd_party_libraries_src_sources} PARENT_SCOPE)
set(3rd_party_libraries_src_headers ${3rd_party_libraries_src_headers} PARENT_SCOPE)
set(3rd_party_libraries_src_includes ${3rd_party_libraries_src_includes} PARENT_SCOPE)
