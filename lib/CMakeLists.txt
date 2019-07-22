# build library libdunemmesh
dune_add_library(dunemmesh
  _DUNE_TARGET_OBJECTS:mmesh_
  ${PROJECT_SOURCE_DIR}/dune/mmesh/mmesh.cc)

target_compile_options(dunemmesh PRIVATE -include${CMAKE_BINARY_DIR}/config.h)

target_link_libraries(dunemmesh ${DUNE_LIBS})
target_link_libraries(dunemmesh CGAL)
target_link_libraries(dunemmesh mpfr)
