add_subdirectory(dune)

dune_python_install_package(
  PATH "."
  DEPENDS _grid2d _grid3d _interfacegrid2d _interfacegrid3d
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
if(POLICY CMP0087)
  cmake_policy(SET CMP0087 NEW)
endif()
