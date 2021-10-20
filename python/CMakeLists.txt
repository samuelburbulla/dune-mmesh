add_subdirectory(dune)

configure_file(setup.py.in setup.py)

dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/mmesh/metadata.cmake
  DEPENDS _grid2d _grid3d _interfacegrid2d _interfacegrid3d
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
if ("${CMAKE_VERSION}" VERSION_GREATER 3.14)
  cmake_policy(SET CMP0087 NEW)
endif ("${CMAKE_VERSION}" VERSION_GREATER 3.14)
