add_subdirectory(dune)

configure_file(setup.py.in setup.py)

dune_python_install_package(
  PATH "."
  CMAKE_METADATA_FILE dune/mmesh/metadata.cmake
  CMAKE_METADATA_FLAGS DUNE_OPTS_FILE
)
