cmake_minimum_required(VERSION 3.13)
project(dune-mmesh CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

find_package(dune-grid REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-grid_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)
include(GridType)

# start a dune project with information from dune.module
dune_project()

dune_enable_all_packages()

# make sure we have CGAL
find_package(CGAL REQUIRED)

add_subdirectory(examples EXCLUDE_FROM_ALL)
add_subdirectory(dune)
add_subdirectory(doc)
add_subdirectory(cmake/modules)

# if Python bindings are enabled, include necessary sub directories.
if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory("python")
  dune_python_install_package(PATH "python")
endif()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
