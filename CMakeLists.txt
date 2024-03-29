cmake_minimum_required(VERSION 3.13)
project(dune-mmesh CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()
dune_enable_all_packages()

find_package(Boost)
install(DIRECTORY "CGAL"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory(python)
endif()

add_subdirectory(cmake/modules)
add_subdirectory(doc)
add_subdirectory(dune)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
