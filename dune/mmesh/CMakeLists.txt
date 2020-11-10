add_subdirectory(cgal)
add_subdirectory(grid)
add_subdirectory(interface)
add_subdirectory(misc)
add_subdirectory(remeshing)
add_subdirectory(test)

set(HEADERS
  mmesh.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh)

dune_add_library("dunemmesh"
  mmesh.cc)

#if(${CMAKE_VERSION} VERSION_GREATER_EQUAL "3.16")
#    target_precompile_headers(dunemmesh
#      PUBLIC
#        cgal/includes.hh
#    )
#endif()
