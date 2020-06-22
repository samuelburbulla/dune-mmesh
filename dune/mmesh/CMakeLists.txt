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
