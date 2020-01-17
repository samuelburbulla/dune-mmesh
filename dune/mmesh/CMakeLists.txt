add_subdirectory(cgal)
add_subdirectory(interface)
add_subdirectory(grid)
add_subdirectory(remeshing)
add_subdirectory(test)

set(HEADERS
  edgelengthindicator.hh
  mmesh.hh)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh)
