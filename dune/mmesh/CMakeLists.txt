add_subdirectory(grid)
add_subdirectory(test)

set(HEADERS
  mmesh.hh
)

install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh)
