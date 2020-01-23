set(HEADERS
  defaults.hh
  triangulationwrapper.hh
)

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/cgal)
