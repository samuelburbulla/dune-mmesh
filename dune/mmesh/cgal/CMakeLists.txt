set(HEADERS
  defaults.hh
  includes.hh
  triangulationwrapper.hh
)

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/cgal)
