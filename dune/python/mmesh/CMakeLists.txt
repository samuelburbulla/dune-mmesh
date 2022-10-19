set(HEADERS
  distance.hh
  grid.hh
  jacobian_iterative.hh
  jacobian.hh
  pyskeletontrace.hh
  utility.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/python/mmesh)
