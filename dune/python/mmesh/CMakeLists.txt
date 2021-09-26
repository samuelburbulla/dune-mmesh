set(HEADERS
  interfacegrid.hh
  jacobian.hh
  pyskeletontrace.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/python/mmesh)
