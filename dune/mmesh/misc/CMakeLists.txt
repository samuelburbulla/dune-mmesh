set(HEADERS
  boundaryidprovider.hh
  capabilities.hh
  partitionhelper.hh
  persistentcontainer.hh
  twistutility.hh
)

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/misc)
