set(HEADERS
  boundaryidprovider.hh
  capabilities.hh
  communication.hh
  objectstream.hh
  partitionhelper.hh
  persistentcontainer.hh
  twistutility.hh
)

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/misc)
