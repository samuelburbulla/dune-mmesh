set(HEADERS
  boundaryprovider.hh
  pyskeletonfunction.hh
  twistutility.hh
)

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/misc)
