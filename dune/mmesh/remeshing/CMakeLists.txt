set(HEADERS
  cellcenterrefinement.hh
  edgelengthratioindicator.hh
  interfaceindicator.hh
  interfacerrefinement.hh
  newestvertexbisection.hh
  longestedgerefinement.hh
  ratioindicator.hh
  volumeindicator.hh
)

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/remeshing)
