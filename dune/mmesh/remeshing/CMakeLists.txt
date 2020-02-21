set(HEADERS
  distance.hh
  longestedgerefinement.hh
  ratioindicator.hh
)

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/remeshing)
