set(HEADERS
  connectedcomponent.hh
  entity.hh
  entityseed.hh
  geometry.hh
  grid.hh
  gridfactory.hh
  hierarchiciterator.hh
  indexsets.hh
  incidentiterator.hh
  intersectioniterator.hh
  intersections.hh
  leafiterator.hh
  traits.hh
  )

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/interface)
