set(HEADERS
  cachingentity.hh
  common.hh
  connectedcomponent.hh
  curvatureoperator.hh
  dgfparser.hh
  entity.hh
  entityseed.hh
  geometry.hh
  grid.hh
  gridfactory.hh
  hierarchiciterator.hh
  incidentiterator.hh
  indexsets.hh
  intersectioniterator.hh
  intersections.hh
  leafiterator.hh
  traits.hh
  )

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/interface)
