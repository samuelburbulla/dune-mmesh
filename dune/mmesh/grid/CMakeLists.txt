set(HEADERS
  cachingentity.hh
  common.hh
  connectedcomponent.hh
  cutsettriangulation.hh
  declaration.hh
  dgfparser.hh
  entity.hh
  entityseed.hh
  explicitgridfactory.hh
  geometry.hh
  gmshgridfactory.hh
  gmshreader.hh
  gridfactory.hh
  hierarchiciterator.hh
  implicitgridfactory.hh
  incidentiterator.hh
  indexsets.hh
  interfaceiterator.hh
  intersectioniterator.hh
  intersections.hh
  leafiterator.hh
  mmesh.hh
  multiid.hh
  pointfieldvector.hh
  polygoncutting.hh
  rangegenerators.hh
  structuredgridfactory.hh
  )

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/grid)
