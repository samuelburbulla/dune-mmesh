set(HEADERS
  cachingentity.hh
  common.hh
  connectedcomponent.hh
  cutsettriangulation.hh
  dgfparser.hh
  declaration.hh
  entity.hh
  entityseed.hh
  explicitgridfactory.hh
  geometry.hh
  gmshparser.hh
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
  multiid.hh
  pointfieldvector.hh
  rangegenerators.hh
  structeredgridfactory.hh
  )

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/grid)
