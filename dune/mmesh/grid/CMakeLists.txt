set(HEADERS
  cachingentity.hh
  common.hh
  connectedcomponent.hh
  dgfparser.hh
  declaration.hh
  entity.hh
  entityseed.hh
  explicitgridfactory.hh
  geometry.hh
  gmshparser.hh
  gridfactory.hh
  hierarchiciterator.hh
  implicitgridfactory.hh
  incidentiterator.hh
  indexsets.hh
  interfaceiterator.hh
  intersectioniterator.hh
  intersections.hh
  leafiterator.hh
  mmeshdefaults.hh
  multiid.hh
  pointfieldvector.hh
  rangegenerators.hh
  structeredgridfactory.hh
  )

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/grid)
