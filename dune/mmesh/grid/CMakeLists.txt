set(HEADERS
  common.hh
  dgfparser.hh
  entity.hh
  entityseed.hh
  explicitgridfactory.hh
  geometry.hh
  gmshparser.hh
  gridfactory.hh
  hierarchiciterator.hh
  implicitgridfactory.hh
  incidentiterator.hh
  interfaceiterator.hh
  indexsets.hh
  intersectioniterator.hh
  intersections.hh
  leafiterator.hh
  mmeshdefaults.hh
  pointfieldvector.hh
  rangegenerators.hh
  structeredgridfactory.hh
  )

install(FILES ${HEADERS}
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/grid)
