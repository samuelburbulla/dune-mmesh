dune_symlink_to_source_files(FILES
  grids
)

dune_add_test(NAME test-grid-2d SOURCES test-grid.cc)
set_property( TARGET test-grid-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-grid-3d SOURCES test-grid.cc)
set_property( TARGET test-grid-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

dune_add_test(NAME test-mmesh-2d SOURCES test-mmesh-2d.cc)

dune_add_test(NAME test-mmesh-3d SOURCES test-mmesh-3d.cc)

dune_add_test(NAME test-twists-2d SOURCES test-twists.cc)
set_property( TARGET test-twists-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-twists-3d SOURCES test-twists.cc)
set_property( TARGET test-twists-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

# dune_add_test(NAME test-carousel-2d SOURCES test-carousel.cc)
# set_property( TARGET test-carousel-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-mmesh-structured-2d SOURCES test-mmesh-structured.cc)
set_property( TARGET test-mmesh-structured-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-mmesh-structured-3d SOURCES test-mmesh-structured.cc)
set_property( TARGET test-mmesh-structured-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

dune_add_test(NAME test-interfaceiterator-2d SOURCES test-interfaceiterator.cc)
set_property( TARGET test-interfaceiterator-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-interfaceiterator-3d SOURCES test-interfaceiterator.cc)
set_property( TARGET test-interfaceiterator-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

# dune_add_test(NAME test-colliding-2d SOURCES test-colliding-2d.cc)
# set_property( TARGET test-colliding-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-interfacegrid-2d SOURCES test-interfacegrid-2d.cc)
set_property( TARGET test-interfacegrid-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-interfacegrid-3d SOURCES test-interfacegrid-3d.cc)
set_property( TARGET test-interfacegrid-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

# dune_add_test(NAME test-fvoninterfacegrid-2d SOURCES test-fvoninterfacegrid.cc)
# set_property( TARGET test-fvoninterfacegrid-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

if(DEFINED MPFVERSION)
  dune_add_test(NAME test-intersectionvolume-2d SOURCES test-intersectionvolume.cc)
  set_property( TARGET test-intersectionvolume-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
 endif()

dune_add_test(NAME test-curvature-2d SOURCES test-curvature.cc)
set_property( TARGET test-curvature-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-curvature-3d SOURCES test-curvature.cc)
set_property( TARGET test-curvature-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

# dune_add_test(NAME test-femadaptation-dg SOURCES test-femadaptation.cc)
# set_property(TARGET test-femadaptation-dg APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
