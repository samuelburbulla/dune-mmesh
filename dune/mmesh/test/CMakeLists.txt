add_subdirectory(examples EXCLUDE_FROM_ALL)

dune_symlink_to_source_files(FILES
  grids
)

dune_add_test(NAME test-grid-2d SOURCES test-grid.cc)
set_property( TARGET test-grid-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
target_link_libraries(test-grid-2d PUBLIC "dunemmesh")

dune_add_test(NAME test-grid-3d SOURCES test-grid.cc)
set_property( TARGET test-grid-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
target_link_libraries(test-grid-3d PUBLIC "dunemmesh")

dune_add_test(NAME test-mmesh-2d SOURCES test-mmesh-2d.cc)
target_link_libraries(test-mmesh-2d PUBLIC "dunemmesh")

dune_add_test(NAME test-mmesh-3d SOURCES test-mmesh-3d.cc)
target_link_libraries(test-mmesh-3d PUBLIC "dunemmesh")

dune_add_test(NAME test-twists-2d SOURCES test-twists.cc)
set_property( TARGET test-twists-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
target_link_libraries(test-twists-2d PUBLIC "dunemmesh")

dune_add_test(NAME test-twists-3d SOURCES test-twists.cc)
set_property( TARGET test-twists-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
target_link_libraries(test-twists-3d PUBLIC "dunemmesh")

dune_add_test(NAME test-mmesh-structured-2d SOURCES test-mmesh-structured.cc)
set_property( TARGET test-mmesh-structured-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
target_link_libraries(test-mmesh-structured-2d PUBLIC "dunemmesh")

dune_add_test(NAME test-mmesh-structured-3d SOURCES test-mmesh-structured.cc)
set_property( TARGET test-mmesh-structured-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
target_link_libraries(test-mmesh-structured-3d PUBLIC "dunemmesh")

dune_add_test(NAME test-interfaceiterator-2d SOURCES test-interfaceiterator.cc)
set_property( TARGET test-interfaceiterator-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
target_link_libraries(test-interfaceiterator-2d PUBLIC "dunemmesh")

dune_add_test(NAME test-interfaceiterator-3d SOURCES test-interfaceiterator.cc)
set_property( TARGET test-interfaceiterator-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
target_link_libraries(test-interfaceiterator-3d PUBLIC "dunemmesh")

dune_add_test(NAME test-interfacegrid-2d SOURCES test-interfacegrid-2d.cc)
set_property( TARGET test-interfacegrid-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
target_link_libraries(test-interfacegrid-2d PUBLIC "dunemmesh")

dune_add_test(NAME test-interfacegrid-3d SOURCES test-interfacegrid-3d.cc)
set_property( TARGET test-interfacegrid-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
target_link_libraries(test-interfacegrid-3d PUBLIC "dunemmesh")

if(DEFINED MPFVERSION)
  dune_add_test(NAME test-intersectionvolume-2d SOURCES test-intersectionvolume.cc)
  set_property( TARGET test-intersectionvolume-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2;CGAL_INTERSECTION" )
  target_link_libraries(test-intersectionvolume-2d PUBLIC "dunemmesh")
endif()

dune_add_test(NAME test-curvature-2d SOURCES test-curvature.cc)
set_property( TARGET test-curvature-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
target_link_libraries(test-curvature-2d PUBLIC "dunemmesh")

dune_add_test(NAME test-curvature-3d SOURCES test-curvature.cc)
set_property( TARGET test-curvature-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
target_link_libraries(test-curvature-3d PUBLIC "dunemmesh")

dune_add_test(NAME test-femadaptation-dg SOURCES test-femadaptation.cc)
set_property(TARGET test-femadaptation-dg APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
target_link_libraries(test-femadaptation-dg PUBLIC "dunemmesh")

#dune_add_test(NAME test-threading SOURCES test-threading.cc)
#set_property(TARGET test-threading APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
#target_link_libraries(test-threading PUBLIC "dunemmesh")

#dune_add_test(NAME test-threading-interface SOURCES test-threading.cc)
#set_property(TARGET test-threading-interface APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2;INTERFACE" )
#target_link_libraries(test-threading-interface PUBLIC "dunemmesh")

#dune_add_test(NAME test-threading-3d SOURCES test-threading.cc)
#set_property(TARGET test-threading-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )
#target_link_libraries(test-threading-3d PUBLIC "dunemmesh")

#dune_add_test(NAME test-threading-interface-3d SOURCES test-threading.cc)
#set_property(TARGET test-threading-interface-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3;INTERFACE" )
#target_link_libraries(test-threading-interface-3d PUBLIC "dunemmesh")
