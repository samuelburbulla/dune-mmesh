add_subdirectory(examples EXCLUDE_FROM_ALL)

dune_symlink_to_source_files(FILES
  grids
  parameter
)

install(FILES massoperator.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/mmesh/test)

dune_add_test(NAME test-grid-2d SOURCES test-grid.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property( TARGET test-grid-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-grid-3d SOURCES test-grid.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property( TARGET test-grid-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

dune_add_test(NAME test-mmesh-2d SOURCES test-mmesh-2d.cc)

dune_add_test(NAME test-mmesh-3d SOURCES test-mmesh-3d.cc)

dune_add_test(NAME test-twists-2d SOURCES test-twists.cc)
set_property( TARGET test-twists-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-twists-3d SOURCES test-twists.cc)
set_property( TARGET test-twists-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

dune_add_test(NAME test-mmesh-structured-2d SOURCES test-mmesh-structured.cc)
set_property( TARGET test-mmesh-structured-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-mmesh-structured-3d SOURCES test-mmesh-structured.cc)
set_property( TARGET test-mmesh-structured-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

dune_add_test(NAME test-interfaceiterator-2d SOURCES test-interfaceiterator.cc)
set_property( TARGET test-interfaceiterator-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-interfaceiterator-3d SOURCES test-interfaceiterator.cc)
set_property( TARGET test-interfaceiterator-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

dune_add_test(NAME test-interfacegrid-2d SOURCES test-interfacegrid-2d.cc)
set_property( TARGET test-interfacegrid-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-interfacegrid-3d SOURCES test-interfacegrid-3d.cc)
set_property( TARGET test-interfacegrid-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

if(DEFINED MPFVERSION)
  dune_add_test(NAME test-intersectionvolume-2d SOURCES test-intersectionvolume.cc)
  set_property( TARGET test-intersectionvolume-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2;CGAL_INTERSECTION" )
endif()

dune_add_test(NAME test-curvature-2d SOURCES test-curvature.cc)
set_property( TARGET test-curvature-2d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-curvature-3d SOURCES test-curvature.cc)
set_property( TARGET test-curvature-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

# Workaround to fix linking issue in dune-fem
if(OPENMP_FOUND)
  if(APPLE)
    get_filename_component(OPENMP_LIB_PATH ${OpenMP_CXX_LIBRARIES} DIRECTORY)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp -L${OPENMP_LIB_PATH}")
  endif()
endif()
# End of workaround

dune_add_test(NAME test-femadaptation-dg SOURCES test-femadaptation.cc)
set_property(TARGET test-femadaptation-dg APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-femadaptation SOURCES test-femadaptation.cc)
set_property(TARGET test-femadaptation APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2; LAGRANGE" )

#dune_add_test(NAME test-threading SOURCES test-threading.cc)
#set_property(TARGET test-threading APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

#dune_add_test(NAME test-threading-interface SOURCES test-threading.cc)
#set_property(TARGET test-threading-interface APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2;INTERFACE" )

#dune_add_test(NAME test-threading-3d SOURCES test-threading.cc)
#set_property(TARGET test-threading-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

#dune_add_test(NAME test-threading-interface-3d SOURCES test-threading.cc)
#set_property(TARGET test-threading-interface-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3;INTERFACE" )

dune_add_test(NAME test-adapt-tuplespace SOURCES test-adapt-tuplespace.cc)
set_property(TARGET test-adapt-tuplespace APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-distance SOURCES test-distance.cc)
set_property(TARGET test-distance APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-l2projection SOURCES test-l2projection.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property(TARGET test-l2projection APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-l2projection-interface SOURCES test-l2projection.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property(TARGET test-l2projection-interface APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2;INTERFACE" )

dune_add_test(NAME test-comm SOURCES test-comm.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property(TARGET test-comm APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )

dune_add_test(NAME test-comm-3d SOURCES test-comm.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property(TARGET test-comm-3d APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3" )

dune_add_test(NAME test-comm-interface SOURCES test-comm.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property(TARGET test-comm-interface APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2;INTERFACE" )

dune_add_test(NAME test-comm-3d-interface SOURCES test-comm.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property(TARGET test-comm-3d-interface APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=3;INTERFACE" )

dune_add_test(NAME test-mpi SOURCES test-mpi.cc MPI_RANKS 1 2 4 8 TIMEOUT 300)
set_property(TARGET test-mpi APPEND PROPERTY COMPILE_DEFINITIONS "GRIDDIM=2" )
