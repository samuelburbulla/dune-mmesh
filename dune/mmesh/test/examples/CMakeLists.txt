dune_symlink_to_source_files(FILES
  grids
)

dune_add_test(NAME "cgal" SOURCES cgal.cc)

dune_add_test(NAME "dgf" SOURCES dgf.cc)

dune_add_test(NAME "dgf-interval" SOURCES dgf-interval.cc)

dune_add_test(NAME "gmsh-2d" SOURCES gmsh.cc)
target_compile_definitions("gmsh-2d" PRIVATE GRIDDIM=2)

dune_add_test(NAME "gmsh-3d" SOURCES gmsh.cc)
target_compile_definitions("gmsh-3d" PRIVATE GRIDDIM=3)

dune_add_test(NAME "gmsh-implicit-2d" SOURCES gmsh-implicit.cc)
target_compile_definitions("gmsh-implicit-2d" PRIVATE GRIDDIM=2)

dune_add_test(NAME "gmsh-implicit-3d" SOURCES gmsh-implicit.cc)
target_compile_definitions("gmsh-implicit-3d" PRIVATE GRIDDIM=3)

dune_add_test(NAME "interfacegrid-2d" SOURCES interfacegrid.cc)
target_compile_definitions("interfacegrid-2d" PRIVATE GRIDDIM=2)

dune_add_test(NAME "interfacegrid-3d" SOURCES interfacegrid.cc)
target_compile_definitions("interfacegrid-3d" PRIVATE GRIDDIM=3)

dune_add_test(NAME "movement-2d" SOURCES movement.cc)
target_compile_definitions("movement-2d" PRIVATE GRIDDIM=2)

dune_add_test(NAME "interaction-2d" SOURCES interaction.cc)
target_compile_definitions("interaction-2d" PRIVATE GRIDDIM=2)

dune_add_test(NAME "interaction-3d" SOURCES interaction.cc)
target_compile_definitions("interaction-3d" PRIVATE GRIDDIM=3)
