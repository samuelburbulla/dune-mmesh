dune_symlink_to_source_files(FILES
  grids
)

dune_add_test(NAME "cgal" SOURCES cgal.cc)
target_link_libraries("cgal" PUBLIC "dunemmesh")

dune_add_test(NAME "dgf" SOURCES dgf.cc)
target_link_libraries("dgf" PUBLIC "dunemmesh")

dune_add_test(NAME "dgf-interval" SOURCES dgf-interval.cc)
target_link_libraries("dgf-interval" PUBLIC "dunemmesh")

dune_add_test(NAME "gmsh-2d" SOURCES gmsh.cc)
target_compile_definitions("gmsh-2d" PRIVATE GRIDDIM=2)
target_link_libraries("gmsh-2d" PUBLIC "dunemmesh")

dune_add_test(NAME "gmsh-3d" SOURCES gmsh.cc)
target_compile_definitions("gmsh-3d" PRIVATE GRIDDIM=3)
target_link_libraries("gmsh-3d" PUBLIC "dunemmesh")

dune_add_test(NAME "gmsh-implicit-2d" SOURCES gmsh-implicit.cc)
target_compile_definitions("gmsh-implicit-2d" PRIVATE GRIDDIM=2)
target_link_libraries("gmsh-implicit-2d" PUBLIC "dunemmesh")

dune_add_test(NAME "gmsh-implicit-3d" SOURCES gmsh-implicit.cc)
target_compile_definitions("gmsh-implicit-3d" PRIVATE GRIDDIM=3)
target_link_libraries("gmsh-implicit-3d" PUBLIC "dunemmesh")

dune_add_test(NAME "interfacegrid-2d" SOURCES interfacegrid.cc)
target_compile_definitions("interfacegrid-2d" PRIVATE GRIDDIM=2)
target_link_libraries("interfacegrid-2d" PUBLIC "dunemmesh")

dune_add_test(NAME "interfacegrid-3d" SOURCES interfacegrid.cc)
target_compile_definitions("interfacegrid-3d" PRIVATE GRIDDIM=3)
target_link_libraries("interfacegrid-3d" PUBLIC "dunemmesh")

dune_add_test(NAME "movement-2d" SOURCES movement.cc)
target_compile_definitions("movement-2d" PRIVATE GRIDDIM=2)
target_link_libraries("movement-2d" PUBLIC "dunemmesh")

dune_add_test(NAME "interaction-2d" SOURCES interaction.cc)
target_compile_definitions("interaction-2d" PRIVATE GRIDDIM=2)
target_link_libraries("interaction-2d" PUBLIC "dunemmesh")

dune_add_test(NAME "interaction-3d" SOURCES interaction.cc)
target_compile_definitions("interaction-3d" PRIVATE GRIDDIM=3)
target_link_libraries("interaction-3d" PUBLIC "dunemmesh")
