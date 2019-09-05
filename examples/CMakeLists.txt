dune_symlink_to_source_files(FILES
  grids
)

add_executable("cgal-example" cgal-example.cc)
target_link_dune_default_libraries("cgal-example")
target_link_libraries("cgal-example" dunemmesh)

add_executable("dgf-example" dgf-example.cc)
target_link_dune_default_libraries("dgf-example")
target_link_libraries("dgf-example" dunemmesh)

add_executable("dgf-interval-example" dgf-interval-example.cc)
target_link_dune_default_libraries("dgf-interval-example")
target_link_libraries("dgf-interval-example" dunemmesh)

add_executable("gmsh-example-2d" gmsh-example.cc)
target_compile_definitions("gmsh-example-2d" PRIVATE GRIDDIM=2)
target_link_dune_default_libraries("gmsh-example-2d")
target_link_libraries("gmsh-example-2d" dunemmesh)

add_executable("gmsh-example-3d" gmsh-example.cc)
target_compile_definitions("gmsh-example-3d" PRIVATE GRIDDIM=3)
target_link_dune_default_libraries("gmsh-example-3d")
target_link_libraries("gmsh-example-3d" dunemmesh)

add_executable("gmsh-implicit-example-2d" gmsh-implicit-example.cc)
target_compile_definitions("gmsh-implicit-example-2d" PRIVATE GRIDDIM=2)
target_link_dune_default_libraries("gmsh-implicit-example-2d")
target_link_libraries("gmsh-implicit-example-2d" dunemmesh)

add_executable("gmsh-implicit-example-3d" gmsh-implicit-example.cc)
target_compile_definitions("gmsh-implicit-example-3d" PRIVATE GRIDDIM=3)
target_link_dune_default_libraries("gmsh-implicit-example-3d")
target_link_libraries("gmsh-implicit-example-3d" dunemmesh)

add_executable("interfacegrid-example-2d" interfacegrid-example.cc)
target_compile_definitions("interfacegrid-example-2d" PRIVATE GRIDDIM=2)
target_link_dune_default_libraries("interfacegrid-example-2d")
target_link_libraries("interfacegrid-example-2d" dunemmesh)

add_executable("interfacegrid-example-3d" interfacegrid-example.cc)
target_compile_definitions("interfacegrid-example-3d" PRIVATE GRIDDIM=3)
target_link_dune_default_libraries("interfacegrid-example-3d")
target_link_libraries("interfacegrid-example-3d" dunemmesh)
