add_subdirectory(grids)
add_subdirectory(python)

add_executable("cgal" cgal.cc)
target_link_dune_default_libraries("cgal")

add_executable("dgf" dgf.cc)
target_link_dune_default_libraries("dgf")

add_executable("dgf-interval" dgf-interval.cc)
target_link_dune_default_libraries("dgf-interval")

add_executable("gmsh-2d" gmsh.cc)
target_compile_definitions("gmsh-2d" PRIVATE GRIDDIM=2)
target_link_dune_default_libraries("gmsh-2d")

add_executable("gmsh-3d" gmsh.cc)
target_compile_definitions("gmsh-3d" PRIVATE GRIDDIM=3)
target_link_dune_default_libraries("gmsh-3d")

add_executable("gmsh-implicit-2d" gmsh-implicit.cc)
target_compile_definitions("gmsh-implicit-2d" PRIVATE GRIDDIM=2)
target_link_dune_default_libraries("gmsh-implicit-2d")

add_executable("gmsh-implicit-3d" gmsh-implicit.cc)
target_compile_definitions("gmsh-implicit-3d" PRIVATE GRIDDIM=3)
target_link_dune_default_libraries("gmsh-implicit-3d")

add_executable("interfacegrid-2d" interfacegrid.cc)
target_compile_definitions("interfacegrid-2d" PRIVATE GRIDDIM=2)
target_link_dune_default_libraries("interfacegrid-2d")

add_executable("interfacegrid-3d" interfacegrid.cc)
target_compile_definitions("interfacegrid-3d" PRIVATE GRIDDIM=3)
target_link_dune_default_libraries("interfacegrid-3d")

add_executable("movement-2d" movement.cc)
target_compile_definitions("movement-2d" PRIVATE GRIDDIM=2)
target_link_dune_default_libraries("movement-2d")

add_executable("interaction-2d" interaction.cc)
target_compile_definitions("interaction-2d" PRIVATE GRIDDIM=2)
target_link_dune_default_libraries("interaction-2d")

add_executable("interaction-3d" interaction.cc)
target_compile_definitions("interaction-3d" PRIVATE GRIDDIM=3)
target_link_dune_default_libraries("interaction-3d")
