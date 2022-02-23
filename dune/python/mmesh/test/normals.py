from dune.grid import reader
from dune.mmesh import mmesh
from dune.mmesh import normals
from dune.fem.function import integrate


# 2D
file = "line2d.msh"
grid = mmesh((reader.gmsh, file), 2)
igrid = grid.hierarchicalGrid.interfaceGrid

normal = normals(igrid)

igrid.writeVTK("normals-2d", pointdata={"n": normal})

int = integrate(igrid, normal, order=0)
assert(abs(int[0]) < 1e-8)
assert(abs(int[1] - 1.) < 1e-8)


# 3D
file = "flat3d.msh"
grid  = mmesh((reader.gmsh, file), 3)
igrid = grid.hierarchicalGrid.interfaceGrid

normal = normals(igrid)

igrid.writeVTK("normals-3d", pointdata={"n": normal})

int = integrate(igrid, normal, order=0)
assert(abs(int[0]) < 1e-8)
assert(abs(int[1] - (-1.)) < 1e-8)
assert(abs(int[2]) < 1e-8)
