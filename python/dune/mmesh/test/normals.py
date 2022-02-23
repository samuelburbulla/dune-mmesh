from dune.grid import reader
from dune.mmesh import mmesh
from dune.mmesh import normals
from dune.fem.function import integrate


# 2D
from dune.mmesh.test.grids import line
grid = mmesh((reader.gmsh, line.filename), 2)
igrid = grid.hierarchicalGrid.interfaceGrid

normal = normals(igrid)

igrid.writeVTK("normals-2d", pointdata={"n": normal})

int = integrate(igrid, [abs(normal[0]), abs(normal[1])], order=0)
assert(abs(int[0]) < 1e-8)
assert(abs(int[1] - 1.) < 1e-8)


# 3D
from dune.mmesh.test.grids import plane
grid  = mmesh((reader.gmsh, plane.filename), 3)
igrid = grid.hierarchicalGrid.interfaceGrid

normal = normals(igrid)

igrid.writeVTK("normals-3d", pointdata={"n": normal})

int = integrate(igrid, [abs(normal[0]), abs(normal[1]), abs(normal[2])], order=0)
assert(abs(int[0]) < 1e-8)
assert(abs(int[1]) < 1e-8)
assert(abs(int[2] - 1.) < 1e-8)
