from dune.grid import reader
from dune.mmesh import mmesh
from dune.mmesh import distance
from ufl import grad
from dune.fem.function import integrate


# 2D
file = "line2d.msh"
grid = mmesh((reader.gmsh, file), 2)
dist = distance(grid)

normal = grad(dist)
grid.writeVTK("distance-2d", pointdata={"d": dist, "n": [normal[0], normal[1], 0]}, nonconforming=True)

int = integrate(grid, dist, order=5)
nor = integrate(grid, normal, order=5)

assert(abs(int - 0.25) < 1e-8)
assert((nor * nor) < 1e-8)


# 3D
file = "flat3d.msh"
grid  = mmesh((reader.gmsh, file), 3)
dist = distance(grid)

normal = grad(dist)
grid.writeVTK("distance-3d", pointdata={"d": dist, "n": normal}, nonconforming=True)

from dune.fem.function import integrate
int = integrate(grid, dist, order=5)
nor = integrate(grid, normal, order=5)

assert(abs(int - 0.25) < 1e-3)
assert((nor * nor) < 1e-8)
