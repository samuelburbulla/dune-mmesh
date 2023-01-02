"""Test distance function."""

from ufl import grad
from dune.grid import reader
from dune.mmesh import mmesh
from dune.mmesh import distance
from dune.fem.function import integrate
from dune.mmesh.test.grids import line
from dune.mmesh.test.grids import plane


# 2D
grid = mmesh((reader.gmsh, line.filename), 2)
igrid = grid.hierarchicalGrid.interfaceGrid

dist = distance(grid)

normal = grad(dist)
grid.writeVTK("distance-2d", pointdata={"d": dist, "n": [normal[0], normal[1], 0]}, nonconforming=True)
igrid.writeVTK("distance-2d-interface")

integral = integrate(grid, dist, order=5)
nor = integrate(grid, normal, order=5)

assert abs(integral - 0.25) < 1e-8
assert (nor * nor) < 1e-8


# 3D
grid  = mmesh((reader.gmsh, plane.filename), 3)
igrid = grid.hierarchicalGrid.interfaceGrid

dist = distance(grid)

normal = grad(dist)
grid.writeVTK("distance-3d", pointdata={"d": dist, "n": normal}, nonconforming=True)
igrid.writeVTK("distance-3d-interface")

integral = integrate(grid, dist, order=5)
nor = integrate(grid, normal, order=5)

assert abs(integral - 0.25) < 1e-3
assert (nor * nor) < 1e-8
