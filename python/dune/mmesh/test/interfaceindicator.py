from dune.grid import reader
from dune.mmesh import mmesh
from dune.mmesh import interfaceIndicator
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.function import integrate
from dune.fem.view import geometryGridView
from ufl import *

# 2D
from dune.mmesh.test.grids import line
grid2 = mmesh((reader.gmsh, line.filename), 2)

# 3D
from dune.mmesh.test.grids import plane
grid3  = mmesh((reader.gmsh, plane.filename), 3)

# 2D GeometryGrid
space = lagrange(grid2, order=1, dimRange=2)
x = SpatialCoordinate(space)
disp = space.interpolate(x, name="disp")
geoGrid = geometryGridView(disp)


for grid in [grid2, grid3, geoGrid]:
  igrid = grid.hierarchicalGrid.interfaceGrid
  I = interfaceIndicator(igrid, grid=grid)

  space = lagrange(grid, order=1)
  u  = TrialFunction(space)
  uu = TestFunction(space)
  n = FacetNormal(space)
  x = SpatialCoordinate(space)
  topbottom = conditional(x[1] < 1e-6, 1, 0) + conditional(x[1] > 1-1e-6, 1, 0)

  A  = inner(grad(u), grad(uu)) * dx
  A += 1e6 * (avg(u) - 1) * avg(uu) * I*dS
  A += 1e6 * (u - 0) * uu * topbottom * ds

  scheme = galerkin([A == 0])
  uh = space.function(name="uh")
  scheme.solve(uh)

  grid.writeVTK("interfaceindicator", pointdata=[uh])

  intU = integrate(grid, uh, order=1)
  print(intU)
  assert(abs(intU - 0.5) < 1e-2)
