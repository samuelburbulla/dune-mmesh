from dune.grid import reader
from dune.mmesh import mmesh
from dune.mmesh import interfaceIndicator
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.function import integrate
from ufl import *


# 2D
file = "line.msh"
grid = mmesh((reader.gmsh, file), 2)
igrid = grid.hierarchicalGrid.interfaceGrid

I = interfaceIndicator(igrid)

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

grid.writeVTK("interfaceindicator-2d", pointdata=[uh])

intU = integrate(grid, uh, order=1)
assert(abs(intU - 0.5) < 1e-3)


# 3D
file = "plane.msh"
grid  = mmesh((reader.gmsh, file), 3)
igrid = grid.hierarchicalGrid.interfaceGrid

I = interfaceIndicator(igrid)

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

grid.writeVTK("interfaceindicator-2d", pointdata=[uh])

intU = integrate(grid, uh, order=1)
assert(abs(intU - 0.5) < 1e-3)
