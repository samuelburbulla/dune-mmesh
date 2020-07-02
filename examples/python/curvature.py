import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from dune.grid import reader
from dune.mmesh import mmesh

dim = 2
file = "../grids/junction.msh"

gridView = mmesh((reader.gmsh, file), dim)
hgrid = gridView.hierarchicalGrid

gridView.writeVTK("domain")

from ufl import *
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.ufl import DirichletBC

igridView = hgrid.interfaceGrid
space = lagrange(igridView, dimRange=dim, order=1)
u = TrialFunction(space)
phi = TestFunction(space)

x = SpatialCoordinate(space)
a = inner(u, phi) * dx - inner(grad(x), grad(phi)) * dx
bc = DirichletBC(space, [0]*dim)

solution = space.interpolate(x, name="solution")
scheme = galerkin([a == 0, bc], space)
scheme.solve(target=solution)

vectorspace = lagrange(igridView, dimRange=3, order=1)
curv = vectorspace.interpolate( as_vector([solution[0], solution[1], 0 if dim == 2 else solution[2]]), name="curvatureTimesNormal" )
igridView.writeVTK("curvature", pointdata=[curv])

from dune.fem.function import integrate
print("error =", integrate(igridView, sqrt(inner(curv, curv))-1, space.order+1))
