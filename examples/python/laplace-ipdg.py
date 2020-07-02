## @example laplace-ipdg.py
#  Laplace with Interior-Penalty Discontinuous Galerkin method

import io
from dune.grid import reader
from dune.mmesh import mmesh, trace, skeleton
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

dim = 2
file = "../grids/horizontal2d.msh"

# MMesh
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid


from ufl import *
import dune.ufl
from dune.fem import parameter
from dune.fem.function import integrate
from dune.fem.space import dglagrange
from dune.fem.scheme import galerkin
parameter.append({"fem.verboserank": 0})
solverParameters =\
   {"newton.tolerance": 1e-8,
    "newton.linear.tolerance": 1e-12,
    "newton.linear.preconditioning.method": "jacobi",
    "newton.linear.maxiterations": 10000,
    "preconditioning.method": "jacobi",
    "newton.verbose": True,
    "newton.linear.verbose": False}


space = dglagrange(gridView, order=3)
x = SpatialCoordinate(space)
u = TrialFunction(space)
v = TestFunction(space)
n = FacetNormal(space.cell())

# create an interface indicator
ispace = dglagrange(igridView, order=0)
interface = ispace.interpolate(1.0, name="interface")
I = skeleton(interface)('+')

exact = -x[1]
mu = 1e3

# bulk problem without interface edges
a = inner(grad(u),grad(v)) * dx
a += mu * jump(u) * jump(v) * (1-I)*dS
a -= inner(avg(grad(u)), n('+')) * jump(v) * (1-I)*dS
a -= inner(avg(grad(v)), n('+')) * jump(u) * (1-I)*dS
a += mu * u * v * ds
a -= u * inner(grad(v), n) * ds
a -= v * inner(grad(u), n) * ds

# source and boundary terms
b = -div(grad(exact)) * v * dx
b += mu * exact * v * ds
b -= exact * inner(grad(v), n) * ds

# additional jump condition at the interface
b += inner(as_vector([1.0]*dim), n('+')) * jump(v) * I*dS

uh = space.interpolate(0,name="solution")
scheme = galerkin([a==b], parameters=solverParameters)
scheme.solve(target=uh)
gridView.writeVTK("laplace-ipdg",
    pointdata={"solution":uh,}, nonconforming=True, subsampling=2)
