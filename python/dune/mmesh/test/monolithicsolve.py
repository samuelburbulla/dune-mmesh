from time import time

from dune.mmesh.test.grids import tjunction

from dune.grid import reader
from dune.mmesh import mmesh
gridView  = mmesh((reader.gmsh, tjunction.name), 2)
igridView = gridView.hierarchicalGrid.interfaceGrid

from ufl import *
from dune.fem.space import lagrange, dglagrange
space = dglagrange(gridView, order=1)
ispace = dglagrange(igridView, order=1)
u = TrialFunction(space)
v = TestFunction(space)
iu = TrialFunction(ispace)
iv = TestFunction(ispace)
uh = space.interpolate(0, name="uh")
iuh = ispace.interpolate(0, name="iuh")

from dune.ufl import Constant
from dune.mmesh import interfaceIndicator
q     = Constant(1, name="q")
omega = Constant(1e-6, name="omega")
beta  = Constant(1, name="beta")
n = FacetNormal(space)
ni = FacetNormal(ispace)
I = interfaceIndicator(igridView)

a  = inner(grad(u), grad(v)) * dx
a += beta * inner(jump(u), jump(v)) * (1-I)*dS
a -= dot(dot(avg(grad(u)), n('+')), jump(v)) * (1-I)*dS
a += beta * inner(u - 0, v) * ds
a -= dot(dot(grad(u), n), v) * ds

ia  = inner(grad(iu), grad(iv)) * dx
ia += beta * inner(jump(iu), jump(iv)) * dS
ia -= dot(dot(avg(grad(iu)), ni('+')), jump(iv)) * dS
ib  = q * iv * dx

from dune.mmesh import skeleton, trace
omega = Constant(1e-6, name="omega")
a -= (skeleton(iuh)('+') - u('+')) / omega * v('+') * I*dS
a -= (skeleton(iuh)('-') - u('-')) / omega * v('-') * I*dS
ia += (iu - trace(uh)('+')) / omega * iv * dx
ia += (iu - trace(uh)('-')) / omega * iv * dx

from dune.fem.scheme import galerkin
scheme  = galerkin([a == 0])
ischeme = galerkin([ia == ib])
from dune.mmesh import monolithicSolve
dt = -time()
monolithicSolve(schemes=(scheme, ischeme), targets=(uh, iuh), verbose=True)
dt += time()
print(f"Took {dt:.6f}")

from dune.fem.function import integrate
intBulk = integrate(gridView, uh, order=1)
intInterface = integrate(igridView, iuh, order=1)

print(intBulk, intInterface)
assert(abs(intBulk - 0.066) < 1e-3)
assert(abs(intInterface - 0.192) < 1e-3)
