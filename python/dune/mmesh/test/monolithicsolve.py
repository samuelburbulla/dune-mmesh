from time import time
from mpi4py import MPI

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
beta  = Constant(3.01, name="beta")
n = FacetNormal(space)
ni = FacetNormal(ispace)
I = interfaceIndicator(igridView)

a  = inner(grad(u), grad(v)) * dx
a += beta * inner(jump(u), jump(v)) * (1-I)*dS
a -= dot(dot(avg(grad(u)), n('+')), jump(v)) * (1-I)*dS
a -= dot(dot(avg(grad(v)), n('+')), jump(u)) * (1-I)*dS
a += beta * inner(u - 0, v) * ds
a -= dot(dot(grad(u), n), v) * ds
a -= dot(dot(grad(v), n), u - 0) * ds

ia  = inner(grad(iu), grad(iv)) * dx
ia += beta * inner(jump(iu), jump(iv)) * dS
ia -= dot(dot(avg(grad(iu)), ni('+')), jump(iv)) * dS
ia -= dot(dot(avg(grad(iv)), ni('+')), jump(iu)) * dS
ib  = q * iv * dx

from dune.mmesh import skeleton, trace
omega = Constant(1e-6, name="omega")
a -= (skeleton(iuh)('+') - u('+')) / omega * v('+') * I*dS
a -= (skeleton(iuh)('-') - u('-')) / omega * v('-') * I*dS
ia += (iu - trace(uh)('+')) / omega * iv * dx
ia += (iu - trace(uh)('-')) / omega * iv * dx

from dune.fem.scheme import galerkin
scheme  = galerkin([a == 0], solver='cg')
ischeme = galerkin([ia == ib], solver='cg')
from dune.mmesh import monolithicSolve
dt = -time()
monolithicSolve(schemes=(scheme, ischeme), targets=(uh, iuh), verbose=True)
dt += time()

rank = MPI.COMM_WORLD.Get_rank()
if rank == 0:
  print(f"Took {dt:.6f}")

from dune.fem.function import integrate
intBulk = integrate(gridView, uh, order=1)
intInterface = integrate(igridView, iuh, order=1)

assert(abs(intBulk - 0.06473048252528972) < 1e-6)
assert(abs(intInterface - 0.19198977421898397) < 1e-6)
