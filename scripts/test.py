# Grid creation
from dune.grid import reader
from dune.mmesh import mmesh
dim = 2
file = "horizontal.msh"
gridView  = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid



# Bulk problem
from ufl import *
from dune.ufl import DirichletBC
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.function import integrate
import time

space = lagrange(gridView, order=3)
u = TrialFunction(space)
v = TestFunction(space)

x = SpatialCoordinate(space)
exact = sin(x[0]*x[1]*4*pi)
f = -div(grad(exact))

a = inner(grad(u), grad(v)) * dx
b = f * v * dx

start_time = time.time()
scheme = galerkin([a == b, DirichletBC(space, exact)], solver=("suitesparse", "umfpack"))
print("bulk:", time.time() - start_time)
uh = space.interpolate(0, name="solution")
scheme.solve(target=uh)

def L2(u1, u2):
    return sqrt(integrate(u1.grid, (u1-u2)**2, order=5))

assert( L2(uh, exact) < 1e-6 )



# Interface problem
ispace = lagrange(igridView, order=3)
iuh = ispace.interpolate(0, name="isolution")

iu = TrialFunction(ispace)
iv = TestFunction(ispace)

ix = SpatialCoordinate(ispace)
iexact = sin(0.5*ix[dim-2]*4*pi)
iF = -div(grad(iexact))

ia = inner(grad(iu), grad(iv)) * dx
ib = iF * iv * dx

start_time = time.time()
ischeme = galerkin([ia == ib, DirichletBC(ispace, iexact)])
print("interface:", time.time() - start_time)
ischeme.solve(target=iuh)
assert( L2(iuh, iexact) < 1e-6 )



# Couple bulk to surface
from dune.mmesh import trace
tr = avg(trace(uh))
ib = inner(grad(tr), grad(iv)) * dx

iuh.interpolate(0)
ischeme = galerkin([ia == ib, DirichletBC(ispace, avg(trace(uh)))])
ischeme.solve(target=iuh)
assert( L2(iuh, iexact) < 1e-6 )



# Couple surface to bulk
from dune.mmesh import skeleton
sk = skeleton(iuh)
b = avg(sk) * avg(v) * dS

uh.interpolate(0)
scheme = galerkin([a == b, DirichletBC(space, 0)])
scheme.solve(target=uh)



# Coupled solve
from dune.mmesh import monolithicSolve, interfaceIndicator
from dune.fem.space import dglagrange
from dune.ufl import Constant

from dune.fem.view import geometryGridView
from dune.fem.function import uflFunction
x = SpatialCoordinate(triangle)

vectorSpace = dglagrange(gridView, dimRange=2, order=1)
position = vectorSpace.interpolate(x, name="position")
geoGrid = geometryGridView(position)

space = dglagrange(geoGrid, order=3)
u = TrialFunction(space)
v = TestFunction(space)

ispace = lagrange(igridView, order=3)
iu = TrialFunction(ispace)
iv = TestFunction(ispace)

uh = space.interpolate(0, name="solution")
iuh = ispace.interpolate(0, name="isolution")

I = interfaceIndicator(igridView, grid=geoGrid)
n = FacetNormal(space)
n_g = FacetNormal(ispace)
beta = Constant(1e2, name="beta")
omega = Constant(1e-6, name="omega")

a  = inner(grad(u), grad(v)) * dx
a += beta * inner(jump(u), jump(v)) * (1-I)*dS
a -= dot(dot(avg(grad(u)), n('+')), jump(v)) * (1-I)*dS

a += beta * inner(u - 0, v) * ds
a -= dot(dot(grad(u), n), v) * ds

sk = skeleton(iuh, grid=geoGrid)
a -= (sk('+') - u('+')) / omega * v('+') * I*dS
a -= (sk('-') - u('-')) / omega * v('-') * I*dS


ia  = inner(grad(iu), grad(iv)) * dx
ia += beta * inner(jump(iu), jump(iv)) * dS
ia -= inner(avg(grad(iu)), n_g('+')) * jump(iv) * dS

tr  = trace(uh)
ia += (iu - tr('+')) / omega * iv * dx
ia += (iu - tr('-')) / omega * iv * dx

ib  = 1 * iv * dx

scheme = galerkin([a == 0])
ischeme = galerkin([ia == ib])

start_time = time.time()
monolithicSolve(schemes=(scheme, ischeme), targets=(uh, iuh), verbose=True)
print("monolithicSolve:", time.time() - start_time)
assert( (max(iuh.as_numpy) - 0.14063) < 1e-3 )
