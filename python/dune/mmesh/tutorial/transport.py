from dune.grid import reader, cartesianDomain
from dune.mmesh import mmesh
from dune.ufl import Constant
import ufl
from dune.fem.space import *
from dune.fem.scheme import galerkin

dim = 2
domain = cartesianDomain([0]*dim, [1]*dim, [10]*dim)
gridView = mmesh(domain)

space = dglagrange(gridView, order=1)
u = ufl.TrialFunction(space)
v = ufl.TestFunction(space)

x = ufl.SpatialCoordinate(space.cell())
n = ufl.FacetNormal(space.cell())

def speed():
    return ufl.as_vector([1.0]*dim)

def f(u):
    return speed() * u

def g(u):
    return ufl.inner( ufl.conditional( ufl.inner(speed(), n('+')) > 0, f(u('+')), f(u('-')) ), n('+') )

tau = Constant(0.01, name="timeStep")
uh_old = space.interpolate(0, name="uh_old")

a = (u - uh_old) / tau * v * ufl.dx
a -= ufl.inner( f(u), ufl.grad(v) ) * ufl.dx
a += g(u) * ufl.jump(v) * ufl.dS
a += ufl.inner(f(u), n) * v * ufl.ds

scheme = galerkin([a==0], solver=("suitesparse","umfpack"))

uh = space.interpolate(x[0], name="uh")

gridView.writeVTK("transport-"+str(0), pointdata={"u": uh}, nonconforming=True)
for step in range(1, 11):
  uh_old.assign(uh)
  scheme.solve(target=uh)
  gridView.writeVTK("transport-"+str(step), pointdata={"u": uh}, nonconforming=True)
