from dune.grid import reader, cartesianDomain
from dune.mmesh import mmesh
from dune.alugrid import aluSimplexGrid
from dune.fem.space import dglagrange
from dune.fem.function import integrate
from ufl import *
from dune.ufl import Constant, DirichletBC
from dune.fem.scheme import galerkin
from time import time
from dune.fem import parameter
parameter.append({"fem.verboserank": 0})

domain = cartesianDomain([0, 0], [1, 1], [10, 10])
gridView = mmesh(domain)
#gridView = aluSimplexGrid(domain)

space = dglagrange(gridView, order=1)
x = SpatialCoordinate(space)
uh = space.interpolate(0, name="uh")

print("uh.size =", uh.size)

u = TrialFunction(space)
v = TestFunction(space)

exact = sin(pi*x[0]) * sin(pi*x[1])
beta = Constant(10, name="beta")
n = FacetNormal(space)
h = avg(FacetArea(space))
hBnd = FacetArea(space)

A = inner(grad(u), grad(v)) * dx
A += beta / h * inner(jump(u), jump(v)) * dS
A -= dot(dot(avg(grad(u)), n('+')), jump(v)) * dS
A -= dot(dot(avg(grad(v)), n('+')), jump(u)) * dS
A += beta / hBnd * inner(u - exact, v) * ds
A -= dot(dot(grad(u), n), v) * ds
A -= dot(dot(grad(v), n), u - exact) * ds

b = -div( grad(exact) ) * v * dx

scheme = galerkin([A == b],
  solver='cg',
  parameters={
    "newton.linear.maxiterations": 1000,
    "newton.linear.verbose": "true"
  })

start = time()
scheme.solve(uh)
print(f"Took {time()-start:.6f}")

gridView.writeVTK("laplace", pointdata={"uh": uh, "exact": exact})

l2 = integrate(gridView, sqrt(dot(uh-exact, uh-exact)), order=5)
print(l2)
assert(l2 < 1e-3)
