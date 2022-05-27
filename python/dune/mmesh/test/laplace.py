from dune.grid import reader, cartesianDomain
from dune.mmesh import mmesh
from dune.alugrid import aluSimplexGrid, aluCubeGrid
from dune.fem.space import lagrange
from dune.fem.function import integrate
from ufl import *
from dune.ufl import DirichletBC
from dune.fem.scheme import galerkin
from time import time

domain = cartesianDomain([0,0],[1,1],[16,16])
gridView = mmesh(domain)
#from dune.mmesh.test.grids import line
#gridView = mmesh((reader.gmsh, line.filename), 2)
#igridView = hgrid.interfaceGrid

#gridView = aluSimplexGrid((reader.gmsh, line.filename), 2)\
hgrid = gridView.hierarchicalGrid


space = lagrange(gridView, order=1)
x = SpatialCoordinate(space)
uh = space.interpolate(0, name="uh")

print("UH:", uh.size)

#ispace = lagrange(igridView, order=1)
#iuh = ispace.interpolate(x[0], name="iuh")

u = TrialFunction(space)
v = TestFunction(space)

exact = sin(pi*x[0]) * sin(pi*x[1])

A = inner(grad(u), grad(v)) * dx
b = -div( grad(exact) ) * v * dx
bc = DirichletBC(space, exact)

scheme = galerkin([A == b, bc], solver='cg',
  parameters={
    "newton.linear.maxiterations": 1000,
    "newton.linear.verbose": "true"
  })

start = time()
scheme.solve(uh)
print(f"Took {time()-start:.2f}")

gridView.writeVTK("laplace", pointdata={"uh": uh, "exact": exact})
#igridView.writeVTK("adaptation-interface", pointdata=[iuh], nonconforming=True)

l2 = integrate(gridView, sqrt(dot(uh-exact, uh-exact)), order=5)
print(l2)
#intiuh = integrate(igridView, iuh, order=1)
#print(intuh, intiuh)
#assert(l1 < 1e-3)
#assert(abs(intiuh - 0.5) < 1e-6)
