# Grid creation
from dune.grid import reader
from dune.mmesh import mmesh
dim = 2
file = "horizontal.msh"
gridView  = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

from dune.mmesh import interfaceIndicator
I = interfaceIndicator(igridView)

from dune.mmesh import distance
dist = distance(gridView)

from ufl import grad
normal = grad(dist)
gridView.writeVTK("distance", pointdata={"d": dist, "n": [normal[0], normal[1], 0]}, nonconforming=True)

from dune.fem.function import integrate
int = integrate(gridView, dist, order=5)
jac = integrate(gridView, grad(dist), order=5)

print(int, jac)
assert(abs(int - 0.273323) < 1e-6)
assert(abs(jac[0]) < 1e-6 and abs(jac[1]) < 1e-6)
