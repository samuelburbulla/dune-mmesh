from dune.grid import reader
from dune.mmesh import mmesh
from dune.fem import adapt
from dune.fem.view import adaptiveLeafGridView as adaptive
from dune.fem.space import dglagrange
from dune.fem.function import integrate

from ufl import SpatialCoordinate, conditional

from dune.mmesh.test.grids import line
grid = mmesh((reader.gmsh, line.filename), 2)
hgrid = grid.hierarchicalGrid
gridView = adaptive(grid)
igridView = adaptive(hgrid.interfaceGrid)

space = dglagrange(gridView, order=1)
x = SpatialCoordinate(space)
uh = space.interpolate(x[0]+x[1], name="uh")

ispace = dglagrange(igridView, order=1)
iuh = ispace.interpolate(x[0], name="iuh")

def addSomeInterface(number):
  for e in gridView.elements:
    for i in gridView.intersections(e):
      if hgrid.isInterface(i) or i.boundary:
        continue
      if number == 0:
        print("Add interface at", i.geometry.center)
        hgrid.addInterface(i)
        return
      number -= 1

addSomeInterface(0)
addSomeInterface(10)
addSomeInterface(20)
addSomeInterface(42)
adapt([uh])
adapt([iuh])

# Initialize new elements with 0
iuh.interpolate(conditional(abs(iuh - x[0]) < 1e-6, iuh, 0))

gridView.writeVTK("adaptation", pointdata=[uh], nonconforming=True)
igridView.writeVTK("adaptation-interface", pointdata=[iuh], nonconforming=True)

intuh = integrate(gridView, uh, order=1)
intiuh = integrate(igridView, iuh, order=1)
print(intuh, intiuh)
assert(abs(intuh - 1) < 1e-6)
assert(abs(intiuh - 0.5) < 1e-6)
