## @example mmesh.py
#  This is an example of how to use MMesh with dune-python

from dune.grid import reader
from dune.mmesh import mmesh

dim = 2
file = "../grids/horizontal2d.msh"

# MMesh
print("= Read " + str(dim) + "d grid from file '" + file + "'")
gridView = mmesh((reader.gmsh, file), dim)

assert gridView.size(0) == 270
print("Found", gridView.size(0), "entities")

print("Write grid to 'test-python-mmesh'")
gridView.writeVTK("test-python-mmesh")

# InterfaceGrid
igridView = gridView.hierarchicalGrid.interfaceGrid

assert igridView.size(0) == 5
print("Found", igridView.size(0), "interface entities")

print("Write interface grid to 'test-python-mmesh-interface'")
igridView.writeVTK("test-python-mmesh-interface")


# Interaction
hgrid = gridView.hierarchicalGrid

print("\n= Obtain interface elements by intersections")
for e in gridView.elements:
    for i in gridView.intersections(e):
        if hgrid.isInterface(i):
            ie = hgrid.asInterfaceEntity(i)
            print("Interface entity centered at:", ie.geometry.center)

print("\n= Obtain intersections by interface elements")

for ie in igridView.elements:
    inters = hgrid.asIntersection(ie)
    print("Intersection with normal:", inters.centerUnitOuterNormal)
