## @example mmesh.py
#  This is an example of how to use MMesh with dune-python

from dune.grid import reader
from dune.mmesh import mmesh

dim = 2
file = "../grids/horizontal2d.msh"

# MMesh
gridView = mmesh((reader.gmsh, file), dim)

assert gridView.size(0) == 270

gridView.writeVTK("test-python-mmesh")

# InterfaceGrid
igridView = gridView.hierarchicalGrid.interfaceGrid

assert igridView.size(0) == 5

igridView.writeVTK("test-python-mmesh-interface")
