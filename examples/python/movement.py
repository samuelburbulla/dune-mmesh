## @example mmesh.py
#  This is an example of how to use MMesh with dune-python

import numpy as np
import dune.geometry
from dune.grid import reader, gridFunction
from dune.mmesh import mmesh

dim = 2
file = "../grids/horizontal2d.msh"

# Create MMesh and InterfaceGrid
grid = mmesh((reader.gmsh, file), dim)
igrid = grid.hierarchicalGrid.interfaceGrid

# Define the movement of the interface
def movement(x):
  m = np.array(x)
  m -= np.full(dim, 0.5)
  m *= 2e-2 * 3.1415
  return np.array([ m[1], -m[0] ])

def getShifts():
  mapper = igrid.mapper({dune.geometry.vertex: 1})
  shifts = np.zeros((mapper.size, dim))
  for v in igrid.vertices:
    shifts[ mapper.index(v) ] = movement( v.geometry.center )
  return shifts

hgrid = grid.hierarchicalGrid

vtk = grid.sequencedVTK("movement")
ivtk = igrid.sequencedVTK("movement-interface")

# Start the time loop
for t in range(0, 101):
  print("t =", t)
  hgrid.preAdapt()
  hgrid.ensureInterfaceMovement( getShifts() )
  hgrid.markElements()
  hgrid.adapt()
  # transfer data some data...
  hgrid.moveInterface( getShifts() )
  hgrid.postAdapt()
  vtk()
  ivtk()
