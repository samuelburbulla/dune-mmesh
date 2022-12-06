from dune.grid import reader
from dune.mmesh import mmesh

from dune.mmesh.test.grids import line
grid = mmesh((reader.gmsh, line.filename), 2)
igrid = grid.hierarchicalGrid.interfaceGrid

mpi_size = 1
try:
  from mpi4py import MPI
  mpi_size = MPI.COMM_WORLD.Get_size()
except:
  pass

size = [[2684], [1375], [], [703, 735]]
isize = [[100], [51], [], [27, 26]]

assert(grid.size(0) in size[mpi_size-1])
assert(igrid.size(0) in isize[mpi_size-1])
