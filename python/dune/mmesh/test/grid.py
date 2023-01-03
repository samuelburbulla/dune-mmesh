"""Test parallel grid."""

# Determine MPI size and rank
try:
  from mpi4py import MPI
  mpi_size = MPI.COMM_WORLD.Get_size()
  mpi_rank = MPI.COMM_WORLD.Get_rank()
except ImportError:
  mpi_size = 1
  mpi_rank = 0

from dune.grid import reader
from dune.mmesh import mmesh

# Generate mesh file on rank 0
if mpi_rank == 0:
  from dune.mmesh.test.grids import line
  filename = line.filename
else:
  filename = "line.msh"

if mpi_size > 0:
  MPI.COMM_WORLD.Barrier()

grid = mmesh((reader.gmsh, filename), 2)
igrid = grid.hierarchicalGrid.interfaceGrid

size = [[2684], [1375], [], [703, 735]]
isize = [[100], [51], [], [27, 26]]

assert grid.size(0) in size[mpi_size-1]
assert igrid.size(0) in isize[mpi_size-1]
