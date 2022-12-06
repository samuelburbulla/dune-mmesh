from dune.grid import reader
from dune.mmesh import mmesh

from dune.mmesh.test.grids import line
grid = mmesh((reader.gmsh, line.filename), 2)
igrid = grid.hierarchicalGrid.interfaceGrid

print(grid.size(0), igrid.size(0))
assert(grid.size(0) == 2684)
assert(igrid.size(0) == 100)
