## @example mmesh.py
#  This is an example of how to use MMesh with dune-python

import dune.grid
import dune.mmesh

domain = dune.grid.cartesianDomain([0, 0], [1, 1], [10, 10])
gridView = dune.mmesh.mmesh(domain)
assert gridView.size(0) == 200
gridView.writeVTK("test-python-mmesh")
