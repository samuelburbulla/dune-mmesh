## @example mmesh.py
#  This is an example of how to use MMesh with dune-python

import dune.grid
import dune.mmesh

domain = dune.grid.cartesianDomain([0, 0], [1, 1], [100, 100])
gridView = dune.mmesh.mmesh(domain)
