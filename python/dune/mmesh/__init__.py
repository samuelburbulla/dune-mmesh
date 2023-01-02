"""The main dune mmesh module file.

This file is imports all internal implementations.
"""

from ._grids import *
from ._move import *
from ._skeletontrace import *
from ._solve import *
from ._utility import *

registry = {}
registry["grid"] = grid_registry

try:
  from dune.fem import parameter
  parameter.append({"fem.adaptation.method": "callback"})
except ImportError:
  pass
