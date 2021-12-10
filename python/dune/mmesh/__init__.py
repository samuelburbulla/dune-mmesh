from ._grids import *
from ._move import *
from ._skeletontrace import *
from ._solve import *
from .utility import *

registry = dict()
registry["grid"] = grid_registry

from dune.fem import parameter
parameter.append({"fem.adaptation.method": "callback"})
