"""
.. module:: _grids
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, ConfigurationError
from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

mmGenerator = SimpleGenerator("HierarchicalGrid", "Dune::Python::MMGrid")
mmifGenerator = SimpleGenerator("HierarchicalGrid", "Dune::Python::MMIFGrid")
def mmesh(constructor, dimgrid=None, **parameters):
    """Create an MMesh grid.

    Args:
        constructor: Grid constructor, e.g. (dune.grid.reader.gmsh, 'grid.msh').
        dimgrid (int, optional): The world dimension (must be 2 or 3). Might be guessed by constructor.
        **parameters: Additional parameters.

    Returns:
        An MMesh grid instance.
    """

    from dune.grid.grid_generator import module, getDimgrid

    if not dimgrid:
        dimgrid = getDimgrid(constructor)

    if not (2 <= dimgrid and dimgrid <= 3):
        raise KeyError("Parameter error in MMesh with dimgrid=" + str(dimgrid) + ": dimgrid has to be either 2 or 3")

    includes = ["dune/mmesh/mmesh.hh", "dune/python/mmesh/interfacegrid.hh"]
    typeName = "Dune::MovingMesh< " + str(dimgrid) + " >"
    gridModule = module(includes, typeName, generator=mmGenerator)

    typeName = "typename Dune::MovingMesh< " + str(dimgrid) + " >::InterfaceGrid"
    igridModule = module(includes, typeName, generator=mmifGenerator)

    gridView = gridModule.LeafGrid(gridModule.reader(constructor))

    return gridView

grid_registry = {
    "MMesh" : mmesh
}

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
