from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, ConfigurationError
from dune.generator.generator import SimpleGenerator
from dune.grid.grid_generator import addAttr
from dune.common.hashit import hashIt

try:
    assertHave("HAVE_DUNE_MMESH")

    igGenerator = SimpleGenerator("InterfaceGrid", "Dune::Python")
    def interfaceGrid(grid):
        includes = grid._includes + ["dune/mmesh/interface/grid.hh", "dune/python/mmesh/interfacegrid.hh"]
        typeName = "typename " + grid._typeName + "::Grid"
        moduleName = "interfacegrid_" + hashIt(typeName)
        module = igGenerator.load(includes, typeName, moduleName)
        addAttr( module, module.LeafGrid );

        return grid.hierarchicalGrid._interfaceGrid

    def mmesh(constructor, dimgrid=None, **parameters):
        from dune.grid.grid_generator import module, getDimgrid

        if not dimgrid:
            dimgrid = getDimgrid(constructor)

        if not (2 <= dimgrid and dimgrid <= 3):
            raise KeyError("Parameter error in MMesh with dimgrid=" + str(dimgrid) + ": dimgrid has to be either 2 or 3")

        includes = ["dune/mmesh/mmesh.hh"]
        typeName = "Dune::MovingMesh< " + str(dimgrid) + " >"
        gridModule = module(includes, typeName)

        gridModule.LeafGrid.interfaceGrid = property( interfaceGrid )

        return gridModule.LeafGrid(gridModule.reader(constructor))

    grid_registry = {
        "MMesh" : mmesh
    }
except ConfigurationError:
    grid_registry = {}
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
