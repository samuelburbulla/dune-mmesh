from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, ConfigurationError

try:
    assertHave("HAVE_DUNE_MMESH")

    def mmesh(constructor, dimgrid=None, **parameters):
        from dune.grid.grid_generator import module, getDimgrid

        if not dimgrid:
            dimgrid = getDimgrid(constructor)

        if not (2 <= dimgrid and dimgrid <= 3):
            raise KeyError("Parameter error in MMesh with dimgrid=" + str(dimgrid) + ": dimgrid has to be either 2 or 3")

        typeName = "Dune::MovingMesh< " + str(dimgrid) + " >"
        includes = ["dune/mmesh/mmesh.hh"]
        gridModule = module(includes, typeName)

        return gridModule.LeafGrid(gridModule.reader(constructor))

    def mmeshInterfaceGrid(constructor, dimgrid=None, **parameters):
        from dune.grid.grid_generator import module, getDimgrid

        if not dimgrid:
            dimgrid = getDimgrid(constructor)

        if not (2 <= dimgrid and dimgrid <= 3):
            raise KeyError("Parameter error in MMesh with dimgrid=" + str(dimgrid) + ": dimgrid has to be either 2 or 3")

        typeName = "Dune::MMeshInterfaceGrid< Dune::MovingMesh< " + str(dimgrid) + " >, " + str(dimgrid) + " >"
        includes = ["dune/mmesh/interface/grid.hh"]
        gridModule = module(includes, typeName)

        return gridModule.LeafGrid(gridModule.reader(constructor))

    grid_registry = {
        "MMesh" : mmesh,
        "MMeshInterfaceGrid" : mmeshInterfaceGrid,
    }
except ConfigurationError:
    grid_registry = {}
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
