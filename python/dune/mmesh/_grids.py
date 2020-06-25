from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, ConfigurationError
from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

try:
    assertHave("HAVE_DUNE_MMESH")

    mmGenerator = SimpleGenerator("HierarchicalGrid", "Dune::Python::MMGrid")
    mmifGenerator = SimpleGenerator("HierarchicalGrid", "Dune::Python::MMIFGrid")
    def mmesh(constructor, dimgrid=None, **parameters):
        from dune.grid.grid_generator import module, getDimgrid, addAttr

        if not dimgrid:
            dimgrid = getDimgrid(constructor)

        if not (2 <= dimgrid and dimgrid <= 3):
            raise KeyError("Parameter error in MMesh with dimgrid=" + str(dimgrid) + ": dimgrid has to be either 2 or 3")

        includes = ["dune/mmesh/mmesh.hh", "dune/python/mmesh/interfacegrid.hh"]
        typeName = "Dune::MovingMesh< " + str(dimgrid) + " >"
        gridModule = module(includes, typeName, generator=mmGenerator)

        typeName = "typename Dune::MovingMesh< " + str(dimgrid) + " >::InterfaceGrid"
        igridModule = module(includes, typeName, generator=mmifGenerator)

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
