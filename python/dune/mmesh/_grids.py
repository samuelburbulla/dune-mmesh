from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

from dune.common.checkconfiguration import assertHave, ConfigurationError
from dune.generator.generator import SimpleGenerator
from dune.common.hashit import hashIt

try:
    assertHave("HAVE_DUNE_MMESH")

    # export normals to the interface elements
    def normal(igridView):
        code="""
        #include <functional>
        template <class IGV>
        auto normal(const IGV &igv, bool minus=false) {
          return [&igv, minus] (const auto& entity, const auto& xLocal) mutable -> auto {
            return igv.grid().getMMesh().asIntersection( entity ).centerUnitOuterNormal() * (minus ? -1.0 : 1.0);
          };
        }
        """

        import io
        import dune.ufl
        from dune.fem.function import cppFunction
        cppFunc_p = cppFunction(igridView, name="normal_p", order=0, fctName="normal", includes=io.StringIO(code), args=[igridView, True])
        cppFunc_m = cppFunction(igridView, name="normal_m", order=0, fctName="normal", includes=io.StringIO(code), args=[igridView, False])
        n_p = dune.ufl.GridFunction( cppFunc_p )
        n_m = dune.ufl.GridFunction( cppFunc_m )
        predefined = {}
        predefined[n_p('+')] = n_p
        predefined[n_p('-')] = n_m
        n_p.predefined = predefined
        return n_p


    mmGenerator = SimpleGenerator("HierarchicalGrid", "Dune::Python::MMGrid")
    mmifGenerator = SimpleGenerator("HierarchicalGrid", "Dune::Python::MMIFGrid")
    def mmesh(constructor, dimgrid=None, **parameters):
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
        igridModule.LeafGrid.normal = normal(gridView.hierarchicalGrid.interfaceGrid)

        return gridView

    grid_registry = {
        "MMesh" : mmesh
    }
except ConfigurationError:
    grid_registry = {}
    pass

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
