import logging, traceback
logger = logging.getLogger(__name__)

import io
from dune.generator import algorithm, builder
from dune.common.hashit import hashIt
from dune.generator.exceptions import CompileError

def interfaceIndicator(igrid, grid=None, restrict=True):
    """Return indicator of interface edges.

    Args:
        igrid: The interfacegrid.
        grid: The bulk grid. Pass this parameter if the bulk grid is wrapped, e.g., by a geometryGridView.
        restrict (bool, optional): If True, the returned UFL expression is restricted.

    Returns:
        Skeleton function that is one at interface edges.
    """

    # Pre-compiled version
    if grid is None:
        grid = igrid.hierarchicalGrid.bulkGrid
        if igrid.dimension == 1:
            import dune.mmesh._utility2d as module
        else:
            import dune.mmesh._utility3d as module

    # JIT version for wrapped grid
    else:
        moduleName = "interfaceindicator_" + hashIt(grid.cppTypeName)
        signature = "Compiling InterfaceIndicator"
        source = """
#include <config.h>
#include <dune/mmesh/mmesh.hh>
#include <dune/python/mmesh/utility.hh>
#include <dune/python/grid/hierarchical.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>
"""
        includes = sorted(set(grid.cppIncludes))
        source += "".join(["#include <" + i + ">\n" for i in includes])
        source += """
PYBIND11_MODULE( """ + moduleName + """, module )
{
  // InterfaceIndicator
  auto clsIndicator = Dune::Python::insertClass< Dune::Fem::InterfaceIndicator< typename """ + grid.cppTypeName + """ > >(
    module,
    "InterfaceIndicator",
    Dune::Python::GenerateTypeName("Dune::Fem::InterfaceIndicator< typename """ + grid.cppTypeName + """ >"),
    Dune::Python::IncludeFiles{"dune/mmesh/mmesh.hh", "dune/python/grid/hierarchical.hh", "dune/python/mmesh/utility.hh"}
  ).first;
  Dune::Fem::registerInterfaceIndicator( module, clsIndicator );
}
"""
        module = builder.load(moduleName, source, signature)

    indicator = module.InterfaceIndicator(grid)
    from dune.ufl import GridFunction
    indicator = GridFunction(indicator)

    if restrict:
        from ufl import avg
        return avg(indicator)
    else:
        return indicator


def normals(igrid):
    """Return normal vectors to the interface grid elements.

    Args:
        igrid: The interface grid.

    Returns:
        Grid function on the interface grid. Coincides with n('+') of the bulk facet normal.
    """
    igrid = igrid.hierarchicalGrid.leafView

    if igrid.dimension == 1:
        import dune.mmesh._utility2d as module
    else:
        import dune.mmesh._utility3d as module

    normals = module.Normals(igrid)
    from dune.ufl import GridFunction
    normals = GridFunction(normals)
    return normals


def distance(grid):
    """Return function representing the distance to the interface

    Args:
        grid: The grid.

    Returns:
        Piecewise linear grid function.
    """

    if grid.dimension == 2:
        import dune.mmesh._utility2d as module
    else:
        import dune.mmesh._utility3d as module

    distance = module.Distance(grid.hierarchicalGrid.leafView)
    from dune.ufl import GridFunction
    distance = GridFunction(distance)
    return distance


def domainMarker(grid, wrapped=False):
    """Return domain markers passed by .msh grid file.

    Args:
       grid: The grid.
       wrapped: Set this to true if grid is wrapped.

    Returns:
       Grid function with cell-wise markers from .msh file.
    """
    code = """
    #include <functional>
    template <class GV>
    auto domainMarker(const GV &gv) {
      auto ret = [] (const auto& entity, const auto& xLocal) mutable -> auto {
        return entity.impl().domainMarker();
      };
      return ret;
    }
    """
    from dune.fem.space import finiteVolume
    from dune.fem.function import cppFunction
    if not wrapped:
        cppfunc = cppFunction(grid, name="domainMarker", order=0, fctName="domainMarker", includes=io.StringIO(code), args=[grid])
    else:
        code = code.replace("impl()", "impl().hostEntity().impl()")
        cppfunc = cppFunction(grid, name="domainMarkerImpl", order=0, fctName="domainMarker", includes=io.StringIO(code), args=[grid])
    fvspace = finiteVolume(grid)
    return fvspace.interpolate(cppfunc, name="domainMarker")


def interfaceDomainMarker(igrid):
    """Return interface domain markers passed by .msh grid file.

    Args:
       igrid: The interface grid.

    Returns:
       Interface grid function with cell-wise markers from .msh file.
    """
    code = """
    #include <iostream>
    #include <functional>
    template <class GV>
    auto interfaceDomainMarker(const GV &igv) {
      auto ret = [&igv] (const auto& entity, const auto& xLocal) mutable -> auto {
        return igv.grid().domainMarker(entity);
      };
      return ret;
    }
    """
    from dune.fem.function import cppFunction
    return cppFunction(igrid, name="interfaceDomainMarker", order=0, fctName="interfaceDomainMarker", includes=io.StringIO(code), args=[igrid])


def edgeMovement(grid, shifts):
    """Return linear interpolation of vertex shifts.

    Args:
       grid: The grid.
       shifts:   List of shifts.

    Returns:
       Vector-valued grid function.
    """
    code="""
    #include <functional>
    #include <dune/geometry/multilineargeometry.hh>
    #include <dune/python/pybind11/pybind11.h>
    #include <dune/python/pybind11/numpy.h>

    template <class GV>
    auto edgeMovement(const GV &gv, const pybind11::array_t<double>& input) {
        static constexpr int dim = GV::dimension;
        using GlobalCoordinate = Dune::FieldVector<double, dim>;

        // obtain shifts from buffer protocol
        pybind11::buffer_info buffer = input.request();
        double *ptr = (double *) buffer.ptr;
        int n = buffer.shape[0];

        std::vector<GlobalCoordinate> shifts( n );
        for ( std::size_t i = 0; i < n; ++i )
          for ( int d = 0; d < dim; ++d )
            shifts[i][d] = ptr[dim*i+d];

        const auto& igrid = gv.grid().interfaceGrid();
        const auto& iindexSet = igrid.leafIndexSet();

        auto ret = [shifts, &igrid, &iindexSet] (const auto& entity, const auto& xLocal) mutable -> auto {
            // return a linear interpolation of the vertex shift values
            std::vector<GlobalCoordinate> shift( dim+1 );
            for ( std::size_t i = 0; i < entity.subEntities( dim ); ++i )
            {
                const auto& vertex = entity.template subEntity<dim>( i );
                if ( vertex.impl().isInterface() )
                {
                    // cast to interface grid vertex to obtain index
                    const auto ivertex = igrid.entity( vertex.impl().hostEntity() );
                    shift[i] = shifts[ iindexSet.index( ivertex ) ];
                }
            }
            const Dune::MultiLinearGeometry<double, dim, dim> interpolation(entity.type(), shift);
            return interpolation.global(xLocal);
        };
        return ret;
    }
    """
    from dune.fem.function import cppFunction
    cppFunc = cppFunction(grid, name="edgeMovement", order=1, fctName="edgeMovement", includes=io.StringIO(code), args=[grid, shifts])
    from dune.fem.space import lagrange
    space = lagrange(grid, dimRange=grid.dimWorld, order=1)
    return space.interpolate(cppFunc, name="edgeMovement")


def interfaceEdgeMovement(igrid, shifts):
    """Return linear interpolation of vertex shifts on the interface.

    Args:
       igrid: The interface grid.
       shifts:    List of shifts.

    Returns:
       Vector-valued interface grid function.
    """
    code="""
    #include <functional>
    #include <dune/geometry/multilineargeometry.hh>
    #include <dune/python/pybind11/pybind11.h>
    #include <dune/python/pybind11/numpy.h>

    template <class GV>
    auto interfaceEdgeMovement(const GV &gv, const pybind11::array_t<double>& input) {
        static constexpr int dim = GV::dimension;
        static constexpr int dimw = GV::dimensionworld;
        using GlobalCoordinate = Dune::FieldVector<double, dimw>;

        // obtain shifts from buffer protocol
        pybind11::buffer_info buffer = input.request();
        double *ptr = (double *) buffer.ptr;
        int n = buffer.shape[0];

        std::vector<GlobalCoordinate> shifts( n );
        for ( std::size_t i = 0; i < n; ++i )
          for ( int d = 0; d < dimw; ++d )
            shifts[i][d] = ptr[dimw*i+d];

        const auto& igrid = gv.grid();
        const auto& iindexSet = igrid.leafIndexSet();

        auto ret = [shifts, &igrid, &iindexSet] (const auto& entity, const auto& xLocal) mutable -> auto {
            // return a linear interpolation of the vertex shift values
            std::vector<GlobalCoordinate> shift( dim+1 );
            for ( std::size_t i = 0; i < entity.subEntities( dim ); ++i )
            {
                const auto& vertex = entity.template subEntity<dim>( i );
                shift[i] = shifts[ iindexSet.index( vertex ) ];
            }
            const Dune::MultiLinearGeometry<double, dim, dimw> interpolation(entity.type(), shift);
            return interpolation.global(xLocal);
        };
        return ret;
    }
    """
    from dune.fem.function import cppFunction
    cppFunc = cppFunction(igrid, name="interfaceEdgeMovement", order=1, fctName="interfaceEdgeMovement", includes=io.StringIO(code), args=[igrid, shifts])
    from dune.fem.space import lagrange
    space = lagrange(igrid, dimRange=igrid.dimWorld, order=1)
    return space.interpolate(cppFunc, name="interfaceEdgeMovement")
