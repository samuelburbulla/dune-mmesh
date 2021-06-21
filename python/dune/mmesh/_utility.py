import logging, traceback
logger = logging.getLogger(__name__)

import io
from dune.generator import algorithm
from dune.generator.exceptions import CompileError

def interfaceIndicator(igrid, grid=None, restrict=True):
    """Return indicator of interface edges.

    Args:
        igrid: The interfacegrid.
        grid (Grid, optional): The bulk grid. Necessary, if wrapped.
        restrict (bool, optional): If True, the returned UFL expression is restricted.

    Returns:
        Skeleton function that is one at interface edges.
    """
    from ufl import avg
    from dune.mmesh import skeleton
    try:
        one = igrid.hierarchicalGrid.one
    except:
        from dune.fem.space import finiteVolume
        space = finiteVolume(igrid)
        one = space.interpolate(1, name="one")
        igrid.hierarchicalGrid.one = one

    if restrict:
        return avg(skeleton(one, grid=grid))
    else:
        return skeleton(one, grid=grid)


def domainMarker(grid):
    """Return domain markers passed by .msh grid file.

    Args:
       grid: The grid.

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
    try:
        cppfunc = cppFunction(grid, name="domainMarker", order=0, fctName="domainMarker", includes=io.StringIO(code), args=[grid])
    except CompileError as e:
        code = code.replace("impl()", "impl().hostEntity().impl()")
        cppfunc = cppFunction(grid, name="domainMarker", order=0, fctName="domainMarker", includes=io.StringIO(code), args=[grid])
    fvspace = finiteVolume(grid)
    return fvspace.interpolate(cppfunc, name="domainMarker")

def normals(igrid):
    """Return normal vectors to the interface grid elements.

    Args:
       igrid: The interface grid.

    Returns:
       Grid function on the interface grid.
    """
    code="""
    #include <functional>
    template <class IGV>
    auto normal(const IGV &igv, bool minus=false) {
      return [&igv, minus] (const auto& entity, const auto& xLocal) mutable -> auto {
        return igv.grid().getMMesh().asIntersection( entity ).centerUnitOuterNormal() * (minus ? -1.0 : 1.0);
      };
    }
    """
    import dune.ufl
    from dune.fem.function import cppFunction
    n_p = cppFunction(igrid, name="normal_p", order=0, fctName="normal", includes=io.StringIO(code), args=[igrid, True])
    n_m = cppFunction(igrid, name="normal_m", order=0, fctName="normal", includes=io.StringIO(code), args=[igrid, False])
    predefined = {}
    predefined[n_p('+')] = n_p
    predefined[n_p('-')] = n_m
    n_p.predefined = predefined
    return n_p


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
