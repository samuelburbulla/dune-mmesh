import io
from dune.generator import algorithm
from dune.fem.function import cppFunction

################################################################################
# Compute trace function
################################################################################
def traceFunction(igridView, bulkGridFunction, name="trace", order=1, inside=True):
    code="""
    #include <functional>
    template <class GV, class IGV, class BulkGF>
    auto trace(const IGV &igv, const BulkGF &bulkgf, bool inside) {
      using IElement = typename IGV::template Codim<0>::Entity;
      using LocalCoordinate = typename IElement::Geometry::LocalCoordinate;
      if (inside)
      {
        std::function<typename BulkGF::RangeType(const IElement&,LocalCoordinate)>
          ret = [&bulkgf, lbulk=localFunction(bulkgf)]
          (const auto& interfaceEntity,const auto& xLocal) mutable -> auto {
          const auto intersection =
                bulkgf.space().gridPart().grid().asIntersection( interfaceEntity );
          lbulk.bind(intersection.inside()); // lambda must be mutable so that the non const function can be called
          return lbulk(intersection.geometryInInside().global(xLocal));
        };
        return ret;
      }
      else
      {
        std::function<typename BulkGF::RangeType(const IElement&,LocalCoordinate)>
          ret = [&bulkgf, lbulk=localFunction(bulkgf)]
          (const auto& interfaceEntity,const auto& xLocal) mutable -> auto {
          typename BulkGF::RangeType ret(0);
          const auto intersection =
                bulkgf.space().gridPart().grid().asIntersection( interfaceEntity );
          if (intersection.neighbor())
          {
            lbulk.bind(intersection.outside()); // lambda must be mutable so that the non const function can be called
            ret = lbulk(intersection.geometryInOutside().global(xLocal));
          }
          return ret;
        };
        return ret;
      }
    }
    """
    return cppFunction(igridView, name=name, order=order, fctName="trace", includes=io.StringIO(code), args=[igridView, bulkGridFunction, inside])
################################################################################


################################################################################
# Obtain domain markers
################################################################################
def domainMarker(gridView):
    code="""
    #include <functional>
    template <class GV>
    auto domainMarker(const GV &gv) {
      auto
        ret = [] (const auto& entity, const auto& xLocal) mutable -> auto {
        return entity.impl().domainMarker();
      };
      return ret;
    }
    """
    return cppFunction(gridView, name="domainMarker", order=0, fctName="domainMarker", includes=io.StringIO(code), args=[gridView])
################################################################################


################################################################################
# Move interface and restore domain markers
################################################################################
def moveInterface(hgrid, movedf):
    assert hgrid.dimension == hgrid.dimensionworld and hasattr(hgrid, 'leafView'), \
        "please pass the bulk hierarchical grid as first argument of moveInterface()"

    code="""
    #include <functional>
    #include <iostream>
    #include <unordered_map>
    #include <dune/mmesh/mmesh.hh>
    #include <dune/grid/utility/persistentcontainer.hh>

    template <class Grid, class DF>
    void moveInterface(Grid &grid, const DF &movedf) {
      using Element = typename Grid::LeafGridView::template Codim<0>::Entity;
      using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
      using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;
      using PersistentContainer = Dune::PersistentContainer<Grid, std::size_t>;

      static constexpr int dim = Grid::dimension;
      const auto& gridView = grid.leafGridView();
      const auto& igrid = grid.interfaceGrid();
      const auto& igridView = igrid.leafGridView();
      const auto& idSet = igrid.globalIdSet();

      using ShiftContainer = std::unordered_map<std::size_t, GlobalCoordinate>;
      ShiftContainer shiftContainer;

      auto getShifts = [&](bool init = false) {
          std::vector<GlobalCoordinate> shifts( igridView.size(dim-1) );
          if (init)
          {
              for( const auto& e : elements( igridView ) )
                for( std::size_t i = 0; i < e.subEntities(dim-1); ++i )
                {
                  const auto& v = e.template subEntity<dim-1>( i );
                  GlobalCoordinate m;
                  movedf.evaluate( v.geometry().center(), m );
                  shifts[ igridView.indexSet().index(v) ] = m;
                  shiftContainer[ idSet.id(v) ] = m;
                }
          }
          else
          {
              for( const auto v : vertices( igridView ) )
              {
                if (shiftContainer.count( idSet.id(v) ) > 0)
                {
                  shifts[ igridView.indexSet().index(v) ] = shiftContainer[ idSet.id(v) ];
                }
                else
                {
                  // interpolate from neighbors
                  GlobalCoordinate shift;
                  int count = 0;
                  for ( const auto vi : incidentInterfaceVertices( v ) )
                    if (shiftContainer.count( idSet.id(vi) ) > 0 )
                    {
                      shift += shiftContainer[ idSet.id(vi) ];
                      count++;
                    }
                  assert( count > 0 );
                  shift /= count;
                  shifts[ igridView.indexSet().index(v) ] = shift;
                }
              }
          }
          return shifts;
      };
      getShifts( true );

      bool adapt = true;
      int count = 25;
      do
      {
          if (count > 0)
            grid.markElements( getShifts() );

          grid.ensureInterfaceMovement( getShifts() );

          if (!grid.preAdapt())
            break;

          PersistentContainer container( grid, 0 );

          for( const auto& e : elements( gridView ) )
            // if( e.mightVanish() )
            {
              container[ e ] = e.impl().domainMarker();
            }

          adapt = grid.adapt();

          for( const auto& e : elements( gridView ) )
            // if( e.isNew() )
            if( grid.hasConnectedComponent(e) )
            {
              const auto evol = e.geometry().volume();
              double marker = 0;

              auto component = grid.getConnectedComponent(e);

              if (evol < 1e-14)
                marker = container[ *component.children().begin() ];
              else
              {
                for ( const auto elem : component.children() )
                  marker += elem.intersectionVolume( e ) / evol * container[ elem ];
              }

              e.impl().hostEntity()->info().domainMarker = (std::size_t) std::round(marker);
            }

          grid.postAdapt();

          count--;
      } while( adapt );

      grid.moveInterface( getShifts() );
    }
    """
    algorithm.run("moveInterface", io.StringIO(code), hgrid, movedf)
################################################################################


################################################################################
# Obtain curvature
################################################################################
def interfaceCurvature(igrid, atVertex=True):
    code="""
    #include <functional>
    #include <dune/grid/common/mcmgmapper.hh>
    #include <dune/mmesh/interface/curvatureoperator.hh>
    #include <dune/geometry/multilineargeometry.hh>
    template <class GV>
    auto curvature(const GV &gv, const bool atVertex) {
        static constexpr int dim = GV::dimension;
        using Element = typename GV::template Codim<0>::Entity;
        using LocalCoordinate = typename Element::Geometry::LocalCoordinate;
        using GlobalCoordinate = typename Element::Geometry::GlobalCoordinate;

        using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GV>;
        Mapper mapper(gv, atVertex ? Dune::mcmgVertexLayout() : Dune::mcmgElementLayout());

        std::vector<double> curvatures(mapper.size());
        std::vector<GlobalCoordinate> centers(mapper.size());

        if (atVertex) {
            using CurvOp = Dune::CurvatureOperator<GV, Mapper, Dune::CurvatureLayout::Vertex>;
            CurvOp curvOp(gv, mapper);
            curvOp(curvatures, centers);

            std::function<double(const Element&, LocalCoordinate)>
                ret = [mapper, curvatures] (const Element& entity, const LocalCoordinate& xLocal) mutable -> double {
                    std::vector<Dune::FieldVector<double, 1>> values;
                    for ( std::size_t i = 0; i < entity.subEntities(dim); ++i )
                      values.push_back( curvatures[ mapper.index( entity.template subEntity<dim>(i) ) ] );

                    const Dune::MultiLinearGeometry<double, dim, 1> interpolation(entity.type(), values);
                    return interpolation.global(xLocal);
                };
            return ret;
        } else {
            using CurvOp = Dune::CurvatureOperator<GV, Mapper, Dune::CurvatureLayout::Element>;
            CurvOp curvOp(gv, mapper);
            curvOp(curvatures, centers);

            std::function<double(const Element&, LocalCoordinate)>
                ret = [mapper, curvatures] (const Element& entity, const LocalCoordinate& xLocal) mutable -> double {
                    return curvatures[ mapper.index( entity ) ];
                };
            return ret;
        }
    }
    """
    cppFunc = cppFunction(igrid, name="curvature", order=(1 if atVertex else 0), fctName="curvature", includes=io.StringIO(code), args=[igrid, atVertex])

    if atVertex:
        from dune.fem.space import lagrange
        space = lagrange(igrid)
    else:
        from dune.fem.space import finiteVolume
        space = finiteVolume(igrid)

    return space.interpolate(cppFunc, name="curvature")
################################################################################
