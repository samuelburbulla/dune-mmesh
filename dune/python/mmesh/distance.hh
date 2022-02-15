#ifndef DUNE_PYTHON_MMESH_DISTANCE
#define DUNE_PYTHON_MMESH_DISTANCE

#if HAVE_DUNE_FEM

#include <dune/common/exceptions.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dune/python/common/typeregistry.hh>
#include <dune/fem/common/intersectionside.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/grid/function.hh>

namespace Dune
{

  namespace Fem
  {

    //////////////
    // Distance //
    //////////////

    template <class GV>
    struct Distance
    : public BindableGridFunctionWithSpace<FemPy::GridPart<GV>, Dune::FieldVector<double, 1>>
    {
      using GridView = GV;
      static constexpr int dim = GridView::dimensionworld;
      using GridPartType = FemPy::GridPart<GridView>;
      using Base = BindableGridFunctionWithSpace<GridPartType, Dune::FieldVector<double, 1>>;
      using RangeType = typename Base::RangeType;
      static constexpr bool scalar = true;

      Distance(const GridView &gridView)
      : Base(FemPy::gridPart<GridView>(gridView), "distance", 0),
        interpolation_(GeometryTypes::simplex(dim), std::vector<Dune::FieldVector<double, 1>>())
      {}

      void bind(const typename Base::EntityType &entity)
      {
        Base::bind(entity);

        const auto& distance = this->gridPart().grid().indicator().distance();

        std::vector<Dune::FieldVector<double, 1>> distances;
        distances.resize(dim+1);
        for ( std::size_t i = 0; i < entity.subEntities( dim ); ++i )
        {
            const auto& vertex = entity.template subEntity<dim>( i );
            distances[i] = distance(vertex);
        }

        interpolation_ = Dune::MultiLinearGeometry<double, dim, 1> (entity.type(), distances);
      }

    public:
      template <class Point>
      void evaluate(const Point &x, RangeType &ret) const
      {
        auto xLocal = Dune::Fem::coordinate(x);
        ret = interpolation_.global(xLocal);
      }

      template <class Point>
      void jacobian(const Point &x, typename Base::JacobianRangeType &ret) const
      {
        const auto geo = this->entity().geometry();

        Dune::FieldVector<double, dim> p;
        RangeType dt0;
        auto t0 = geo.local(p);
        evaluate(t0, dt0);

        RangeType dt;
        for (std::size_t i = 0; i < dim; ++i)
        {
          p[i] = 1.;
          auto t = geo.local(p);
          p[i] = 0.;
          evaluate(t, dt);
          ret[0][i] = dt - dt0;
        }
      }

      template <class Point>
      void hessian(const Point &x, typename Base::HessianRangeType &ret) const
      {
        ret = typename Base::HessianRangeType(0);
      }

    private:
      Dune::MultiLinearGeometry<double, dim, 1> interpolation_;
    };

    template< class GridView >
    inline static void registerDistance(pybind11::module module, pybind11::class_< Distance<GridView> > cls)
    {
      using pybind11::operator""_a;

      cls.def( pybind11::init( [] ( const GridView &gridView ) {
          return Distance<GridView> ( gridView );
        } ), "gridView"_a, pybind11::keep_alive< 1, 2 >() );

      cls.def_property_readonly( "scalar", [] ( Distance<GridView> &self ) { return true; } );

      Dune::FemPy::registerGridFunction( module, cls );
    }

  }  // namespace Fem

} // namespace Dune

#endif // HAVE_DUNE_FEM

#endif // DUNE_PYTHON_MMESH_DISTANCE
