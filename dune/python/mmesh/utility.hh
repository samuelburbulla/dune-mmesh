#ifndef DUNE_PYTHON_MMESH_UTILITY
#define DUNE_PYTHON_MMESH_UTILITY

#if HAVE_DUNE_FEM

#include <dune/common/exceptions.hh>
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

    /////////////////////////
    // Interface indicator //
    /////////////////////////

    template <class BulkGV>
    struct InterfaceIndicator
    : public BindableGridFunctionWithSpace<FemPy::GridPart<BulkGV>, Dune::FieldVector<double, 1>>
    {
      using GridView = BulkGV;
      using GridPartType = FemPy::GridPart<BulkGV>;
      using Base = BindableGridFunctionWithSpace<GridPartType, Dune::FieldVector<double, 1>>;
      using SideGeometry = typename GridPartType::GridType::Intersection::LocalGeometry::Implementation;
      using RangeType = typename Base::RangeType;
      static constexpr bool scalar = true;

      InterfaceIndicator(const BulkGV &bulkGV)
      : Base(FemPy::gridPart<BulkGV>(bulkGV), "interfaceindicator", 0),
        onInterface_(false)
      {}

      void bind(const typename Base::EntityType &entity)
      {
        Base::bind(entity);
      }

      void bind(const typename BulkGV::Intersection &intersection, IntersectionSide side)
      {
        if( Base::gridPart().grid().isInterface( intersection ) )
          onInterface_ = true;
        else
          onInterface_ = false;
      }

      template <class Point>
      void evaluate(const Point &x, typename Base::RangeType &ret) const
      {
        if (onInterface_)
          ret = RangeType(1);
        else
          ret = RangeType(0);
      }

      template <class Point>
      void jacobian(const Point &x, typename Base::JacobianRangeType &ret) const
      {
        ret = typename Base::JacobianRangeType(0);
      }

      template <class Point>
      void hessian(const Point &x, typename Base::HessianRangeType &ret) const
      {
        ret = typename Base::HessianRangeType(0);
      }

    private:
      bool onInterface_;
    };

    template< class BulkGV >
    inline static void registerInterfaceIndicator(pybind11::module module, pybind11::class_< InterfaceIndicator<BulkGV> > cls)
    {
      using pybind11::operator""_a;

      cls.def( pybind11::init( [] ( const BulkGV &bulkGV ) {
          return InterfaceIndicator<BulkGV> ( bulkGV );
        } ), "bulkGV"_a, pybind11::keep_alive< 1, 2 >() );

      cls.def_property_readonly( "scalar", [] ( InterfaceIndicator<BulkGV> &self ) { return true; } );

      Dune::FemPy::registerGridFunction( module, cls );
    }



    /////////////
    // Normals //
    /////////////

    template <class InterfaceGV>
    struct Normals
    : public BindableGridFunctionWithSpace<FemPy::GridPart<InterfaceGV>, Dune::FieldVector<double, InterfaceGV::dimensionworld>>
    {
      static constexpr int dimw = InterfaceGV::dimensionworld;
      using GridView = InterfaceGV;
      using GridPartType = FemPy::GridPart<InterfaceGV>;
      using Base = BindableGridFunctionWithSpace<GridPartType, Dune::FieldVector<double, dimw>>;
      using RangeType = typename Base::RangeType;
      static constexpr bool scalar = false;

      Normals(const InterfaceGV &interfaceGV)
      : Base(FemPy::gridPart<InterfaceGV>(interfaceGV), "normals", 0)
      {}

      void bind(const typename Base::EntityType &entity)
      {
        Base::bind(entity);
        normal_ = Base::gridPart().grid().getMMesh().asIntersection(entity).centerUnitOuterNormal();
      }

    public:
      template <class Point>
      void evaluate(const Point &x, RangeType &ret) const
      {
        ret = normal_;
      }

      template <class Point>
      void jacobian(const Point &x, typename Base::JacobianRangeType &ret) const
      {
        ret = typename Base::JacobianRangeType(0);
      }

      template <class Point>
      void hessian(const Point &x, typename Base::HessianRangeType &ret) const
      {
        ret = typename Base::HessianRangeType(0);
      }

    private:
      Dune::FieldVector<double, dimw> normal_;
    };

    template< class InterfaceGV >
    inline static void registerNormals(pybind11::module module, pybind11::class_< Normals<InterfaceGV> > cls)
    {
      using pybind11::operator""_a;

      cls.def( pybind11::init( [] ( const InterfaceGV &interfaceGV ) {
          return Normals<InterfaceGV> ( interfaceGV );
        } ), "interfaceGV"_a, pybind11::keep_alive< 1, 2 >() );

      cls.def_property_readonly( "scalar", [] ( Normals<InterfaceGV> &self ) { return false; } );

      Dune::FemPy::registerGridFunction( module, cls );
    }



    //////////////////////
    // DistanceFunction //
    //////////////////////

    template <class GV>
    struct DistanceFunction
    : public BindableGridFunctionWithSpace<FemPy::GridPart<GV>, Dune::FieldVector<double, 1>>
    {
      using GridView = GV;
      static constexpr int dim = GridView::dimensionworld;
      using GridPartType = FemPy::GridPart<GridView>;
      using Base = BindableGridFunctionWithSpace<GridPartType, Dune::FieldVector<double, 1>>;
      using RangeType = typename Base::RangeType;
      static constexpr bool scalar = true;

      DistanceFunction(const GridView &gridView)
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
        auto xLocal = Dune::Fem::coordinate(x);
        auto jacT = interpolation_.jacobianTransposed(xLocal);
        for (std::size_t i = 0; i < dim; ++i)
          ret[0][i] = jacT[i][0];
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
    inline static void registerDistanceFunction(pybind11::module module, pybind11::class_< DistanceFunction<GridView> > cls)
    {
      using pybind11::operator""_a;

      cls.def( pybind11::init( [] ( const GridView &gridView ) {
          return DistanceFunction<GridView> ( gridView );
        } ), "gridView"_a, pybind11::keep_alive< 1, 2 >() );

      cls.def_property_readonly( "scalar", [] ( DistanceFunction<GridView> &self ) { return true; } );

      Dune::FemPy::registerGridFunction( module, cls );
    }

  }  // namespace Fem

} // namespace Dune

#endif // HAVE_DUNE_FEM

#endif // DUNE_PYTHON_MMESH_UTILITY
