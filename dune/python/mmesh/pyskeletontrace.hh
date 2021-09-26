#ifndef DUNE_MMESH_PYSKELETONTRACE
#define DUNE_MMESH_PYSKELETONTRACE

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

    ///////////////////////
    // Skeleton function //
    ///////////////////////

    template <class BulkGV, class InterfaceGridFunction>
    struct SkeletonGF
    : public BindableGridFunctionWithSpace<FemPy::GridPart<BulkGV>,
        typename InterfaceGridFunction::RangeType>
    {
      using GridView = BulkGV;
      using GridPart = FemPy::GridPart<BulkGV>;
      using Base = BindableGridFunctionWithSpace<GridPart, typename InterfaceGridFunction::RangeType>;
      using ILocalFunction = ConstLocalFunction<InterfaceGridFunction>;
      using SideGeometry = typename GridPart::GridType::Intersection::LocalGeometry::Implementation;
      static constexpr bool scalar = (InterfaceGridFunction::RangeType::dimension == 1);

      SkeletonGF(const BulkGV &bulkGV, const InterfaceGridFunction &igf)
      : Base(FemPy::gridPart<BulkGV>(bulkGV), "skeleton_"+igf.name(), igf.order()),
        ilf_(igf), onInterface_(false)
      {}

      void bind(const typename Base::EntityType &entity)
      {
        Base::bind(entity);
      }

      void bind(const typename BulkGV::Intersection &intersection,
                IntersectionSide side)
      {
        if( Base::gridPart().grid().isInterface( intersection ) )
        {
          onInterface_ = true;
          ilf_.bind( Base::gridPart().grid().asInterfaceEntity( intersection ) );
          sideGeometry_ = (side == IntersectionSide::in)
            ? intersection.geometryInInside().impl() : intersection.geometryInOutside().impl();
        }
        else onInterface_ = false;
      }

    public:
      template <class Point>
      void evaluate(const Point &x, typename Base::RangeType &ret) const
      {
        if (onInterface_)
        {
          auto xLocal = Dune::Fem::coordinate(x);
          auto ix     = sideGeometry_.local( xLocal );
          ilf_.evaluate(ix,ret);
        } else ret = typename Base::RangeType(0);
      }

      template <class Point>
      void jacobian(const Point &x, typename Base::JacobianRangeType &ret) const
      {
        if (onInterface_)
        {
          auto xLocal = Dune::Fem::coordinate(x);
          auto ix     = sideGeometry_.local( xLocal );
          ilf_.jacobian(ix,ret);
        } else ret = typename Base::JacobianRangeType(0);
      }

      template <class Point>
      void hessian(const Point &x, typename Base::HessianRangeType &ret) const
      {
        if (onInterface_)
        {
          auto xLocal = Dune::Fem::coordinate(x);
          auto ix     = sideGeometry_.local( xLocal );
          ilf_.hessian(ix,ret);
        } else ret = typename Base::HessianRangeType(0);
      }

    private:
      ILocalFunction ilf_;
      bool onInterface_;
      SideGeometry sideGeometry_;
    };

    template< class BulkGV, class InterfaceGridFunction >
    inline static void registerSkeletonGF(pybind11::module module, pybind11::class_< SkeletonGF<BulkGV, InterfaceGridFunction> > cls)
    {
      using pybind11::operator""_a;

      cls.def( pybind11::init( [] ( const BulkGV &bulkGV, const InterfaceGridFunction &bgf ) {
          return SkeletonGF<BulkGV, InterfaceGridFunction> ( bulkGV, bgf );
        } ), "bulkGV"_a, "igf"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );

      cls.def_property_readonly( "scalar", [] ( SkeletonGF<BulkGV, InterfaceGridFunction> &self ) { return self.scalar; } );

      Dune::FemPy::registerGridFunction( module, cls );
    }


    ////////////////////
    // Trace function //
    ////////////////////

    template <class IGV, class BulkGridFunction, Dune::Fem::IntersectionSide side>
    struct TraceGF
    : public Dune::Fem::BindableGridFunctionWithSpace<Dune::FemPy::GridPart<IGV>,
                   typename BulkGridFunction::RangeType>
    {
      using GridView = IGV;
      using GridPart = Dune::FemPy::GridPart<IGV>;
      using Base = Dune::Fem::BindableGridFunctionWithSpace<GridPart,typename BulkGridFunction::RangeType>;
      using BLocalFunction = Dune::Fem::ConstLocalFunction<BulkGridFunction>;
      using BulkGridPart = typename BulkGridFunction::GridPartType;
      using SideGeometry = typename GridPart::GridType::MMeshType::Intersection::LocalGeometry::Implementation;
      static constexpr bool scalar = (BulkGridFunction::RangeType::dimension == 1);

      TraceGF(const IGV &iGV, const BulkGridFunction &bgf)
      : Base(Dune::FemPy::gridPart<IGV>(iGV), "trace_"+bgf.name(), bgf.order() ),
        bulkGridPart_(bgf.gridPart()),
        blf_(bgf)
      {}

      void bind(const typename Base::EntityType &entity)
      {
        Base::bind(entity);
        const auto intersection = Base::gridPart().grid().getMMesh().asIntersection( entity );

        if constexpr (side == Dune::Fem::IntersectionSide::in)
          blf_.bind(bulkGridPart_.convert(intersection.inside()));
        else if (intersection.neighbor())
          blf_.bind(bulkGridPart_.convert(intersection.outside()));
        else // is this the best we can do?
          DUNE_THROW( Dune::NotImplemented, "TraceFunction on boundary can no be used with outside entity" );
        sideGeometry_ = (side == IntersectionSide::in)
          ? intersection.geometryInInside().impl() : intersection.geometryInOutside().impl();
      }

      template <class Point>
      void evaluate(const Point &x, typename Base::RangeType &ret) const
      {
        // again need to transfer the x (in this case on the interface) to an x in bulk
        auto xLocal = Dune::Fem::coordinate(x);
        auto bx     = sideGeometry_.global( xLocal );
        blf_.evaluate(bx,ret);
      }

      // need to implement jacobian and hessian as their tangential components
      template <class Point>
      void jacobian(const Point &x, typename Base::JacobianRangeType &ret) const
      {
        // again need to transfer the x (in this case on the interface) to an x in bulk
        auto xLocal = Dune::Fem::coordinate(x);
        auto bx     = sideGeometry_.global( xLocal );
        typename BulkGridFunction::JacobianRangeType bulkJac;
        blf_.jacobian(bx,bulkJac);
        ret = bulkJac; // could decide to remove normal component - but perhaps don't have to?
      }

      template <class Point>
      void hessian(const Point &x, typename Base::HessianRangeType &ret) const
      {
        // again need to transfer the x (in this case on the interface) to an x in bulk
        auto xLocal = Dune::Fem::coordinate(x);
        auto bx     = sideGeometry_.global( xLocal );
        typename BulkGridFunction::HessianRangeType bulkHessian;
        blf_.hessian(bx,bulkHessian);
        ret = bulkHessian;
      }

    private:
      BLocalFunction blf_;
      const BulkGridPart& bulkGridPart_;
      SideGeometry sideGeometry_;
    };

    template< class IGV, class BulkGridFunction, Dune::Fem::IntersectionSide side >
    inline static void registerTraceGF(pybind11::module module, pybind11::class_< TraceGF<IGV, BulkGridFunction, side> > cls)
    {
      using pybind11::operator""_a;

      cls.def( pybind11::init( [] ( const IGV &iGV, const BulkGridFunction &bgf ) {
          return TraceGF<IGV, BulkGridFunction, side> ( iGV, bgf );
        } ), "iGV"_a, "bgf"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );

      cls.def_property_readonly( "scalar", [] ( TraceGF<IGV, BulkGridFunction, side> &self ) { return self.scalar; } );

      Dune::FemPy::registerGridFunction( module, cls );
    }


  }  // namespace Fem

} // namespace Dune

#endif // HAVE_DUNE_FEM

#endif // DUNE_MMESH_PYSKELETONTRACE
