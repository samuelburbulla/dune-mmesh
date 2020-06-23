#ifndef DUNE_MMESH_PYSKELETONFUNCTION
#define DUNE_MMESH_PYSKELETONFUNCTION

#if HAVE_DUNE_FEM

#include <dune/common/exceptions.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/grid/function.hh>

namespace Dune
{

  namespace Fem
  {

    template <class BulkGV, class InterfaceGridFunction>
    struct SkeletonGF
    : public BindableGridFunctionWithSpace<FemPy::GridPart<BulkGV>,
        typename InterfaceGridFunction::RangeType>
    {
      using GridPart = FemPy::GridPart<BulkGV>;
      using Base = BindableGridFunctionWithSpace<GridPart, typename InterfaceGridFunction::RangeType>;
      using ILocalFunction = ConstLocalFunction<InterfaceGridFunction>;
      using SideGeometry = typename GridPart::GridType::Intersection::LocalGeometry::Implementation;

      SkeletonGF(const BulkGV &bulkGV, const InterfaceGridFunction &igf)
      : Base(FemPy::gridPart<BulkGV>(bulkGV), "skeleton_"+igf.name(), igf.order()),
        ilf_(igf), onInterface_(false), name_(igf.name())
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
        // TODO
        DUNE_THROW( Dune::NotImplemented, "SkeletonFunction::jacobian not implemented" );
      }

      template <class Point>
      void hessian(const Point &x, typename Base::HessianRangeType &ret) const
      {
        // TODO
        DUNE_THROW( Dune::NotImplemented, "SkeletonFunction::hessian not implemented" );
      }

      unsigned int order() const { return order_; }

      std::string name() const { return "skeleton_"+name_; } // needed for output

      ILocalFunction ilf_;
      bool onInterface_;
      const std::string name_;
      int order_;
      SideGeometry sideGeometry_;
    };

  }  // namespace Fem

} // namespace Dune

#endif // HAVE_DUNE_FEM

#endif // DUNE_MMESH_PYSKELETONFUNCTION
