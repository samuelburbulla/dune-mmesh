// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_MMESH_JACOBIAN_HH
#define DUNE_PYTHON_MMESH_JACOBIAN_HH

#include <array>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <type_traits>

#include <dune/fem/operator/linear/spoperator.hh>

#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dune
{

  namespace Python
  {

    namespace MMesh
    {
      /** \class InterfaceNeighborStencil
       *  \brief Stencil contaning the entries (en,ien) for all entities en in the space
       *         and interface edges ien.
       */
      template <class DomainSpace, class RangeSpace, class Partition = Dune::Partitions::InteriorBorder>
      struct InterfaceNeighborStencil : public Fem::Stencil<DomainSpace, RangeSpace>
      {
        typedef Fem::Stencil<DomainSpace, RangeSpace> BaseType;
      public:
        InterfaceNeighborStencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
          : BaseType( dSpace, rSpace )
        {
          const auto& gridPart = dSpace.gridPart();
          for( const auto& entity : elements(gridPart, Partition{}) )
          {
            for( const auto& intersection : intersections(gridPart, entity) )
            {
              if( intersection.neighbor() && gridPart.grid().isInterface(intersection) )
              {
                auto ientity = gridPart.grid().asInterfaceEntity(intersection);
                BaseType::fill(entity, ientity);
              }
            }
          }
        }
      };

      template< class Sch, class ISch, class Sol, class ISol >
      class Jacobian
       : public Fem::Impl::GalerkinOperator< typename Sch::ModelType >
      {
      public:
        using Base = Fem::Impl::GalerkinOperator< typename Sch::ModelType >;
        using typename Base::DomainValueVectorType;

        using Scheme = Sch;
        using IScheme = ISch;
        using Solution = Sol;
        using ISolution = ISol;

        using AType = typename Scheme::JacobianOperatorType;
        using BType = Dune::Fem::SparseRowLinearOperator<ISol, Sol, Dune::Fem::SparseRowMatrix<double, int>>;
        using CType = Dune::Fem::SparseRowLinearOperator<Sol, ISol, Dune::Fem::SparseRowMatrix<double, int>>;
        using DType = typename IScheme::JacobianOperatorType;

        Jacobian( const Scheme &scheme, const IScheme &ischeme, const Solution &uh, const ISolution &th )
         : Base(scheme.gridPart(), scheme.model()),
           scheme_(scheme), ischeme_(ischeme),
           A_( "A", uh.space(), uh.space() ),
           B_( "B", th.space(), uh.space() ),
           C_( "C", uh.space(), th.space() ),
           D_( "D", th.space(), th.space() )
        {}

        void update( const Solution &uh, const ISolution &th )
        {
          scheme_.jacobian(uh, A_);
          ischeme_.jacobian(th, D_);

          // assemble(scheme_, th, uh, B_);
          assemble(ischeme_, uh, th, C_);
        }

        void solve( const Solution &f, const ISolution &g, const Solution& ux, const ISolution& tx ) const
        {
          // use Schur complement and call umfpack
        }

      private:

        template< class Scheme, class GridPart, class Intersection, class U, class T, class J, class K >
        void addLinearizedSkeletonIntegral ( const Scheme &scheme, const GridPart &gridPart, const Intersection &intersection, const U &uIn, const T &uOut,
                                             DomainValueVectorType &phiIn, DomainValueVectorType &phiOut, J &jInIn, K &jOutIn ) const
        {
          const auto &domainBasisIn = jInIn.domainBasisFunctionSet();
          const auto &domainBasisOut = jOutIn.rangeBasisFunctionSet();

          const auto &rangeBasisIn = jInIn.rangeBasisFunctionSet();

          const int order = std::max( std::max( this->maxOrder(uIn), this->maxOrder(uOut) ), this->maxOrder( domainBasisIn, domainBasisOut, rangeBasisIn ) );

          const auto geometry = intersection.geometry();
          typedef typename Base::template QuadratureSelector< typename J::RangeSpaceType >::SurfaceQuadratureType SurfaceQuadratureType;
          const SurfaceQuadratureType quadrature( gridPart, intersection, Base::surfaceQuadratureOrder(order), SurfaceQuadratureType::INSIDE );

          typedef Fem::CachingQuadrature< typename Scheme::GridPartType, 0, Fem::Capabilities::DefaultQuadrature< typename J::RangeSpaceType >::template DefaultQuadratureTraits> InteriorQuadratureType;
          const InteriorQuadratureType iquadrature( uOut.entity(), Base::interiorQuadratureOrder(order) );

          for( std::size_t qp = 0, nop = quadrature.nop(); qp != nop; ++qp )
          {
            const auto weight = quadrature.weight( qp ) * geometry.integrationElement( quadrature.localPoint( qp ) );

            const auto qpIn = quadrature[ qp ];
            const auto qpOut = iquadrature[ qp ];

            this->values( domainBasisIn, qpIn, phiIn );
            this->values( domainBasisOut, qpOut, phiOut );

            typename Base::RangeValueType rangeValue;
            this->value( uOut, qpOut, rangeValue );

            auto integrand = linearizedSkeleton( scheme.model(), qpIn, this->domainValue( uIn, qpIn ), qpOut, rangeValue );
            for( std::size_t col = 0, cols = domainBasisIn.size(); col < cols; ++col )
            {
              Fem::LocalMatrixColumn< J > jInInCol( jInIn, col );
              auto intPhi = integrand.first( this->value( phiIn, col ) );

              Hybrid::forEach( typename Base::RangeValueIndices(), [ &qpIn, &jInInCol, &intPhi, weight ] ( auto i ) {
                  std::get< i >( intPhi.first ) *= weight;
                  jInInCol.axpy( qpIn, std::get< i >( intPhi.first ) );
                } );
            }
            for( std::size_t col = 0, cols = domainBasisOut.size(); col < cols; ++col )
            {
              Fem::LocalMatrixColumn< J > jOutInCol( jOutIn, col );
              auto intPhi = integrand.second( this->value( phiOut, col ) );

              Hybrid::forEach( typename Base::RangeValueIndices(), [ &qpIn, &jOutInCol, &intPhi, weight ] ( auto i ) {
                  std::get< i >( intPhi.first ) *= weight;
                  jOutInCol.axpy( qpIn, std::get< i >( intPhi.first ) );
                } );
            }
          }
        }

        template< class Scheme, class DomainGridFunction, class RangeGridFunction, class JacobianOperator >
        void assemble (const Scheme &scheme, const DomainGridFunction &u, const RangeGridFunction &t, JacobianOperator &jOp) const
        {
          typedef typename JacobianOperator::DomainSpaceType  DomainSpaceType;
          typedef typename JacobianOperator::RangeSpaceType   RangeSpaceType;

          typedef Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;

          // select correct default quadrature orders
          auto defaultInteriorOrder_ = [] (const int order) { return Fem::Capabilities::DefaultQuadrature< RangeSpaceType >::volumeOrder(order); };
          auto defaultSurfaceOrder_  = [] (const int order) { return Fem::Capabilities::DefaultQuadrature< RangeSpaceType >::surfaceOrder(order); };

          InterfaceNeighborStencil< DomainSpaceType, RangeSpaceType > stencil( jOp.domainSpace(), jOp.rangeSpace() );
          jOp.reserve( stencil );
          jOp.clear();

          const std::size_t maxNumLocalDofs = jOp.domainSpace().blockMapper().maxNumDofs() * jOp.domainSpace().localBlockSize;
          auto phiIn = Base::makeDomainValueVector( maxNumLocalDofs );
          auto phiOut = Base::makeDomainValueVector( maxNumLocalDofs );

          TemporaryLocalMatrixType jOpInIn( jOp.domainSpace(), jOp.rangeSpace() ), jOpOutIn( jOp.domainSpace(), jOp.rangeSpace() );
          Dune::Fem::ConstLocalFunction< DomainGridFunction > uIn( u );

          const auto &gridPart = u.gridPart();
          const auto &grid = gridPart.grid();
          for( const auto &inside : elements( gridPart, Partitions::interiorBorder ) )
          {
            uIn.bind( inside );

            for( const auto &intersection : intersections( gridPart, inside ) )
            {
              if( intersection.neighbor() && grid.isInterface( intersection ) )
              {
                const auto &outside = grid.asInterfaceEntity( intersection );

                jOpOutIn.init( inside, outside );
                jOpOutIn.clear();

                Dune::Fem::ConstLocalFunction< RangeGridFunction > tOut( t );
                tOut.bind( outside );

                addLinearizedSkeletonIntegral( scheme, gridPart, intersection, uIn, tOut, phiIn, phiOut, jOpInIn, jOpOutIn );
                jOp.addLocalMatrix( inside, outside, jOpOutIn );
              }
            }
          }
        }

        const Scheme& scheme_;
        const IScheme& ischeme_;
        AType A_;
        BType B_;
        CType C_;
        DType D_;
      };


      template< class Jacobian, class... options >
      inline static auto registerJacobian ( pybind11::handle scope, pybind11::class_< Jacobian, options... > cls )
      {
        using Solution = typename Jacobian::Solution;
        using ISolution = typename Jacobian::ISolution;

        cls.def( "update", [] ( Jacobian &self, const Solution& uh, const ISolution& th ) {
            self.update(uh, th);
          },
          R"doc(
            Update the mixed-dimensional jacobian.
          )doc"
        );

        cls.def( "solve", [] ( Jacobian &self, const Solution& f, const ISolution& g, const Solution& ux, const ISolution& tx ) {
            self.solve(f, g, ux, tx);
          },
          R"doc(
            Solve the mixed-dimensional jacobian.
          )doc"
        );
      }

    } // namespace MMesh

  } // namespace Python

} // namespace Dune

#endif
