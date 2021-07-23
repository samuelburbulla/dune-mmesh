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
      struct InterfaceNeighborStencil : public Fem::Stencil<DomainSpace,RangeSpace>
      {
        typedef Fem::Stencil<DomainSpace,RangeSpace> BaseType;
      public:
        typedef Partition                                 PartitionType;
        typedef typename BaseType::DomainEntityType       DomainEntityType;
        typedef typename BaseType::RangeEntityType        RangeEntityType;
        typedef typename BaseType::DomainGlobalKeyType    DomainGlobalKeyType;
        typedef typename BaseType::RangeGlobalKeyType     RangeGlobalKeyType;
        typedef typename BaseType::LocalStencilType       LocalStencilType;
        typedef typename BaseType::GlobalStencilType      GlobalStencilType;

        InterfaceNeighborStencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
          : BaseType( dSpace, rSpace )
        {
          const auto& gridPart = rSpace.gridPart();
          const auto& grid = gridPart.grid();
          for (const auto & entity: elements(gridPart, PartitionType{}) )
          {
            for (const auto & intersection: intersections(gridPart, entity) )
            {
              if( intersection.neighbor() && grid.isInterface(intersection) )
              {
                auto ientity = grid.asInterfaceEntity(intersection);
                BaseType::fill(ientity, entity);
              }
            }
          }
        }
      };

      template< class Sch, class ISch, class Sol, class ISol >
      class Jacobian
      {
      public:
        using Scheme = Sch;
        using IScheme = ISch;
        using Solution = Sol;
        using ISolution = ISol;

        using AType = typename Scheme::JacobianOperatorType;
        using BType = Dune::Fem::SparseRowLinearOperator<ISol, Sol, Dune::Fem::SparseRowMatrix<double, int>>;
        using CType = Dune::Fem::SparseRowLinearOperator<Sol, ISol, Dune::Fem::SparseRowMatrix<double, int>>;
        using DType = typename IScheme::JacobianOperatorType;

        Jacobian( const Scheme &scheme, const IScheme &ischeme, const Solution &uh, const ISolution &th )
         : scheme_(scheme), ischeme_(ischeme),
           A_( "A", uh.space(), uh.space() ),
           // B_( "B", th.space(), uh.space() ),
           B_( "B", th.space(), uh.space() ),
           C_( "C", uh.space(), th.space() ),
           D_( "D", th.space(), th.space() )
        {}

        void update( const Solution &uh, const ISolution &th )
        {
          scheme_.jacobian(uh, A_);
          ischeme_.jacobian(th, D_);

          assemble(scheme_, uh, th, B_);
          // assemble(ischeme_, th, uh, C_);
        }

        void solve( const Solution &f, const ISolution &g, const Solution& ux, const ISolution& tx ) const
        {
          // use Schur complement and call umfpack
        }

      private:
        typedef typename Scheme::ModelType::DomainValueType DomainValueType;
        typedef typename Scheme::ModelType::RangeValueType RangeValueType;
        typedef std::make_index_sequence< std::tuple_size< DomainValueType >::value > DomainValueIndices;
        typedef std::make_index_sequence< std::tuple_size< RangeValueType >::value > RangeValueIndices;

        template< std::size_t... i >
        static auto makeDomainValueVector ( std::size_t maxNumLocalDofs, std::index_sequence< i... > )
        {
          return std::make_tuple( std::vector< std::tuple_element_t< i, DomainValueType > >( maxNumLocalDofs )... );
        }

        static auto makeDomainValueVector ( std::size_t maxNumLocalDofs )
        {
          return makeDomainValueVector( maxNumLocalDofs, DomainValueIndices() );
        }

        template< class Scheme, class GridFunction, class InterfaceGridFunction, class JacobianOperator >
        void assemble (const Scheme &scheme, const GridFunction &u, const InterfaceGridFunction &t, JacobianOperator &jOp) const
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
          auto phiIn = makeDomainValueVector( maxNumLocalDofs );
          auto phiOut = makeDomainValueVector( maxNumLocalDofs );

          TemporaryLocalMatrixType jOpInIn( jOp.domainSpace(), jOp.rangeSpace() ), jOpOutIn( jOp.domainSpace(), jOp.rangeSpace() );
          Dune::Fem::ConstLocalFunction< GridFunction > uIn( u );

          const auto &grid = scheme.gridPart().grid();
          const auto &indexSet = scheme.gridPart().indexSet();
          for( const auto &inside : elements( scheme.gridPart(), Partitions::interiorBorder ) )
          {
            uIn.bind( inside );

            // jOpInIn.init( inside, inside );
            // jOpInIn.clear();

            for( const auto &intersection : intersections( scheme.gridPart(), inside ) )
            {
              if( intersection.neighbor() && grid.isInterface( intersection ) )
              {
                const auto &outside = grid.asInterfaceEntity( intersection );

                jOpOutIn.init( outside, inside );
                jOpOutIn.clear();

                Dune::Fem::ConstLocalFunction< InterfaceGridFunction > uOut( t );
                uOut.bind( outside );

                scheme.addLinearizedSkeletonIntegral( intersection, uIn, uOut, phiIn, phiOut, jOpInIn, jOpOutIn );
                jOp.addLocalMatrix( outside, inside, jOpOutIn );
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
