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

#include <dune/istl/umfpack.hh>
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

      //! Convert intersection if gridPart is wrapped, e.g. geometryGridPart
      template< class GridPart, class Intersection >
      const auto convert( const GridPart& gridPart, const Intersection& intersection, int )
      {
        const auto inside = gridPart.convert(intersection.inside());
        const auto outside = gridPart.convert(intersection.outside());

        for (auto is : intersections(gridPart, inside))
          if (is.outside() == outside)
            return is;
        DUNE_THROW(InvalidStateException, "Intersection not found!");
      }

      //! Default (trivial) convert
      template< class GridPart, class Intersection >
      const auto convert( const GridPart& gridPart, const Intersection& intersection, char )
      {
        return intersection;
      }

      //! Give preference to the non-trivial version
      template< class GridPart, class Intersection >
      const auto convert( const GridPart& gridPart, const Intersection& intersection )
      {
        return convert(gridPart, intersection, 0);
      }


      /** \class NeighborInterfaceStencil
       *  \brief Stencil contaning the entries (ien,en) for all interface entities ien
       *         and adjacent bulk entities en.
       */
      template <class DomainSpace, class RangeSpace, class Partition = Dune::Partitions::InteriorBorder>
      struct NeighborInterfaceStencil : public Fem::Stencil<DomainSpace, RangeSpace>
      {
        typedef Fem::Stencil<DomainSpace, RangeSpace> BaseType;

        typedef typename DomainSpace::GridPartType DomainGridPart;
        typedef typename RangeSpace::GridPartType RangeGridPart;

      public:
        NeighborInterfaceStencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
          : BaseType( dSpace, rSpace ),
            domainGridPart_( dSpace.gridPart() ),
            rangeGridPart_( rSpace.gridPart() )
        {}

        void setupStencil()
        {
          const auto& mmesh = domainGridPart_.grid().getMMesh();
          for( const auto& entity : elements(domainGridPart_, Partition{}) )
          {
            const auto intersection = convert(rangeGridPart_, mmesh.asIntersection( entity ));

            BaseType::fill(entity, intersection.inside());
            BaseType::fill(entity, intersection.outside());
          }
        }

      private:
        const DomainGridPart& domainGridPart_;
        const RangeGridPart& rangeGridPart_;
      };

      /** \class InterfaceNeighborStencil
       *  \brief Stencil contaning the entries (en,ien) for all entities en in the bulk
       *         and interface edges ien.
       */
      template <class DomainSpace, class RangeSpace, class Partition = Dune::Partitions::InteriorBorder>
      struct InterfaceNeighborStencil : public Fem::Stencil<DomainSpace, RangeSpace>
      {
        typedef Fem::Stencil<DomainSpace, RangeSpace> BaseType;

        typedef typename DomainSpace::GridPartType DomainGridPart;
        typedef typename RangeSpace::GridPartType RangeGridPart;

      public:
        InterfaceNeighborStencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
          : BaseType( dSpace, rSpace ),
            domainGridPart_( dSpace.gridPart() ),
            rangeGridPart_( rSpace.gridPart() )
        {}

        void setupStencil()
        {
          for( const auto& entity : elements(domainGridPart_, Partition{}) )
          {
            for( const auto& intersection : intersections(domainGridPart_, entity) )
            {
              if( intersection.neighbor() && domainGridPart_.grid().isInterface(intersection) )
              {
                auto ientity = domainGridPart_.grid().asInterfaceEntity(intersection);
                BaseType::fill(entity, ientity);
              }
            }
          }
        }

      private:
        const DomainGridPart& domainGridPart_;
        const RangeGridPart& rangeGridPart_;
      };

      template< class Sch, class ISch, class Sol, class ISol >
      class Jacobian
      {
      public:
        using Scheme = Sch;
        using IScheme = ISch;
        using Solution = Sol;
        using ISolution = ISol;

        using SparseMatrix = Dune::Fem::SparseRowMatrix<double, int>;

        using AType = typename Scheme::JacobianOperatorType;
        using BType = Dune::Fem::SparseRowLinearOperator<ISolution, Solution, SparseMatrix>;
        using CType = Dune::Fem::SparseRowLinearOperator<Solution, ISolution, SparseMatrix>;
        using DType = typename IScheme::JacobianOperatorType;

        typedef typename AType::DomainSpaceType BulkSpaceType;
        typedef typename DType::DomainSpaceType InterfaceSpaceType;

        Jacobian( const Scheme &scheme, const IScheme &ischeme, const Solution &uh, const ISolution &th,
          const double eps, const std::function<void()> &callback )
         : scheme_(scheme), ischeme_(ischeme),
           A_( "A", uh.space(), uh.space() ),
           B_( "B", th.space(), uh.space() ),
           C_( "C", uh.space(), th.space() ),
           D_( "D", th.space(), th.space() ),
           eps_(eps),
           callback_(callback)
        {}

        void init()
        {
          NeighborInterfaceStencil< InterfaceSpaceType, BulkSpaceType > stencilB( B_.domainSpace(), B_.rangeSpace() );
          stencilB.setupStencil();
          B_.clear();
          B_.reserve( stencilB );

          InterfaceNeighborStencil< BulkSpaceType, InterfaceSpaceType > stencilC( C_.domainSpace(), C_.rangeSpace() );
          stencilC.setupStencil();
          C_.clear();
          C_.reserve( stencilC );
        }

        void update( const Solution &uh, const ISolution &th )
        {
          scheme_.jacobian(uh, A_);
          ischeme_.jacobian(th, D_);

          B_.clear();
          assembleB(scheme_, th, uh);

          C_.clear();
          assembleC(ischeme_, uh, th);
        }

        void solve( const Solution &f, const ISolution &g, Solution& u, ISolution& t )
        {
          using Matrix = Dune::BCRSMatrix<double>;
          using Vector = Dune::BlockVector<double>;

          std::size_t n = f.size();
          std::size_t m = g.size();

          std::size_t nz = A_.exportMatrix().maxNzPerRow() + B_.exportMatrix().maxNzPerRow();
          Matrix M (n+m, n+m, nz, 0.1, Matrix::implicit);
          Vector x(n+m), b(n+m);

          auto addBlock = [&M](const auto& block, std::size_t row0, std::size_t col0)
          {
            auto crs = block.exportMatrix().exportCRS();
            const auto& val = std::get<0>(crs);
            const auto& col = std::get<1>(crs);
            const auto& row = std::get<2>(crs);

            for (std::size_t i = 0; i < row.size()-1; ++i)
              for (std::size_t j = row[i]; j < row[i+1]; ++j)
                if (std::abs(val[j]) > std::numeric_limits<double>::epsilon())
                  M.entry(row0 + i, col0 + col[j]) = val[j];
          };

          addBlock(A_, 0, 0);
          addBlock(B_, 0, n);
          addBlock(C_, n, 0);
          addBlock(D_, n, n);

          M.compress();

          for (std::size_t i = 0; i < n; ++i)
            b[i] = f.leakPointer()[i];
          for (std::size_t i = 0; i < m; ++i)
            b[i+n] = g.leakPointer()[i];

          Dune::UMFPack<Matrix> solver(M);
          Dune::InverseOperatorResult res;
          solver.apply(x, b, res);
          solver.free();

          for (std::size_t i = 0; i < u.size(); ++i)
            u.leakPointer()[i] = x[i];
          for (std::size_t i = 0; i < t.size(); ++i)
            t.leakPointer()[i] = x[i+n];
        }

      private:

        template< class Scheme, class DomainGridFunction, class RangeGridFunction >
        void assembleB (const Scheme &scheme, const DomainGridFunction &t, const RangeGridFunction &u)
        {
          typedef InterfaceSpaceType DomainSpaceType;
          typedef BulkSpaceType      RangeSpaceType;

          typedef Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;
          typedef typename RangeGridFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

          const auto& gridPart = t.gridPart();
          const auto& grid = gridPart.grid();
          const auto& mmesh = grid.getMMesh();

          auto dFIn = u;
          auto dFOut = u;

          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > FTmpIn( u.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > FTmpOut( u.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > dFTmpIn( u.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > dFTmpOut( u.space() );

          Dune::Fem::MutableLocalFunction< DomainGridFunction > tLocal( t );

          Dune::Fem::ConstLocalFunction< RangeGridFunction > uInside( u );
          Dune::Fem::ConstLocalFunction< RangeGridFunction > uOutside( u );

          Dune::Fem::ConstLocalFunction< RangeGridFunction > dFLocalIn( dFIn );
          Dune::Fem::ConstLocalFunction< RangeGridFunction > dFLocalOut( dFOut );

          TemporaryLocalMatrixType localMatrixIn( B_.domainSpace(), B_.rangeSpace() );
          TemporaryLocalMatrixType localMatrixOut( B_.domainSpace(), B_.rangeSpace() );

          for( const auto &interface : elements( gridPart, Partitions::interiorBorder ) )
          {
            tLocal.bind( interface );
            auto& tDof = tLocal.localDofVector();

            const auto intersection = convert(u.gridPart(), mmesh.asIntersection( interface ));

            const auto& inside = intersection.inside();
            const auto& outside = intersection.outside();

            FTmpIn.bind( inside );
            FTmpOut.bind( outside );

            dFTmpIn.bind( inside );
            dFTmpOut.bind( outside );

            uInside.bind( inside );
            uOutside.bind( outside );

            localMatrixIn.init( interface, inside );
            localMatrixOut.init( interface, outside );

            localMatrixIn.clear();
            localMatrixOut.clear();

            FTmpIn.clear();
            FTmpOut.clear();
            scheme.fullOperator().impl().addSkeletonIntegral( intersection, uInside, uOutside, FTmpIn, FTmpOut );

            for (std::size_t i = 0; i < tDof.size(); ++i)
            {
              dFIn.clear();
              dFOut.clear();

              dFTmpIn.clear();
              dFTmpOut.clear();

              dFIn.addLocalDofs( inside, FTmpIn.localDofVector() );
              dFOut.addLocalDofs( outside, FTmpOut.localDofVector() );

              dFIn *= -1.;
              dFOut *= -1.;

              tDof[ i ] += eps_;
              callback_();
              scheme.fullOperator().impl().addSkeletonIntegral( intersection, uInside, uOutside, dFTmpIn, dFTmpOut );
              tDof[ i ] -= eps_;

              dFIn.addLocalDofs( inside, dFTmpIn.localDofVector() );
              dFOut.addLocalDofs( outside, dFTmpOut.localDofVector() );

              dFIn /= eps_;
              dFOut /= eps_;

              dFLocalIn.bind( inside );
              dFLocalOut.bind( outside );

              for (std::size_t j = 0; j < dFLocalIn.localDofVector().size(); ++j)
                localMatrixIn.set(j, i, dFLocalIn[j]);

              for (std::size_t j = 0; j < dFLocalOut.localDofVector().size(); ++j)
                localMatrixOut.set(j, i, dFLocalOut[j]);
            }

            B_.addLocalMatrix( interface, inside, localMatrixIn );
            B_.addLocalMatrix( interface, outside, localMatrixOut );
          }
          B_.compress();
        }

        template< class Scheme, class DomainGridFunction, class RangeGridFunction >
        void assembleC (const Scheme &ischeme, const DomainGridFunction &u, const RangeGridFunction &t)
        {
          typedef BulkSpaceType      DomainSpaceType;
          typedef InterfaceSpaceType RangeSpaceType;

          typedef Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;
          typedef typename RangeGridFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

          const auto& gridPart = u.gridPart();
          const auto& grid = gridPart.grid();

          auto dG = t;

          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > GTmp( t.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > dGTmp( t.space() );

          Dune::Fem::MutableLocalFunction< DomainGridFunction > uLocal( u );

          Dune::Fem::ConstLocalFunction< RangeGridFunction > tInterface( t );
          Dune::Fem::ConstLocalFunction< RangeGridFunction > dGLocal( dG );

          TemporaryLocalMatrixType localMatrix( C_.domainSpace(), C_.rangeSpace() );

          for( const auto &element : elements( gridPart, Partitions::interiorBorder ) )
          {
            for( const auto& intersection : intersections( gridPart, element ) )
            {
              if( grid.isInterface( intersection ) )
              {
                const auto& interface = grid.asInterfaceEntity( intersection );

                uLocal.bind( element );
                auto& uDof = uLocal.localDofVector();

                GTmp.bind( interface );
                dGTmp.bind( interface );
                tInterface.bind( interface );

                localMatrix.init( element, interface );
                localMatrix.clear();

                GTmp.clear();
                ischeme.fullOperator().impl().addInteriorIntegral( tInterface, GTmp );

                for (std::size_t i = 0; i < uDof.size(); ++i)
                {
                  dG.clear();
                  dGTmp.clear();

                  dG.addLocalDofs( interface, GTmp.localDofVector() );
                  dG *= -1.;

                  uDof[ i ] += eps_;
                  callback_();
                  ischeme.fullOperator().impl().addInteriorIntegral( tInterface, dGTmp );
                  uDof[ i ] -= eps_;

                  dG.addLocalDofs( interface, dGTmp.localDofVector() );
                  dG /= eps_;

                  dGLocal.bind( interface );

                  for (std::size_t j = 0; j < dGLocal.localDofVector().size(); ++j)
                    localMatrix.set(j, i, dGLocal[j]);
                }

                C_.addLocalMatrix( element, interface, localMatrix );
              }
            }
          }
          C_.compress();
        }

        const Scheme& scheme_;
        const IScheme& ischeme_;
        AType A_;
        BType B_;
        CType C_;
        DType D_;
        const double eps_;
        const std::function<void()> callback_;
      };

      template< class Jacobian, class... options >
      inline static auto registerJacobian ( pybind11::handle scope, pybind11::class_< Jacobian, options... > cls )
      {
        using Solution = typename Jacobian::Solution;
        using ISolution = typename Jacobian::ISolution;

        cls.def( "init", [] ( Jacobian &self ) {
            self.init();
          },
          R"doc(
            Initialize the mixed-dimensional jacobian.
          )doc"
        );

        cls.def( "update", [] ( Jacobian &self, const Solution& uh, const ISolution& th ) {
            self.update(uh, th);
          },
          R"doc(
            Update the mixed-dimensional jacobian.
          )doc"
        );

        cls.def( "solve", [] ( Jacobian &self, const Solution& f, const ISolution& g, Solution& ux, ISolution& tx ) {
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
