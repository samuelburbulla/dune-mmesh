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

      /** \class NeighborInterfaceStencil
       *  \brief Stencil contaning the entries (ien,en) for all interface entities ien
       *         and adjacent bulk entities en.
       */
      template <class DomainSpace, class RangeSpace, class Partition = Dune::Partitions::InteriorBorder>
      struct NeighborInterfaceStencil : public Fem::Stencil<DomainSpace, RangeSpace>
      {
        typedef Fem::Stencil<DomainSpace, RangeSpace> BaseType;
      public:
        NeighborInterfaceStencil(const DomainSpace &dSpace, const RangeSpace &rSpace)
          : BaseType( dSpace, rSpace )
        {
          const auto& gridPart = dSpace.gridPart();
          const auto& mmesh = gridPart.grid().getMMesh();
          for( const auto& entity : elements(gridPart, Partition{}) )
          {
            const auto& intersection = mmesh.asIntersection( entity );

            BaseType::fill(entity, intersection.inside());
            BaseType::fill(entity, intersection.outside());
          }
        }
      };

      /** \class InterfaceNeighborStencil
       *  \brief Stencil contaning the entries (en,ien) for all entities en in the bulk
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

      template<class Matrix>
      void printMatrix(Matrix& M)
      {
        auto crs = M.exportCRS();
        const auto& val = std::get<0>(crs);
        const auto& col = std::get<1>(crs);
        const auto& row = std::get<2>(crs);

        for (std::size_t i = 0; i < row.size()-1; ++i)
          for (std::size_t j = row[i]; j < row[i+1]; ++j)
            std::cout << "(" << i << ", " << col[j] << ") -> " << val[j] << std::endl;
      }

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

        Jacobian( const Scheme &scheme, const IScheme &ischeme, const Solution &uh, const ISolution &th )
         : scheme_(scheme), ischeme_(ischeme),
           A_( "A", uh.space(), uh.space() ),
           B_( "B", th.space(), uh.space() ),
           C_( "C", uh.space(), th.space() ),
           D_( "D", th.space(), th.space() )
        {}

        void update( const Solution &uh, const ISolution &th )
        {
          scheme_.jacobian(uh, A_);
          ischeme_.jacobian(th, D_);

          B_.clear();
          C_.clear();
          assembleB(scheme_, th, uh, B_);
          assembleC(ischeme_, uh, th, C_);
        }

        void solve( const Solution &f, const ISolution &g, Solution& u, ISolution& t )
        {
          using Matrix = Dune::BCRSMatrix<double>;
          using Vector = Dune::BlockVector<double>;

          std::size_t n = f.size();
          std::size_t m = g.size();

          Matrix M (n+m, n+m, 10, 0.4, Matrix::implicit);
          Vector x(n+m), b(n+m);

          auto addBlock = [&M](const auto& block, std::size_t row0, std::size_t col0)
          {
            auto crs = block.exportMatrix().exportCRS();
            const auto& val = std::get<0>(crs);
            const auto& col = std::get<1>(crs);
            const auto& row = std::get<2>(crs);

            for (std::size_t i = 0; i < row.size()-1; ++i)
              for (std::size_t j = row[i]; j < row[i+1]; ++j)
                M.entry(row0 + i, col0 + col[j]) = val[j];
          };

          addBlock(A_, 0, 0);
          addBlock(B_, 0, n);
          addBlock(C_, n, 0);
          addBlock(D_, n, n);

          M.compress();

          Dune::UMFPack<Matrix> solver(M);
          Dune::InverseOperatorResult res;

          for (std::size_t i = 0; i < n; ++i)
            b[i] = f.leakPointer()[i];
          for (std::size_t i = 0; i < m; ++i)
            b[i+n] = g.leakPointer()[i];

          solver.apply(x, b, res);
          solver.free();

          for (std::size_t i = 0; i < u.size(); ++i)
            u.leakPointer()[i] = x[i];
          for (std::size_t i = 0; i < t.size(); ++i)
            t.leakPointer()[i] = x[i+n];

          /*
          std::size_t n = u.size();
          std::size_t m = t.size();

          auto& A = A_.exportMatrix();
          auto& B = B_.exportMatrix();
          auto& C = C_.exportMatrix();

          Fem::UMFPACKInverseOperator<ISolution, typename DType::MatrixType> Dinv;
          Dinv.bind(D_);

          auto ei = u;
          ei.clear();

          auto ci = t;
          auto di = t;

          Fem::SparseRowMatrix<double, int> DinvC;
          DinvC.reserve(m, n, 10);

          for (std::size_t i = 0; i < ei.size(); ++i)
          {
            ei.leakPointer()[i] = 1.;
            C.apply(ei, ci);

            Dinv(ci, di);

            for (std::size_t j = 0; j < C.rows(); ++j)
            {
              const auto& dij = di.leakPointer()[j];
              if (std::abs(dij) > std::numeric_limits<double>::epsilon())
                DinvC.set(j, i, dij);
            }
            ei.leakPointer()[i] = 0.;
          }
          DinvC.compress();

          Fem::SparseRowMatrix<double, int> BDinvC;
          BDinvC.reserve(n, n, 10);

          {
            auto crs = B.exportCRS();
            const auto& val = std::get<0>(crs);
            const auto& col = std::get<1>(crs);
            const auto& row = std::get<2>(crs);

            // BDinvC = B * DinvC
            for (std::size_t i = 0; i < row.size()-1; ++i)
              for (std::size_t j = row[i]; j < row[i+1]; ++j)
              {
                auto d = DinvC.get(col[j], i);
                if (std::abs(val[j] * d) > std::numeric_limits<double>::epsilon())
                  BDinvC.add(i, col[j], val[j] * d);
              }
            BDinvC.compress();
          }


          AType SOp("S", u.space(), u.space());

          auto& S = SOp.exportMatrix();
          S.reserve(n, n, A.maxNzPerRow() + BDinvC.maxNzPerRow());

          {
            auto crs = A.exportCRS();
            const auto& val = std::get<0>(crs);
            const auto& col = std::get<1>(crs);
            const auto& row = std::get<2>(crs);

            // S = A
            for (std::size_t i = 0; i < row.size()-1; ++i)
              for (std::size_t j = row[i]; j < row[i+1]; ++j)
                S.add(i, col[j], val[j]);
          }

          {
            auto crs = BDinvC.exportCRS();
            const auto& val = std::get<0>(crs);
            const auto& col = std::get<1>(crs);
            const auto& row = std::get<2>(crs);

            // S -= BDinvC
            for (std::size_t i = 0; i < row.size()-1; ++i)
              for (std::size_t j = row[i]; j < row[i+1]; ++j)
                S.add(i, col[j], -val[j]);
          }
          S.compress();


          Fem::UMFPACKInverseOperator<Solution, typename AType::MatrixType> Sinv;
          Sinv.bind(SOp);

          auto Dinvg = g;
          Dinv(g, Dinvg);

          auto BDinvg = u;
          B(Dinvg, BDinvg);

          auto fBDinvg = f;
          fBDinvg -= BDinvg;

          Sinv(fBDinvg, u);

          auto Cu = t;
          C.apply(u, Cu);
          auto gCu = g;
          gCu -= Cu;
          Dinv(gCu, t);

          Dinv.finalize();
          Sinv.finalize();
          */
        }

      private:

        template< class Scheme, class DomainGridFunction, class RangeGridFunction, class JacobianOperator >
        void assembleB (const Scheme &scheme, const DomainGridFunction &t, const RangeGridFunction &u, JacobianOperator &B) const
        {
          typedef typename JacobianOperator::DomainSpaceType  DomainSpaceType;
          typedef typename JacobianOperator::RangeSpaceType   RangeSpaceType;

          typedef Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;
          typedef typename RangeGridFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

          const auto& gridPart = t.gridPart();
          const auto& grid = gridPart.grid();
          const auto& mmesh = grid.getMMesh();

          NeighborInterfaceStencil< DomainSpaceType, RangeSpaceType > stencil( B.domainSpace(), B.rangeSpace() );
          B.reserve( stencil );

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

          TemporaryLocalMatrixType localMatrixIn( B.domainSpace(), B.rangeSpace() );
          TemporaryLocalMatrixType localMatrixOut( B.domainSpace(), B.rangeSpace() );

          for( const auto &interface : elements( gridPart, Partitions::interiorBorder ) )
          {
            tLocal.bind( interface );
            auto& tDof = tLocal.localDofVector();

            const auto& intersection = mmesh.asIntersection( interface );

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

            static const double eps = 1e-8;
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

              tDof[ i ] += eps;
              scheme.fullOperator().impl().addSkeletonIntegral( intersection, uInside, uOutside, dFTmpIn, dFTmpOut );
              tDof[ i ] -= eps;

              dFIn.addLocalDofs( inside, dFTmpIn.localDofVector() );
              dFOut.addLocalDofs( outside, dFTmpOut.localDofVector() );

              dFIn /= eps;
              dFOut /= eps;

              dFLocalIn.bind( inside );
              dFLocalOut.bind( outside );

              for (std::size_t j = 0; j < dFLocalIn.localDofVector().size(); ++j)
                localMatrixIn.set(j, i, dFLocalIn[j]);

              for (std::size_t j = 0; j < dFLocalOut.localDofVector().size(); ++j)
                localMatrixOut.set(j, i, dFLocalOut[j]);
            }

            B.addLocalMatrix( interface, inside, localMatrixIn );
            B.addLocalMatrix( interface, outside, localMatrixOut );
          }
          B.compress();
        }

        template< class Scheme, class DomainGridFunction, class RangeGridFunction, class JacobianOperator >
        void assembleC (const Scheme &ischeme, const DomainGridFunction &u, const RangeGridFunction &t, JacobianOperator &C) const
        {
          typedef typename JacobianOperator::DomainSpaceType  DomainSpaceType;
          typedef typename JacobianOperator::RangeSpaceType   RangeSpaceType;

          typedef Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;
          typedef typename RangeGridFunction::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

          const auto& gridPart = u.gridPart();
          const auto& grid = gridPart.grid();

          InterfaceNeighborStencil< DomainSpaceType, RangeSpaceType > stencil( C.domainSpace(), C.rangeSpace() );
          C.reserve( stencil );

          auto dG = t;

          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > GTmp( t.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > dGTmp( t.space() );

          Dune::Fem::MutableLocalFunction< DomainGridFunction > uLocal( u );

          Dune::Fem::ConstLocalFunction< RangeGridFunction > tInterface( t );
          Dune::Fem::ConstLocalFunction< RangeGridFunction > dGLocal( dG );

          TemporaryLocalMatrixType localMatrix( C.domainSpace(), C.rangeSpace() );

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

                static const double eps = 1e-8;
                for (std::size_t i = 0; i < uDof.size(); ++i)
                {
                  dG.clear();
                  dGTmp.clear();

                  dG.addLocalDofs( interface, GTmp.localDofVector() );
                  dG *= -1.;

                  uDof[ i ] += eps;
                  ischeme.fullOperator().impl().addInteriorIntegral( tInterface, dGTmp );
                  uDof[ i ] -= eps;

                  dG.addLocalDofs( interface, dGTmp.localDofVector() );
                  dG /= eps;

                  dGLocal.bind( interface );

                  for (std::size_t j = 0; j < dGLocal.localDofVector().size(); ++j)
                    localMatrix.set(j, i, dGLocal[j]);
                }

                C.addLocalMatrix( element, interface, localMatrix );
              }
            }
          }
          C.compress();
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
