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

          Matrix M (n+m, n+m, 10, 0.1, Matrix::implicit);
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

          const auto& gridPart = t.gridPart();
          const auto& grid = gridPart.grid();
          const auto& mmesh = grid.getMMesh();

          NeighborInterfaceStencil< DomainSpaceType, RangeSpaceType > stencil( B.domainSpace(), B.rangeSpace() );
          B.reserve( stencil );

          auto F = u;
          scheme(u, F);
          auto dF = u;

          Dune::Fem::MutableLocalFunction< DomainGridFunction > tLocal( t );

          for( const auto &outside : elements( gridPart, Partitions::interiorBorder ) )
          {
            tLocal.bind( outside );
            auto& tDof = tLocal.localDofVector();

            double eps = 1e-8;
            for (std::size_t i = 0; i < tDof.size(); ++i)
            {
              tDof[ i ] += eps;
              scheme(u, dF); // TODO: we should do this only locally
              tDof[ i ] -= eps;

              dF -= F;
              dF /= eps;

              const auto& intersection = mmesh.asIntersection( outside );
              std::array<decltype(intersection.inside()), 2> bulkElements {{ intersection.inside(), intersection.outside() }};
              for (const auto& inside : bulkElements)
              {
                Dune::Fem::ConstLocalFunction< RangeGridFunction > dFLocal( dF );
                dFLocal.bind( inside );

                typedef Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;
                TemporaryLocalMatrixType localMatrix( B.domainSpace(), B.rangeSpace() );
                localMatrix.init(outside, inside);
// ---
                // scheme.fullOperator().addInteriorIntegral( tOutside, wInside );
                //

                //
                // Dune::Fem::ConstLocalFunction< GridFunction > uInside( u );
                //
                // addSkeletonIntegral( intersection, uInside, uOutside, wInside, wOutside );
                //
                // f2.addLocalDofs( outside, w.localDofVector() );
                //
// ---
                for (std::size_t j = 0; j < dFLocal.localDofVector().size(); ++j)
                  localMatrix.set(j, i, dFLocal[j]);

                B.addLocalMatrix( outside, inside, localMatrix );
              }
            }
          }
          B.compress();
        }

        template< class Scheme, class DomainGridFunction, class RangeGridFunction, class JacobianOperator >
        void assembleC (const Scheme &ischeme, const DomainGridFunction &u, const RangeGridFunction &t, JacobianOperator &C) const
        {
          typedef typename JacobianOperator::DomainSpaceType  DomainSpaceType;
          typedef typename JacobianOperator::RangeSpaceType   RangeSpaceType;

          const auto& gridPart = u.gridPart();
          const auto& grid = gridPart.grid();

          InterfaceNeighborStencil< DomainSpaceType, RangeSpaceType > stencil( C.domainSpace(), C.rangeSpace() );
          C.reserve( stencil );

          auto G = t;
          ischeme(t, G);
          auto dG = t;

          Dune::Fem::MutableLocalFunction< DomainGridFunction > uLocal( u );

          for( const auto &outside : elements( gridPart, Partitions::interiorBorder ) )
          {
            for( const auto& intersection : intersections( gridPart, outside ) )
            {
              if( intersection.neighbor() && grid.isInterface( intersection ) )
              {
                uLocal.bind( outside );
                auto& uDof = uLocal.localDofVector();

                double eps = 1e-8;
                for (std::size_t i = 0; i < uDof.size(); ++i)
                {
                  uDof[ i ] += eps;
                  ischeme(t, dG); // TODO: we should do this only locally
                  uDof[ i ] -= eps;

                  dG -= G;
                  dG /= eps;

                  const auto& inside = grid.asInterfaceEntity( intersection );

                  Dune::Fem::ConstLocalFunction< RangeGridFunction > dGLocal( dG );
                  dGLocal.bind( inside );

                  typedef Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > TemporaryLocalMatrixType;
                  TemporaryLocalMatrixType localMatrix( C.domainSpace(), C.rangeSpace() );
                  localMatrix.init(outside, inside);

                  for (std::size_t j = 0; j < dGLocal.localDofVector().size(); ++j)
                    localMatrix.set(j, i, dGLocal[j]);

                  C.addLocalMatrix( outside, inside, localMatrix );
                }
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
