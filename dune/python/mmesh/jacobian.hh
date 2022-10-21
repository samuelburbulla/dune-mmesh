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
      template< class GridPart, class Intersection, class Entity >
      const typename GridPart::IntersectionType convert( const GridPart& gridPart, const Intersection& intersection, const Entity& inside, Dune::PriorityTag<0> )
      {
        const Entity outside = gridPart.convert(intersection.outside());

        for (auto is : intersections(gridPart, inside))
          if (is.neighbor())
            if (is.outside() == outside)
              return is;
        DUNE_THROW(InvalidStateException, "Intersection not found!");
      }

      //! Default (trivial) convert
      template< class GridPart, class Intersection, class Entity >
      std::enable_if_t<std::is_convertible_v<Intersection, typename GridPart::IntersectionType>, const typename GridPart::IntersectionType>
      convert( const GridPart& gridPart, const Intersection& intersection, const Entity& inside, Dune::PriorityTag<1> )
      {
        return intersection;
      }

      //! Give preference to the non-trivial version
      template< class GridPart, class Intersection, class Entity >
      const typename GridPart::IntersectionType convert( const GridPart& gridPart, const Intersection& intersection, const Entity& inside )
      {
        return convert(gridPart, intersection, inside, Dune::PriorityTag<1>{});
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

        void setupStencil() const
        {
          const auto& mmesh = domainGridPart_.grid().getMMesh();
          for( const auto& entity : elements(domainGridPart_, Partition{}) )
          {
            const auto mmeshIntersection = mmesh.asIntersection( entity );
            const auto inside = rangeGridPart_.convert( mmeshIntersection.inside() );
            const auto intersection = convert( rangeGridPart_, mmeshIntersection, inside );

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

        void setupStencil() const
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

        using UMFPackMatrix = Dune::BCRSMatrix<double>;
        using UMFPackVector = Dune::BlockVector<double>;

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
          B_.reserve( stencilB );

          InterfaceNeighborStencil< BulkSpaceType, InterfaceSpaceType > stencilC( C_.domainSpace(), C_.rangeSpace() );
          stencilC.setupStencil();
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
          if (M_.buildStage() != UMFPackMatrix::built)
            setup();

          std::size_t n = B_.exportMatrix().rows();
          std::size_t m = B_.exportMatrix().cols();

          auto setBlock = [this](const auto& block, std::size_t row0, std::size_t col0)
          {
            auto crs = block.exportMatrix().exportCRS();
            const auto& val = std::get<0>(crs);
            const auto& col = std::get<1>(crs);
            const auto& row = std::get<2>(crs);

            for (std::size_t i = 0; i < row.size()-1; ++i)
              for (std::size_t j = row[i]; j < row[i+1]; ++j)
                  this->M_[row0 + i][col0 + col[j]] = val[j];
          };

          setBlock(A_, 0, 0);
          setBlock(B_, 0, n);
          setBlock(C_, n, 0);
          setBlock(D_, n, n);

          for (std::size_t i = 0; i < n; ++i)
            b_[i] = f.leakPointer()[i];
          for (std::size_t i = 0; i < m; ++i)
            b_[i+n] = g.leakPointer()[i];

          Dune::UMFPack<UMFPackMatrix> solver(M_);
          Dune::InverseOperatorResult res;
          solver.apply(x_, b_, res);
          solver.free();

          for (std::size_t i = 0; i < u.size(); ++i)
            u.leakPointer()[i] = x_[i];
          for (std::size_t i = 0; i < t.size(); ++i)
            t.leakPointer()[i] = x_[i+n];
        }

      private:

        // Setup UMFPackMatrix
        void setup()
        {
          std::size_t n = B_.exportMatrix().rows();
          std::size_t m = B_.exportMatrix().cols();
          x_.resize(n+m);
          b_.resize(n+m);

          std::size_t nz = std::max( A_.exportMatrix().maxNzPerRow() + B_.exportMatrix().maxNzPerRow(),
                                     C_.exportMatrix().maxNzPerRow() + D_.exportMatrix().maxNzPerRow() );
          M_.setBuildMode(UMFPackMatrix::implicit);
          M_.setImplicitBuildModeParameters(nz, 0.1);
          M_.setSize(n+m, n+m);

          auto addBlock = [this](const auto& block, std::size_t row0, std::size_t col0)
          {
            auto crs = block.exportMatrix().exportCRS();
            const auto& val = std::get<0>(crs);
            const auto& col = std::get<1>(crs);
            const auto& row = std::get<2>(crs);

            for (std::size_t i = 0; i < row.size()-1; ++i)
              for (std::size_t j = row[i]; j < row[i+1]; ++j)
                  this->M_.entry(row0 + i, col0 + col[j]) = val[j];
          };

          addBlock(A_, 0, 0);
          addBlock(B_, 0, n);
          addBlock(C_, n, 0);
          addBlock(D_, n, n);

          M_.compress();
        }


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

          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > FmTmpIn( u.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > FmTmpOut( u.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > FpTmpIn( u.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > FpTmpOut( u.space() );

          Dune::Fem::MutableLocalFunction< DomainGridFunction > tLocal( t );

          Dune::Fem::ConstLocalFunction< RangeGridFunction > uInside( u );
          Dune::Fem::ConstLocalFunction< RangeGridFunction > uOutside( u );

          Dune::Fem::MutableLocalFunction< RangeGridFunction > dFLocalIn( dFIn );
          Dune::Fem::MutableLocalFunction< RangeGridFunction > dFLocalOut( dFOut );

          TemporaryLocalMatrixType localMatrixIn( B_.domainSpace(), B_.rangeSpace() );
          TemporaryLocalMatrixType localMatrixOut( B_.domainSpace(), B_.rangeSpace() );

          for( const auto &interface : elements( gridPart, Partitions::interiorBorder ) )
          {
            tLocal.bind( interface );
            auto& tDof = tLocal.localDofVector();

            const auto mmeshIntersection = mmesh.asIntersection( interface );
            const auto inside = u.gridPart().convert( mmeshIntersection.inside() );
            const auto intersection = convert( u.gridPart(), mmeshIntersection, inside );

            const auto& outside = intersection.outside();

            FmTmpIn.bind( inside );
            FmTmpOut.bind( outside );
            FpTmpIn.bind( inside );
            FpTmpOut.bind( outside );

            uInside.bind( inside );
            uOutside.bind( outside );

            dFLocalIn.bind( inside );
            dFLocalOut.bind( outside );

            localMatrixIn.init( interface, inside );
            localMatrixOut.init( interface, outside );

            localMatrixIn.clear();
            localMatrixOut.clear();

            for (std::size_t i = 0; i < tDof.size(); ++i)
            {
              dFIn.clear();
              dFOut.clear();

              FmTmpIn.clear();
              FmTmpOut.clear();
              FpTmpIn.clear();
              FpTmpOut.clear();

              double h = std::max(tDof[ i ] * eps_, eps_);
              tDof[ i ] -= h;
              callback_();
              scheme.fullOperator().impl().addSkeletonIntegral( intersection, uInside, uOutside, FmTmpIn, FmTmpOut );
              tDof[ i ] += h;

              dFIn.addLocalDofs( inside, FmTmpIn.localDofVector() );
              dFOut.addLocalDofs( outside, FmTmpOut.localDofVector() );

              dFIn *= -1.;
              dFOut *= -1.;

              tDof[ i ] += h;
              callback_();
              scheme.fullOperator().impl().addSkeletonIntegral( intersection, uInside, uOutside, FpTmpIn, FpTmpOut );
              tDof[ i ] -= h;

              dFIn.addLocalDofs( inside, FpTmpIn.localDofVector() );
              dFOut.addLocalDofs( outside, FpTmpOut.localDofVector() );

              dFIn /= 2 * h;
              dFOut /= 2 * h;

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

          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > GmTmp( t.space() );
          Dune::Fem::TemporaryLocalFunction< DiscreteFunctionSpaceType > GpTmp( t.space() );

          Dune::Fem::MutableLocalFunction< DomainGridFunction > uLocal( u );

          Dune::Fem::ConstLocalFunction< RangeGridFunction > tInterface( t );
          Dune::Fem::MutableLocalFunction< RangeGridFunction > dGLocal( dG );

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

                GmTmp.bind( interface );
                GpTmp.bind( interface );
                tInterface.bind( interface );

                dGLocal.bind( interface );
                localMatrix.init( element, interface );
                localMatrix.clear();

                for (std::size_t i = 0; i < uDof.size(); ++i)
                {
                  dG.clear();

                  GmTmp.clear();
                  GpTmp.clear();

                  double h = std::max(uDof[ i ] * eps_, eps_);
                  uDof[ i ] -= h;
                  callback_();
                  ischeme.fullOperator().impl().addInteriorIntegral( tInterface, GmTmp );
                  uDof[ i ] += h;

                  dG.addLocalDofs( interface, GmTmp.localDofVector() );
                  dG *= -1.;

                  uDof[ i ] += h;
                  callback_();
                  ischeme.fullOperator().impl().addInteriorIntegral( tInterface, GpTmp );
                  uDof[ i ] -= h;

                  dG.addLocalDofs( interface, GpTmp.localDofVector() );
                  dG /= 2 * h;

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
        UMFPackMatrix M_;
        UMFPackVector x_, b_;
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
