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
      using Dune::Indices::_0;
      using Dune::Indices::_1;

      //! Convert intersection if gridPart is wrapped, e.g. geometryGridPart
      template< class GridPart, class Intersection, class Entity >
      const typename GridPart::IntersectionType convert( const GridPart& gridPart, const Intersection& intersection, const Entity& inside, Dune::PriorityTag<0> )
      {
        const Entity outside = gridPart.convert(intersection.outside());

        for (auto is : intersections(gridPart, inside))
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
      template <class DomainSpace, class RangeSpace, class Partition = Dune::Partitions::All>
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
      template <class DomainSpace, class RangeSpace, class Partition = Dune::Partitions::All>
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


      template<class X, class BulkDF, class InterfaceDF>
      class ParallelizedScalarProduct : public ScalarProduct<X> {
      public:
        typedef X domain_type;
        typedef typename X::field_type field_type;
        typedef typename FieldTraits<field_type>::real_type real_type;

        ParallelizedScalarProduct(const BulkDF& uh, const InterfaceDF& th)
        : u(uh), v(uh), t(th), s(th)
        {}

        virtual field_type dot (const X& x, const X& y) const
        {
          u.blockVector() = x[_0];
          t.blockVector() = x[_1];
          v.blockVector() = y[_0];
          s.blockVector() = y[_1];
          return u.scalarProductDofs(v) + t.scalarProductDofs(s); // already parallel
        }

        virtual real_type norm (const X& x) const
        {
          u.blockVector() = x[_0];
          t.blockVector() = x[_1];
          return std::sqrt( u.scalarProductDofs(u) + t.scalarProductDofs(t) ); // already parallel
        }

      private:
        mutable BulkDF u, v;
        mutable InterfaceDF t, s;
      };

      template<class M, class X, class Y, class BulkDF, class InterfaceDF>
      class ParallelizedMatrixAdapter : public AssembledLinearOperator<M,X,Y>
      {
      public:
        //! export types
        typedef M matrix_type;
        typedef X domain_type;
        typedef Y range_type;
        typedef typename X::field_type field_type;

        //! constructor: just store a reference to a matrix
        explicit ParallelizedMatrixAdapter (const M& A, const BulkDF& uh, const InterfaceDF& th)
        : _A_(stackobject_to_shared_ptr(A)), u(uh), t(th)
        {}

        //! apply operator to x:  \f$ y = A(x) \f$
        void apply (const X& x, Y& y) const override
        {
          _A_->mv(x, y);
          communicate(y);
        }

        //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
        void applyscaleadd (field_type alpha, const X& x, Y& y) const override
        {
          _A_->usmv(alpha, x, y);
          communicate(y);
        }

        //! get matrix via *
        const M& getmat () const override
        {
          return *_A_;
        }

        //! Category of the solver (see SolverCategory::Category)
        SolverCategory::Category category() const override
        {
          return SolverCategory::sequential;
        }

      private:
        void communicate(Y& y) const
        {
          u.blockVector() = y[_0];
          t.blockVector() = y[_1];
          u.communicate();
          t.communicate();
          y[_0] = u.blockVector();
          y[_1] = t.blockVector();
        }

        const std::shared_ptr<const M> _A_;
        mutable BulkDF u;
        mutable InterfaceDF t;
      };


      template< class Sch, class ISch, class Sol, class ISol >
      class Jacobian
      {
        using ThisType = Jacobian<Sch, ISch, Sol, ISol>;
      public:
        using Scheme = Sch;
        using IScheme = ISch;
        using Solution = Sol;
        using ISolution = ISol;

        using AType = typename Scheme::JacobianOperatorType;
        using BType = Dune::Fem::ISTLLinearOperator<ISolution, Solution>;
        using CType = Dune::Fem::ISTLLinearOperator<Solution, ISolution>;
        using DType = typename IScheme::JacobianOperatorType;

        // Block Matrix
        using AMatrixBlock = Dune::BCRSMatrix<typename AType::MatrixBlockType>;
        using BMatrixBlock = Dune::BCRSMatrix<typename BType::MatrixBlockType>;
        using CMatrixBlock = Dune::BCRSMatrix<typename CType::MatrixBlockType>;
        using DMatrixBlock = Dune::BCRSMatrix<typename DType::MatrixBlockType>;

        using ABRowType = Dune::MultiTypeBlockVector<AMatrixBlock, BMatrixBlock>;
        using CDRowType = Dune::MultiTypeBlockVector<CMatrixBlock, DMatrixBlock>;
        using BlockMatrix = Dune::MultiTypeBlockMatrix<ABRowType, CDRowType>;

        // Block Vector
        using AVectorBlock = typename AType::ColumnBlockVectorType;
        using DVectorBlock = typename DType::ColumnBlockVectorType;
        using BlockVector = Dune::MultiTypeBlockVector<AVectorBlock, DVectorBlock>;

        typedef typename AType::DomainSpaceType BulkSpaceType;
        typedef typename DType::DomainSpaceType InterfaceSpaceType;
        using ParameterType = Fem::NewtonParameter<typename IScheme::LinearInverseOperatorType::SolverParameterType>;

        Jacobian( const Scheme &scheme, const IScheme &ischeme, const Solution &uh, const ISolution &th,
          const double eps, const std::function<void()> &callback )
         : scheme_(scheme), ischeme_(ischeme),
           A_( "A", uh.space(), uh.space() ),
           B_( "B", th.space(), uh.space() ),
           C_( "C", uh.space(), th.space() ),
           D_( "D", th.space(), th.space() ),
           eps_(eps),
           callback_(callback),
           n_( uh.size() ),
           m_( th.size() )
        {
          x_[_0] = uh.blockVector();
          x_[_1] = th.blockVector();
        }

        const Scheme& scheme() const { return scheme_; };
        const BlockMatrix& M() const { return M_; };

        void init()
        {
          if (n_ == 0 || m_ == 0)
            return;

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

          if (n_ > 0 && m_ > 0)
          {
            B_.clear();
            C_.clear();
            assembleB(scheme_, th, uh);
            assembleC(ischeme_, uh, th);
          }
        }

        void solve( const Solution &f, const ISolution &g, Solution& u, ISolution& t )
        {
          auto setBlock = [this](const auto& block, const auto blockrow, const auto blockcol)
          {
            this->M_[blockrow][blockcol] = block.exportMatrix();
          };

          if (n_ > 0)
            setBlock(A_, _0, _0);
          if (n_ > 0 && m_ > 0)
          {
            setBlock(B_, _0, _1);
            setBlock(C_, _1, _0);
          }
          if (m_ > 0)
            setBlock(D_, _1, _1);

          b_[_0] = f.blockVector();
          b_[_1] = g.blockVector();

          auto params = std::make_shared<Fem::ISTLSolverParameter>( ParameterType( scheme_.parameter() ).linear() );
          const size_t numIterations = params->preconditionerIteration();
          const double relaxFactor = params->relaxation();

          using SolverAdapterType = Fem::ISTLSolverAdapter<-1, BlockVector>;
          using ReductionType = typename SolverAdapterType::ReductionType;
          SolverAdapterType solver( ReductionType( params ), params );
          solver.setMaxIterations( params->maxIterations() );

          ParallelizedMatrixAdapter<BlockMatrix, BlockVector, BlockVector, Solution, ISolution> linop(M_, u, t);
          ParallelizedScalarProduct<BlockVector, Solution, ISolution> scp (u, t);
          Dune::Richardson<BlockVector, BlockVector> prec(1.0);
          InverseOperatorResult res;

          solver( linop, scp, prec, b_, x_, res );

          u.blockVector() = x_[_0];
          t.blockVector() = x_[_1];
          u.communicate();
          t.communicate();
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

          for( const auto &interface : elements( gridPart, Partitions::all ) )
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

          for( const auto &element : elements( gridPart, Partitions::all ) )
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
        BlockMatrix M_;
        BlockVector x_, b_;
        std::size_t n_, m_;
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
