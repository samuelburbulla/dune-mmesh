// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/misc/l2norm.hh>

#include <dune/mmesh/mmesh.hh>
#include <dune/alugrid/grid.hh>

using namespace Dune;
using namespace Fem;

static constexpr int dim = GRIDDIM;
using BulkGridType = MovingMesh<dim>;
#ifdef INTERFACE
  using GridType = typename MovingMesh<dim>::InterfaceGrid;
#else
  using GridType = BulkGridType;
#endif
static constexpr bool verbose = false;

// Custom data handle
// ==================

template <class Data, class IndexSet>
struct DataHandle
: CommDataHandleIF< DataHandle<Data, IndexSet>, Data >
{
  using DataType = typename Data::value_type;

  DataHandle (Data& data, Data& dataCheck, const IndexSet& indexSet, int codim = 0)
  : data_(data), dataCheck_(dataCheck), indexSet_(indexSet), mycodim_(codim)
  {}

  bool contains( int dimension, int codim ) const
  {
    return codim == mycodim_;
  }

  template <class Entity>
  std::size_t size( const Entity &entity ) const
  {
    return 1;
  }

  template <class Buffer, class Entity>
  void gather( Buffer& buffer, const Entity& entity ) const
  {
    DataType d = data_[ indexSet_.index(entity) ];
    buffer.write( d );
  }

  template <class Buffer, class Entity>
  void scatter( Buffer& buffer, const Entity& entity, std::size_t size )
  {
    DataType d;
    buffer.read( d );
    auto i = indexSet_.index(entity);
    data_[ i ] = d;

    if(data_[i] != dataCheck_[i])
      DUNE_THROW(InvalidStateException, "Element at (" << entity.geometry().center() << ") is wrong: " << data_[i] << " != " << dataCheck_[i]);
  }

  static bool fixedSize( int dim, int codim )
  {
    return false;
  }

  Data& data_, dataCheck_;
  const IndexSet& indexSet_;
  const int mycodim_;
};



// Types and classes for testing communication of discrete functions
// =================================================================

static constexpr int polOrd = 1;
using GridPartType = AdaptiveLeafGridPart< GridType >;
using SpaceType = FunctionSpace<double, double, dim, 1>;
using DiscreteFunctionSpaceType = LagrangeDiscontinuousGalerkinSpace<SpaceType, GridPartType, polOrd>;
using DiscreteFunctionType = AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >;

struct ExactSolution
: public Fem::Function< SpaceType, ExactSolution >
{
  typedef SpaceType::RangeType RangeType;
  typedef SpaceType::RangeFieldType RangeFieldType;
  typedef SpaceType::DomainType DomainType;

  //! f(x,y) = x*(1-x)*y*(1-y)
  void evaluate ( const DomainType &x, RangeType &ret ) const
  {
    ret = 2.; // maximum of function is 2
    for( int i = 0; i < DomainType::dimension; ++i )
      ret *= sin( x[ i ]*(1.0 -x[ i ])*4.);
  }
};

struct ProjectionAllPartitionNoComm
{
  template <class FunctionType, class DiscreteFunctionType, class Partition = Partitions::All>
  static void project (const FunctionType &f, DiscreteFunctionType &u)
  {
    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;

    const auto& space = u.space();
    u.clear();

    TemporaryLocalFunction< DiscreteFunctionSpaceType > dfLocal(space);
    RangeType ret, phi;


    for(const auto& entity : elements(space.gridPart(), Partition{}))
    {
      const auto& geo = entity.geometry();
      dfLocal.bind(entity);
      dfLocal.clear();

      CachingQuadrature<GridPartType, 0> quad(entity, 2 * space.order() + 2);
      for(size_t qP = 0; qP < quad.nop(); ++qP)
      {
        auto x = geo.global(quad.point(qP));
        f.evaluate(x, ret);
        ret *= quad.weight(qP);
        dfLocal.axpy(quad[qP], ret);
      }

      u.addLocalDofs(entity, dfLocal);
      dfLocal.unbind();
    }
  }
};


int main(int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);

  GridPtr< BulkGridType > gridptr_ ((dim == 2) ? "grids/junction2d.msh" : "grids/plane3d.msh");
  BulkGridType& bulkGrid = *gridptr_;

#ifdef INTERFACE
  GridType& grid = bulkGrid.interfaceGrid();
#else
  GridType& grid = bulkGrid;
#endif

  // Communicate data with custom data handle
  // ========================================

  // Initialize data
  using IdType = typename GridType::GlobalIdSet::IdType;
  using DataType = std::vector<IdType>;
  DataType data (grid.size(0)), dataCheck (grid.size(0));
  for (const auto& e : elements(grid.leafGridView(), Partitions::all))
  {
    auto idx = grid.leafIndexSet().index(e);
    auto id = grid.globalIdSet().id(e);
    dataCheck[idx] = id;

    if (e.partitionType() == InteriorEntity)
      data[idx] = id;
  }

  static constexpr int dimgrid = GridType::dimension;
  DataType dataVertex(grid.size(dimgrid)), dataVertexCheck(grid.size(dimgrid));
  for (const auto& v : vertices(grid.leafGridView(), Partitions::all))
  {
    auto idx = grid.leafIndexSet().index(v);
    auto id = grid.globalIdSet().id(v);
    dataVertexCheck[idx] = id;

    if (v.partitionType() == InteriorEntity || v.partitionType() == BorderEntity)
      dataVertex[idx] = id;
  }

  // Print
  auto printData = [&](bool before = true)
  {
    if (!verbose) return;
    std::cout << "Rank " << grid.comm().rank() << " (" << (before ? "before" : "after ") << "): ";
    for (const auto d : data)
      std::cout << d << "|  ";
    for (const auto d : dataVertex)
      std::cout << d << "|  ";
    std::cout << std::flush << std::endl;
  };

  // Communicate
  using DataHandleType = DataHandle< DataType, typename GridType::LeafIndexSet >;
  DataHandleType dataHandle (data, dataCheck, grid.leafIndexSet());
  DataHandleType dataHandleVertex (dataVertex, dataVertexCheck, grid.leafIndexSet(), dimgrid);

  // Print
  printData();

  Dune::Timer timer; timer.start();
  grid.communicate( dataHandle, InteriorBorder_All_Interface, ForwardCommunication );
  grid.communicate( dataHandleVertex, InteriorBorder_All_Interface, ForwardCommunication );
  auto maxT = grid.comm().max( timer.elapsed() );
  if (grid.comm().rank() == 0)
    std::cout << "Comm took " << maxT << std::endl;

  // Print
  printData(false);

  // Communicate discrete function
  // =============================
  MPIManager::initialize( argc, argv );

  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType space( gridPart );
  DiscreteFunctionType solution( "solution", space );
  ExactSolution f;

  // create L2 Norm, no communication
  L2Norm< GridPartType > l2norm( solution.space().gridPart(), 2*solution.space().order(), false );
  solution.clear();

  ProjectionAllPartitionNoComm::project( f, solution );

  // compute l2 error on all elements
  double old_error = l2norm.distance( f, solution );
  if (verbose)
    std::cout << "P[" << grid.comm().rank() << "]  start comm: " << old_error << std::endl;

  solution.clear();
  ProjectionAllPartitionNoComm::project<ExactSolution, DiscreteFunctionType, Partitions::InteriorBorder>( f, solution );

  // communicate
  timer.reset();
  timer.start();
  space.communicate(solution, DFCommunicationOperation::Add());
  timer.stop();

  // calculate l2 error again on all elements
  double error = l2norm.distance( f, solution );
  if (verbose)
    std::cout << "P[" << grid.comm().rank() << "]  done comm: " << error << std::endl;

  if( std::abs( old_error - error ) > 1e-10 )
    DUNE_THROW(InvalidStateException, "Communication not working correctly on rank " << grid.comm().rank() << ": " << old_error << " -> " << error);
  else
    if (grid.comm().rank() == 0)
      std::cout << "Communication was successful and took " << timer.elapsed() << std::endl;

  return EXIT_SUCCESS;
}
