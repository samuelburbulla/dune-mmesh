// include configurations options

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Includes from the IOStream Library
// ----------------------------------

#include <iostream>
#include <sstream>

// Includes from DUNE-FEM
// ----------------------

// include Lagrange discrete function space
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/function/localfunction/bindable.hh>
#include <dune/fem/io/file/vtkio.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>

#include <dune/fem/space/common/interpolate.hh>
#include <dune/common/timer.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>

#include "massoperator.hh"
#include <dune/mmesh/mmesh.hh>

// Global Type Definitions
// -----------------------

using BulkGridType = Dune::MovingMesh<GRIDDIM>;
#ifdef INTERFACE
  using GridType = typename BulkGridType::InterfaceGrid;
#else
  using GridType = BulkGridType;
#endif
static constexpr int polOrder = 1;

typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
typedef Dune::Fem::FunctionSpace< typename GridType::ctype, typename GridType::ctype, GridType::dimensionworld, 1 > SpaceType;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, polOrder > DiscreteSpaceType;

typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteSpaceType > DiscreteFunctionType;
typedef Dune::Fem::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType > LinearOperatorType;
typedef Dune::Fem::KrylovInverseOperator< DiscreteFunctionType > InverseOperatorType;

// TestGrid
// --------

class TestGrid
{
  typedef TestGrid ThisType;
public:
  typedef BulkGridType HGridType;

protected:
  TestGrid ( const std::string& name )
  : gridptr_( name )
  {
    gridptr_->loadBalance();
  }

public:
  TestGrid ( const ThisType & ) = delete;

  ThisType &operator= ( const ThisType & ) = delete;

  static ThisType &instance ( const std::string& name )
  {
    static ThisType staticInstance( name );
    return staticInstance;
  }

  static HGridType &grid ( const std::string name = macroGridName() )
  {
    return *(instance( name ).gridptr_);
  }

  static int refineStepsForHalf ()
  {
    return Dune::DGFGridInfo< HGridType >::refineStepsForHalf();
  }

protected:
  static std::string macroGridName ()
  {
    std::ostringstream s;
    s << "grids/junction2d.msh";
    return s.str();
  }

  Dune::GridPtr< HGridType > gridptr_;
};

typedef MassOperator< DiscreteFunctionType, LinearOperatorType > MassOperatorType;

// Function to Project
// -------------------
template <class GridPart, class RangeType>
struct Function : public Dune::Fem::BindableGridFunction< GridPart, RangeType >
{
  typedef Dune::Fem::BindableGridFunction<GridPart, RangeType > Base;
  using Base::Base;

  template <class Point>
  void evaluate ( const Point &p, RangeType &value ) const
  {
    auto x = Base::global(p);
    value = 1.;
    for( int k = 0; k < x.dimension; ++k )
      value *= sin( M_PI * x[k] );
  }
  template <class Point>
  void jacobian ( const Point &p, typename Base::JacobianRangeType &jacobian ) const
  {
    auto x = Base::global(p);
    for( int j = 0; j < x.dimension; ++j )
    {
      // jacobian has only one row, calc j-th column
      jacobian[0][j] = M_PI;
      for( int k = 0; k < x.dimension; ++k )
        jacobian[0][j] *= (j == k ? cos( M_PI*x[k] ) : sin( M_PI*x[k] ));
    }
  }
  template <class Point>
  void hessian ( const Point &p, typename Base::HessianRangeType &h ) const
  {
  }

  unsigned int order() const { return 4; }
  std::string name() const { return "MyFunction"; } // needed for output
};

// Algorithm
// ---------

struct Algorithm
{
  typedef Dune::FieldVector< double, 2 > ErrorType;
  typedef Function< GridPartType, SpaceType::RangeType > FunctionType;

  explicit Algorithm ( GridType &grid )
  : grid_( grid )
  {}

  ErrorType operator() ( int step );

  void nextMesh ()
  {
    grid_.globalRefine( 2 );
    grid_.loadBalance();
  }

private:
  GridType& grid_;
};

inline Algorithm::ErrorType Algorithm::operator() ( int step )
{
  GridPartType gridPart( grid_ );
  DiscreteSpaceType space( gridPart );
  DiscreteFunctionType solution( "solution", space );
  solution.clear();

  MassOperatorType massOperator( space );
  DiscreteFunctionType rhs( "rhs", space );
  FunctionType function_( gridPart );

  InverseOperatorType inverseOperator;
  inverseOperator.bind( massOperator );

  // apply solver
  Dune::Timer timer;
  timer.start();
  massOperator.assembleRHS( function_, rhs );
  if ( grid_.comm().rank() == 0 )
    std::cout << "assmbleRHS took " << timer.elapsed() << std::endl;

  timer.reset();
  inverseOperator( rhs, solution );
  if ( grid_.comm().rank() == 0 )
    std::cout << "inverseOperator took " << timer.elapsed() << std::endl;

  // write vtk
  Dune::Fem::VTKIO< GridPartType > vtkIO( gridPart, Dune::VTK::nonconforming );
  vtkIO.addVertexData( function_ );
  vtkIO.addVertexData( solution );
  vtkIO.write( "test-l2projection-" + std::to_string(step) );

  // return error
  ErrorType error;
  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  Dune::Fem::H1Norm< GridPartType > h1norm( gridPart );

  error[ 0 ] = l2norm.distance( function_, solution );
  error[ 1 ] = h1norm.distance( function_, solution );
  if ( grid_.comm().rank() == 0 )
    std::cout << error << std::endl;
  return error;
}


void run( GridType& grid )
{
  const int nrSteps = 1;

  Algorithm algorithm( grid );
  std::vector< typename Algorithm::ErrorType > error( nrSteps );
  error[ 0 ] = algorithm( 0 );
  for( int step = 1; step < nrSteps; ++step )
  {
    algorithm.nextMesh();
    error[ step ] = algorithm( step );
  }

  for( int step = 1; step < nrSteps; ++step )
  {
    double l2eoc = log( error[ step ][ 0 ] / error[ step -1 ][ 0 ] ) / log( 0.5 );
    double h1eoc = log( error[ step ][ 1 ] / error[ step -1 ][ 1 ] ) / log( 0.5 );
    if ( grid.comm().rank() == 0 )
      std::cout << "EOC: " << l2eoc << "   " << h1eoc << std::endl;

    if( std::abs( l2eoc-1 - polOrder ) > 0.2 )
    {
      DUNE_THROW(Dune::InvalidStateException, "EOC check of solving mass matrix system failed L2eoc " << l2eoc << " " << polOrder);
    }

    if( std::abs( h1eoc - polOrder ) > 0.2 )
    {
      DUNE_THROW(Dune::InvalidStateException, "EOC check of solving mass matrix system failed H1eoc " << h1eoc << " " << polOrder);
    }
  }
}


// Main Program
// ------------

int main ( int argc, char **argv )
try
{
  // initialize MPI manager and PETSc
  Dune::Fem::MPIManager::initialize( argc, argv );

  // add command line parameters to global parameter table
  Dune::Fem::Parameter::append( argc, argv );
  // append parameters from the parameter file
  Dune::Fem::Parameter::append( (argc < 2) ? "parameter" : argv[ 1 ] );

  BulkGridType &bulkGrid = TestGrid::grid();

#ifdef INTERFACE
  GridType& grid = bulkGrid.interfaceGrid();
#else
  GridType& grid = bulkGrid;
#endif

  if ( grid.comm().rank() == 0 )
    std::cout << grid.size(0) << " cells" << std::endl;

  run( grid );

  Dune::Fem::Parameter::write( "parameter.log" );

  return 0;
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;
}
