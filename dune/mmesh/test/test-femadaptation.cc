#include <config.h>

#include <iostream>
#include <sstream>
#include <string>
#include <tuple>

#if HAVE_DUNE_FEM

#include <dune/mmesh/mmesh.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/dataoutput.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/misc/l2norm.hh>

#ifndef POLORDER
#define POLORDER 1
#endif


// DataOutputParameters
// --------------------

struct DataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, DataOutputParameters >
{
  DataOutputParameters ( const int step )
  : step_( step )
  {}

  DataOutputParameters ( const DataOutputParameters &other )
  : step_( other.step_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "femadaptation-" << step_ << "-";
    return s.str();
  }

private:
  int step_;
};

// Scheme
// ------

template < class GridPart, class FunctionSpace >
struct Scheme
{
  typedef GridPart GridPartType;
  typedef typename GridPartType::GridType GridType;
  typedef FunctionSpace FunctionSpaceType;

#ifdef LAGRANGE
  typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;
#else
  typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, POLORDER > DiscreteFunctionSpaceType;
#endif

  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

  typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType >  RestrictionProlongationType;

  typedef Dune::Fem::AdaptationManager< GridType, RestrictionProlongationType > AdaptationManagerType;

  typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;

  Scheme( GridPartType &gridPart, const int step  = 0 )
    : gridPart_( gridPart ),
      grid_( gridPart_.grid() ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      restrictProlong_( solution_ ),
      adaptationManager_( gridPart_.grid(), restrictProlong_ ),
      step_( step )
  {
    if( discreteSpace_.begin() != discreteSpace_.end() )
      solution_.localFunction( *(discreteSpace_.begin()) )[0] = 0.;
    solution_.clear();
  }

  const DiscreteFunctionType &solution() const
  {
    return solution_;
  }

  template <class GF>
  void initialize( const GF& gridFunction )
  {
    Dune::Fem::interpolate( gridFunction, solution_ );
  }

  //! mark elements for adaptation
  bool mark ( double time ) const
  {
    int refine = 0;
    int coarse = 0;
    int total = 0;

    // loop over all elements
    for( const auto& entity : discreteSpace_ )
    {
      // find center
      auto center = entity.geometry().center();
      DomainType x = DomainType(-time);
      x += center;

      if( grid_.leafIndexSet().index( entity ) == 3 )
      {
        grid_.mark( 1, entity );
        refine++;
      }
      else if( grid_.leafIndexSet().index( entity ) == 17 )
      {
        grid_.mark( -1, entity );
        coarse++;
      }

      total++;
    }

    std::cout << "Marked " << refine << " for refine and " << coarse << " for coarse of total " << total << std::endl;

    return (coarse + refine > 0);
  }

  //! do the adaptation for a given marking
  void adapt()
  {
    // apply adaptation and load balancing
    adaptationManager_.adapt();
  }

protected:
  GridPartType &gridPart_; // grid part(view), e.g. here the leaf grid the discrete space is build with
  GridType &grid_;
  DiscreteFunctionSpaceType discreteSpace_; // discrete function space
  DiscreteFunctionType solution_;   // the unknown
  RestrictionProlongationType restrictProlong_ ; // local restriction/prolongation object
  AdaptationManagerType  adaptationManager_ ;    // adaptation manager handling adaptation
  const int step_;
};

template< class FunctionSpace >
struct Function : Dune::Fem::Function< FunctionSpace, Function< FunctionSpace > >
{
  void evaluate( const typename FunctionSpace::DomainType &x, typename FunctionSpace::RangeType &y ) const
  {
    y[ 0 ] = x[ 0 ] + x[ 1 ];
  }
};

// algorithm
// ---------

template <class HGridType>
double algorithm ( HGridType &grid, const int step )
{
  // we want to solve the problem on the leaf elements of the grid
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > GridPartType;
  GridPartType gridPart(grid);

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, HGridType :: dimensionworld, 1 > FunctionSpaceType;

  typedef Scheme< GridPartType, FunctionSpaceType > SchemeType;
  SchemeType scheme( gridPart, step );

  typedef Function< FunctionSpaceType > FunctionType;
  FunctionType f;
  typedef Dune::Fem::GridFunctionAdapter< FunctionType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", f, gridPart, 5 );

  scheme.initialize( gridExactSolution );

  // output
  typedef std::tuple< const typename SchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution);
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );
  dataOutput.write();

  for( int time = 0; time <= 5; time++ )
  {
    if( Dune::Fem::MPIManager::rank() == 0 )
      std::cout << "step: " << time << std::endl;

    // mark element for adaptation
    scheme.mark( time );

    // adapt grid
    scheme.adapt();

    // data I/O
    dataOutput.write();
  }

  // compute error
  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  return l2norm.distance( gridExactSolution, scheme.solution() );
}


int main ( int argc, char **argv )
try
{
  // initialize MPI, if necessary
  Dune::Fem::MPIManager::initialize( argc, argv );

  // append overloaded parameters from the command line
  Dune::Fem::Parameter::append( argc, argv );

  // append possible given parameter files
  for( int i = 1; i < argc; ++i )
    Dune::Fem::Parameter::append( argv[ i ] );
  // make sure output parameters are added
  Dune::Fem::Parameter::append( "fem.prefix","output" );
  Dune::Fem::Parameter::append( "fem.io.savecount", "1" );
  Dune::Fem::Parameter::append( "fem.io.outputformat", "vtk-cell" );
  Dune::Fem::Parameter::append( "fem.io.partitioning", "rank" );
  Dune::Fem::Parameter::append( "fem.loadbalancing.step", "1" );
  Dune::Fem::Parameter::append( "fem.adaptation.method", "callback" );

  // type of hierarchical grid
  typedef Dune::MovingMesh<2> HGridType;

  // create grid from DGF file
  std::stringstream gridfilestr;
  gridfilestr << "grids/" << HGridType::dimension << "dgrid.dgf";

  std::string gridfile;
  Dune::Fem::Parameter::get( "fem.io.macrogrid", gridfilestr.str(), gridfile );

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  auto error = algorithm( grid, 0 );
  std::cout << "Error: " << error << std::endl;

  assert( error < 1e-12 );

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
#endif
