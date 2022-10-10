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
#define POLORDER 2
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

struct IDataOutputParameters
: public Dune::Fem::LocalParameter< Dune::Fem::DataOutputParameters, IDataOutputParameters >
{
  IDataOutputParameters ( const int step )
  : step_( step )
  {}

  IDataOutputParameters ( const IDataOutputParameters &other )
  : step_( other.step_ )
  {}

  std::string prefix () const
  {
    std::stringstream s;
    s << "femadaptation-interface-" << step_ << "-";
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

  typedef typename GridType :: template Codim< 0 > :: Entity ElementType;
  typedef typename DiscreteFunctionSpaceType :: DomainType DomainType;

  Scheme( GridPartType &gridPart, const int step = 0 )
    : gridPart_( gridPart ),
      grid_( gridPart_.grid() ),
      discreteSpace_( gridPart_ ),
      solution_( "solution", discreteSpace_ ),
      restrictProlong_( solution_ ),
      adaptationManager_( gridPart_.grid(), restrictProlong_ ),
      step_( step )
  {
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
    return grid_.markElements();
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
    y[ 0 ] = x[ 0 ] * x[ 0 ] + x[ 1 ] * x[ 1 ];
  }
};

// algorithm
// ---------

template <class HGridType, class IGridType>
std::array<double, 2> algorithm ( HGridType &grid, IGridType &igrid, const int step )
{
  // we want to solve the problem on the leaf elements of the grid
  typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > GridPartType;
  GridPartType gridPart(grid);

  // and on the interface grid
  typedef Dune::Fem::AdaptiveLeafGridPart< IGridType > IGridPartType;
  IGridPartType igridPart(igrid);

  // use a scalar function space
  typedef Dune::Fem::FunctionSpace< double, double, HGridType :: dimensionworld, 1 > FunctionSpaceType;
  typedef Dune::Fem::FunctionSpace< double, double, IGridType :: dimensionworld, 1 > IFunctionSpaceType;

  typedef Scheme< GridPartType, FunctionSpaceType > SchemeType;
  SchemeType scheme( gridPart, step );
  typedef Scheme< IGridPartType, IFunctionSpaceType > ISchemeType;
  ISchemeType ischeme( igridPart, step );

  typedef Function< FunctionSpaceType > FunctionType;
  FunctionType f;
  typedef Function< IFunctionSpaceType > IFunctionType;
  IFunctionType i_f;

  typedef Dune::Fem::GridFunctionAdapter< FunctionType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution("exact solution", f, gridPart, 5 );
  typedef Dune::Fem::GridFunctionAdapter< IFunctionType, IGridPartType > IGridExactSolutionType;
  IGridExactSolutionType igridExactSolution("interface exact solution", i_f, igridPart, 5 );

  scheme.initialize( gridExactSolution );
  ischeme.initialize( igridExactSolution );

  // output
  typedef std::tuple< const typename SchemeType::DiscreteFunctionType *, GridExactSolutionType * > IOTupleType;
  typedef Dune::Fem::DataOutput< HGridType, IOTupleType > DataOutputType;
  IOTupleType ioTuple( &(scheme.solution()), &gridExactSolution);
  DataOutputType dataOutput( grid, ioTuple, DataOutputParameters( step ) );
  dataOutput.write();

  typedef std::tuple< const typename ISchemeType::DiscreteFunctionType *, IGridExactSolutionType * > IIOTupleType;
  typedef Dune::Fem::DataOutput< IGridType, IIOTupleType > IDataOutputType;
  IIOTupleType iioTuple( &(ischeme.solution()), &igridExactSolution);
  IDataOutputType idataOutput( igrid, iioTuple, IDataOutputParameters( step ) );
  idataOutput.write();

  for( int time = 0; time <= 2; time++ )
  {
    if( Dune::Fem::MPIManager::rank() == 0 )
      std::cout << "step: " << time << std::endl;

    // mark element for adaptation
    scheme.mark( time );

    // adapt grid
    scheme.adapt();
    ischeme.adapt();

    // data I/O
    dataOutput.write();
    idataOutput.write();

    // make range smaller
    grid.indicator().maxH() *= 0.8;
    grid.indicator().minH() *= 1.2;
  }

  // compute error
  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  Dune::Fem::L2Norm< IGridPartType > il2norm( igridPart );
  std::array<double, 2> ret {{ l2norm.distance( gridExactSolution, scheme.solution() ),
           il2norm.distance( igridExactSolution, ischeme.solution() ) }};

  // test addInterface
  auto add = [&]()
  {
    for (const auto& e : elements(gridPart))
      for (const auto& i : intersections(gridPart, e))
         if (!grid.isInterface(i))
         {
           std::cout << "Add interface at " << i.geometry().center() << std::endl;
           grid.addInterface(i);
           return;
         }
  };
  add();
  ischeme.adapt();
  idataOutput.write();

  return ret;
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
  gridfilestr << "grids/line2d.msh";

  std::string gridfile;
  Dune::Fem::Parameter::get( "fem.io.macrogrid", gridfilestr.str(), gridfile );

  // construct macro using the DGF Parser
  Dune::GridPtr< HGridType > gridPtr( gridfile );
  HGridType& grid = *gridPtr ;

  grid.indicator().maxH() = 0.75;
  grid.indicator().minH() = 0.05;
  grid.indicator().factor() = 1.0;

  auto errors = algorithm( grid, grid.interfaceGrid(), 0 );
  std::cout << "Error bulk: " << errors[0] << std::endl;
  std::cout << "Error interface: " << errors[1] << std::endl;
  assert( errors[0] < 1e-12 );
  assert( errors[1] < 1e-12 );

  return 0;
}
catch( const Dune::Exception &exception )
{
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
#else
int main ( int argc, char **argv )
{
  return 0;
}
#endif
