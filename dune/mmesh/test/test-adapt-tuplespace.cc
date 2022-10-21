#include <config.h>

#include <iostream>


#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/blockvectordiscretefunction.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/io/file/dataoutput.hh>

#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fempy/grid/adaptation.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/capabilities.hh>

#include <dune/fem/io/parameter.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;
using namespace Fem;

typedef Dune::MovingMesh<2> MyGridType;
typedef AdaptiveLeafGridPart< MyGridType > GridPartType;

static const std::string usingSpaceName("Using TupleDiscreteFunctionSpace");
typedef Dune::Fem::FunctionSpace< MyGridType::ctype, double, 2, 2 > FuncSpace1;
typedef Dune::Fem::FunctionSpace< MyGridType::ctype, double, 2, 1 > FuncSpace2;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FuncSpace1, GridPartType, 2 > DiscreteFunctionSpaceType1;
typedef Dune::Fem::DiscontinuousGalerkinSpace< FuncSpace2, GridPartType, 1, CachingStorage > DiscreteFunctionSpaceType2;
typedef Dune::Fem::TupleDiscreteFunctionSpace< DiscreteFunctionSpaceType1, DiscreteFunctionSpaceType2 > DiscreteFunctionSpaceType;

typedef typename DiscreteFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
typedef DofManager< MyGridType > DofManagerType;
typedef AdaptationManager< MyGridType, RestrictProlongDefault< DiscreteFunctionType > > AdaptationManagerType;

// ***********************************************************
// the exact solution to the problem for EOC calculation
struct ExactSolution
: public Fem::Function< FunctionSpaceType, ExactSolution >
{
  typedef FunctionSpaceType::RangeType RangeType;
  typedef FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef FunctionSpaceType::DomainType DomainType;

  void evaluate ( const DomainType &x, RangeType &ret ) const
  {
    ret = 2.; // maximum of function is 2
    for( int i = 0; i < DomainType::dimension; ++i )
      ret *= sin( x[ i ]*(1.0 -x[ i ])*4.);
  }

  void evaluate ( const DomainType &x, RangeFieldType time, RangeType &ret ) const
  {
    evaluate( x, ret );
  }
};

// ********************************************************************
void adapt( MyGridType &grid, DiscreteFunctionType &solution, int step )
{
  const DiscreteFunctionType::DiscreteFunctionSpaceType &space = solution.space();
  RestrictProlongDefault<DiscreteFunctionType> rp(solution);
  rp.setFatherChildWeight(DGFGridInfo< MyGridType >::refineWeight());

  AdaptationManagerType adop(grid,rp);

  std::string message;

  int mark = 1;
  int count = std::abs(step);

  if(step < 0)
    mark = -1;

  if( Parameter :: verbose() )
  {
    std::cout << "Grid leaf size:             " << grid.size( 0 ) << std::endl;
    std::cout << "AdaptiveLeafIndexSet.size:  " << space.indexSet().size( 0 ) << std::endl;
  }

  for( int i = 0; i < count; ++i )
  {
    for( const auto& entity : space)
      grid.mark( mark, entity );
    adop.adapt();
    std::cout << message << std::endl;
  }

  if( Parameter :: verbose() )
  {
    std::cout << "Grid leaf size:             " << grid.size( 0 ) << std::endl;
    std::cout << "AdaptiveLeafIndexSet.size:  " << space.indexSet().size( 0 ) << std::endl;
  }

}
// ********************************************************************
double algorithm ( MyGridType &grid, DiscreteFunctionType &solution, int step, int turn, int writestep )
{
  const unsigned int order = solution.space().order();

  DiscreteFunctionType tmp ( solution );

  ExactSolution f;
  auto gridFunc = gridFunctionAdapter(f, solution.space().gridPart(), 2);
  {
    Dune::Fem::interpolate( gridFunc, solution );
    Dune :: Fem :: L2Norm< GridPartType > l2norm ( solution.space().gridPart(), 2*order+2 ) ;
    double new_error = l2norm.distance( f, solution );
    std::cout << "Before: " << new_error;
  }

  // output
  Dune::Fem::VTKIO< GridPartType > vtkIO( solution.space().gridPart(), Dune::VTK::nonconforming );
  vtkIO.addVertexData( gridFunc );
  vtkIO.addVertexData( solution );
  if (writestep == 0)
    vtkIO.write( "test-adapt-tuplespace-0" );

  adapt(grid, solution, step);

  // tmp solution should be zero after adapt
  double tmpErr = tmp.normSquaredDofs();
  if( tmpErr > 0 )
  {
    // return big error
    return 1e80;
  }

  // output
  vtkIO.write( "test-adapt-tuplespace-" + std::to_string(writestep+1) );

  // calculation L2 error on refined grid
  // pol ord for calculation the error chould by higher than
  // pol for evaluation the basefunctions
  Dune :: Fem :: L2Norm< GridPartType > l2norm ( solution.space().gridPart(), 2*order+2 ) ;
  double error = l2norm.distance( gridFunc, solution );
  std::cout << "After: " << error << std::endl;

  //! perform l2-projection to refined grid
  Dune::Fem::interpolate( gridFunc, solution );
  double new_error = l2norm.distance( gridFunc, solution );
  std::cout << "Interpolation on new grid: " << new_error << std::endl << std::endl;

  return error;
}

//**************************************************
//
//  main programm, run algorithm twice to calc EOC
//
//**************************************************
int main( int argc, char *argv[] )
try {
  MPIManager :: initialize( argc, argv );

  // std::string paramFile( paramName );
  std::cout << usingSpaceName << std::endl;

  // append parameter
  Parameter :: append( argc , argv );
  // Parameter :: append( paramFile );

  int ml = 2 ; // default value = 2
  ml = Parameter :: getValue ("lagrangeadapt.maxlevel", ml);

  std::vector<double> error(ml);

  // create grid from DGF file
  std::stringstream gridfilestr;
  gridfilestr << "grids/2dcoarse.dgf";

  std::string gridfile;
  Dune::Fem::Parameter::get( "fem.io.macrogrid", gridfilestr.str(), gridfile );

  Dune::Fem::Parameter::append( "fem.verboserank", 0 );
  Dune::Fem::Parameter::append( "fem.adaptation.method", "callback" );
  Dune::Fem::Parameter::append( "fem.prefix","output" );
  Dune::Fem::Parameter::append( "fem.io.savecount", "1" );
  Dune::Fem::Parameter::append( "fem.io.outputformat", "vtk-cell" );

  // construct macro using the DGF Parser
  Dune::GridPtr< MyGridType > gridPtr( gridfile );
  MyGridType& grid = *gridPtr;

  const int step = 1;

  GridPartType part ( grid );
  DiscreteFunctionSpaceType space( part );

  // threshold for EOC difference to predicted value
  const double eocThreshold = Parameter :: getValue("adapt.eocthreshold", double(0.2) );

  const bool isLocallyAdaptive = Dune::Fem::Capabilities::isLocallyAdaptive< MyGridType > :: v ;

  DiscreteFunctionType solution ( "sol", space );
  solution.clear();
  std::cout << "------------   Refining   ------------" << std::endl;
  for(int i=0; i<ml; i+=1)
  {
    error[i] = algorithm ( grid , solution, step, (i==ml-1), i);
    if (i>0)
    {
      if ( isLocallyAdaptive )
      {
        double eoc = log( error[i-1]/error[i]) / M_LN2;
        std::cout << "EOC = " << eoc << std::endl;
        if( std::abs( eoc - (space.order()+eocThreshold) ) < 0 )
        {
          DUNE_THROW(InvalidStateException,"EOC check of refinement failed");
        }
      }
      else
        std::cout << "no EOC for non-adaptive grid" << std::endl;
    }
  }
  std::cout << "------------   Coarsening   ------------" << std::endl;
  for(int i=ml-1; i>=0; i-=1)
  {
    error[i] = algorithm ( grid , solution,-step, 1, ml+1-i);
    if (i<ml-1)
    {
      if( isLocallyAdaptive )
      {
        double eoc = log( error[i+1]/error[i]) / M_LN2;
        std::cout << "EOC = " << eoc << std::endl;
        if( std::abs( eoc + (space.order()+eocThreshold) ) < 0 )
        {
          DUNE_THROW(InvalidStateException,"EOC check of coarsening failed");
        }
      }
      else
        std::cout << "no EOC for non-adaptive grid" << std::endl;
    }
  }
  return 0;
}
catch( const Dune :: Exception &exception )
{
  std :: cerr << exception << std :: endl;
  return 1;
}
