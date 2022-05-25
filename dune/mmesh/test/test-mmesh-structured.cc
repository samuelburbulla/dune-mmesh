// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/math.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

//! Helper method for checking multiple properties given in a list with verbose output
template<class T>
void checkProperties( const std::string &name, std::initializer_list< std::initializer_list<T> > data )
{
  std::cout << "Check " << std::setfill(' ') << std::setw(30) << std::left << name;

  for( std::vector<T> datavector : data )
  {
    const T& resultValue = datavector[0];
    const T& testValue = datavector[1];

    bool correct = ( resultValue == testValue );

    if ( !correct )
    {
      std::cout.precision(17);
      std::cout << " - failed: result [" << resultValue << "] was expected to be [" << testValue << "]" << std::endl;
    }
    assert( resultValue == testValue );
  }

  std::cout << " - correct" << std::endl;
}

//! Helper method for checking a single property
template<class T>
void checkProperty( const std::string &name, const T& resultValue, const T& testValue )
{
  return checkProperties( name, {{ resultValue, testValue }} );
}

/** Test-template main program. Instantiate a single expression
 * template and evaluate it on a simple grid.
 */
int main(int argc, char *argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    std::cout << "-- MMesh implementation test --" << std::endl;

    std::cout << "Build structured grid..." << std::endl;

    // Create MMesh
    // ------------
    const unsigned int dim = GRIDDIM;
    using Grid = MovingMesh< dim >;
    const unsigned int dimworld = Grid::dimension;
    using FieldType = Grid::FieldType;
    using GlobalPosition = FieldVector< FieldType, dimworld >;

    using GridFactory = Dune::MMeshStructuredGridFactory< Grid >;

    GlobalPosition lowerLeft(0.0), upperRight(1.0);
    std::array<unsigned int,dimworld> elements;
    unsigned int num_elements_per_dim = 5;
    elements.fill(num_elements_per_dim);

    GridFactory gridFactory(lowerLeft, upperRight, elements);

    // get grid reference
    Grid& mMesh = *gridFactory.grid();

    // ---------------------------------------

    std::cout << "Build LeafGridView and LeafIndexSet..." << std::endl;
    auto indexSet = mMesh.leafIndexSet();

    // Test MMesh index set
    checkProperty( "size of element index set",
                   indexSet.size(0),
                   std::size_t ( pow(num_elements_per_dim,dim)*factorial(dim) ) );
    checkProperty( "size of vertex index set",
                   indexSet.size(dim),
                   std::size_t ( pow(num_elements_per_dim+1,dim) ) );

    return EXIT_SUCCESS;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return EXIT_FAILURE;
  }
  catch (CGAL::Failure_exception &e){
    std::cerr << "CGAL reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return EXIT_FAILURE;
  }
}
