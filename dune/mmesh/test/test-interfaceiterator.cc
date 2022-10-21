// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

/** Test the interface iterator
 */
int main(int argc, char *argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    static constexpr int dim = GRIDDIM;

    std::cout << "-- Interface iterator test for " + std::to_string(dim) + "D --" << std::endl;

    std::cout << "Build grid from mimesh" << dim << "d.msh file..." << std::endl;

    // Create MMesh
    // ------------
    using Grid = Dune::MovingMesh<dim>;
    using GridFactory = Dune::GmshGridFactory<Grid>;
    GridFactory gridFactory( "grids/mimesh" + std::to_string(dim) + "d.msh" );
    Grid& grid = *gridFactory.grid();

    // ---------------------------------------

    const auto& gridView = grid.leafGridView();

    // Write grid
    VTKWriter<typename Grid::LeafGridView> vtkWriter( gridView );
    vtkWriter.write("test-interfaceiterator-" + std::to_string(dim) + "d-0");

    // ---------------------------------------
    int numberOfInterfaceSegments = 0;
    // Iterate over interface
    for ( const auto& segment : interfaceElements( gridView ) )
    {
      std::cout << segment.geometry().center() << std::endl;
      numberOfInterfaceSegments++;
    }

    int expected = dim == 2 ? 6 : 4;

    if( numberOfInterfaceSegments != expected )
      DUNE_THROW( GridError, "There are " << numberOfInterfaceSegments << " interface segments instead of " << expected << " expected!" << std::endl );

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
