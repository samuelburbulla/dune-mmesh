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

/** Test-template main program. Instantiate a single expression
 * template and evaluate it on a simple grid.
 */
int main(int argc, char *argv[])
{
  try {
    std::cout << "-- MMesh test with interface for 3D --" << std::endl;

    std::cout << "Build grid from mimesh3d.msh file..." << std::endl;

    // Create MMesh
    // ------------
    using Grid = Dune::MovingMesh<3>;
    using GridFactory = Dune::GmshGridFactory< Grid >;
    GridFactory gridFactory( "grids/mimesh3d.msh" );
    Grid& grid = *gridFactory.grid();

    // ---------------------------------------

    const auto& gridView = grid.leafGridView();

    // Write grid
    VTKWriter<typename Grid::LeafGridView> vtkWriter( gridView );
    vtkWriter.write("test-mimesh-3d-0");

    // ---------------------------------------

    std::cout << "Size of vertex index set: " << gridView.size(3) << std::endl;

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
