// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <unordered_set>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

/** Test-template main program.
 */
int main(int argc, char *argv[])
{
  try {
    std::cout << "-- Colliding test for 2D --" << std::endl;

    static constexpr int dim = 2;

    // Create MMesh
    // ------------
    std::cout << "Build grid using .msh file..." << std::endl;
    using Grid = Dune::MovingMesh<dim>;
    using GridFactory = Dune::GmshGridFactory<Grid>;
    GridFactory gridFactory( "grids/colliding2d.msh" );
    Grid& grid = *gridFactory.grid();

    // ---------------------------------------

    const auto& gridView = grid.leafGridView();

    using IGrid = typename Grid::InterfaceGrid;
    const IGrid& igrid = grid.interfaceGrid();
    const auto& igridView = igrid.leafGridView();

    // Write grid
    VTKWriter<typename Grid::LeafGridView> vtkWriter( gridView );
    vtkWriter.write("test-colliding-2d-0");

    VTKWriter<typename IGrid::LeafGridView> ivtkWriter( igridView );
    ivtkWriter.write("test-colliding-interface-2d-0");

    // do some movement
    using GlobalCoordinate = Dune::FieldVector<double, dim>;

    std::unordered_set<std::size_t> upperVertices;
    for( const auto& vertex : vertices( igridView ) )
    {
      const GlobalCoordinate& x = vertex.geometry().center();
      if ( x[1] > 0.5 )
        upperVertices.insert( igrid.globalIdSet().id( vertex ) );
    }

    // define the movement
    auto movement = [upperVertices, &igrid]( const auto& vertex )
    {
      double speed = 1.01e-2;
      GlobalCoordinate m (0.0);

      if ( upperVertices.count( igrid.globalIdSet().id( vertex ) ) > 0 )
        m[1] = -speed;
      else
        m[1] = speed;

      return m;
    };

    for ( int t = 1; t <= 20; t++ )
    {
      // prepare grid
      grid.preAdapt();

      std::vector<GlobalCoordinate> shifts( igridView.size(dim-1) );
      for( const auto& vertex : vertices( igridView ) )
      {
        // obtain the constant shift
        const std::size_t idx = igrid.leafIndexSet().index( vertex );
        shifts[ idx ] = movement( vertex );
      }

      grid.ensureInterfaceMovement( shifts );

      while ( grid.preAdapt() )
      {
        grid.adapt();
        grid.postAdapt();

        shifts.resize( igridView.size(dim-1) );
        for( const auto& vertex : vertices( igridView ) )
        {
          // obtain the constant shift
          const std::size_t idx = igrid.leafIndexSet().index( vertex );
          shifts[ idx ] = movement( vertex );
        }
        grid.ensureInterfaceMovement( shifts );
      }

      grid.moveInterface( shifts );

      // delete markers in grid
      grid.postAdapt();

      // Write grid
      vtkWriter.write("test-colliding-2d-" + std::to_string(t));
      ivtkWriter.write("test-colliding-interface-2d-" + std::to_string(t));
    }

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
  catch (std::exception &e){
    std::cerr << "STD reported error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return EXIT_FAILURE;
  }
}
