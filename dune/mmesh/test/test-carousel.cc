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

/** Test-template main program.
 */
int main(int argc, char *argv[])
{
  static constexpr int dim = GRIDDIM;
  try {
    std::cout << "-- Carousel test for " << dim << "D --" << std::endl;

    // Create MMesh
    // ------------
    std::cout << "Build grid using .msh file..." << std::endl;
    using Grid = Dune::MovingMesh<dim>;
    using GridFactory = Dune::GmshGridFactory<Grid>;
    GridFactory gridFactory( "grids/horizontal" + std::to_string(dim) + "d.msh" );
    Grid& grid = *gridFactory.grid();
    // ---------------------------------------

    const auto& gridView = grid.leafGridView();

    // Write grid
    VTKWriter<typename Grid::LeafGridView> vtkWriter( gridView );
    vtkWriter.write("test-carousel-" + std::to_string(dim) + "d-0");

    // Obtain interfacegrid
    using IGrid = typename Grid::InterfaceGrid;
    const IGrid& igrid = grid.interfaceGrid();
    using IGridView = typename IGrid::LeafGridView;
    const IGridView& igridView = igrid.leafGridView();
    const auto& iindexSet = igrid.leafIndexSet();

    VTKWriter<IGridView> ivtkWriter( igridView );
    ivtkWriter.write("test-carousel-" + std::to_string(dim) + "d-interface-0");

    // do some movement

    using GlobalCoordinate = Dune::FieldVector<double, dim>;

    // define the movement
    auto movement = []( GlobalCoordinate x )
    {
      double speed = 1e-3;
      GlobalCoordinate m = x;
      m -= GlobalCoordinate( 0.5 );
      m *= speed * M_PI;

      GlobalCoordinate s( 0.0 );
      s[0] = m[1];
      s[1] = -m[0];

      // s -= m;

      return s;
    };

    for ( int t = 1; t <= 500; t++ )
    {
      // skip the loop in 3d until the remeshing is implemented in 3d
      if ( dim == 3 && t == 9 )
        return 0;

      std::cout << "t = " << t << std::endl;

      std::vector<GlobalCoordinate> shifts ( igridView.size(dim-1) );
      for( const auto& vertex : vertices( igridView ) )
        shifts[iindexSet.index(vertex)] = movement(vertex.geometry().center());

      grid.preAdapt();

      // 1. ensure movement
      grid.ensureInterfaceMovement( shifts );

      // 2a. mark elements
      grid.markElements();

      // 3a. adapt
      grid.adapt();

      // 4. transfer data ...
      for (const auto& element : elements(gridView))
        if ( element.isNew() )
        {
          double sum = 0.0;
          const auto& component = element.impl().connectedComponent();
          for( const auto& e : component.children() )
            sum += e.intersectionVolume( element );

          if( std::abs( sum - element.geometry().volume() ) > 1e-8 )
          {
            std::cout << "Cell at: " << element.geometry().center() << std::endl;
            std::cout << "Sum of intersection volumes " << sum << " should be " << element.geometry().volume() << std::endl;
            std::cout << "Component size: " << component.size() << std::endl;
            std::cout << "Corners at: " << std::endl;
            for( int i = 0; i < dim+1; ++i )
              std::cout << element.geometry().corner(i) << std::endl;
            assert(false);
          }
        }

      grid.postAdapt();

      shifts.resize( igridView.size(dim-1) );
      for( const auto& vertex : vertices( igridView ) )
        shifts[iindexSet.index(vertex)] = movement(vertex.geometry().center());

      // move vertices
      grid.moveInterface( shifts );

      // Write grid
      vtkWriter.write("test-carousel-" + std::to_string(dim) + "d-" + std::to_string(t));
      ivtkWriter.write("test-carousel-" + std::to_string(dim) + "d-interface-" + std::to_string(t));
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
