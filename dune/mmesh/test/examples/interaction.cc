/** \example movement.cc
 * This is an example of how to interact between entities of the MMesh and the InterfaceGrid.
 *
 * We use the same mesh file as in interfacegrid.cc.
 */
#include <config.h>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-grid includes
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// dune-mmesh includes
#include <dune/mmesh/mmesh.hh>

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  try
  {
    static constexpr int dim = GRIDDIM;

    using Grid = Dune::MovingMesh<dim>;

    using GridFactory = Dune::GmshGridFactory<Grid>;
    GridFactory gridFactory( "grids/horizontal" + std::to_string(dim) + "d.msh" );

    Grid& grid = *gridFactory.grid();

    using GridView = typename Grid::LeafGridView;
    const GridView& gridView = grid.leafGridView();

    using IGrid = typename Grid::InterfaceGrid;
    const IGrid& igrid = grid.interfaceGrid();

    using IGridView = typename IGrid::LeafGridView;
    const IGridView& igridView = igrid.leafGridView();


    // Seen from bulk
    for ( const auto& element : elements( gridView ) )
    {
      for ( const auto& intersection : intersections( gridView, element ) )
      {
        if( grid.isInterface( intersection ) )
        {
          // obtain interface entity coinciding with intersection
          const auto interfaceEntity = grid.asInterfaceEntity( intersection );

          // print some information
          std::cout << interfaceEntity.geometry().center() << std::endl;
        }
      }
    }

    // Seen from interface
    for ( const auto& element : elements( igridView ) )
    {
      // obtain intersection coinciding with interface entity (seen from lower entity index to higher)
      const auto intersection = grid.asIntersection( element );

      // print some information
      std::cout << intersection.centerUnitOuterNormal() << std::endl;
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
