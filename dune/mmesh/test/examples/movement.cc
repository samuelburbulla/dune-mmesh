/** \example movement.cc
 * This is an example of how to move the interface within MMesh and get the connected components for mapping data.
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

    // write grids
    Dune::VTKWriter<GridView> vtkWriter( gridView );
    vtkWriter.write("movement-" + std::to_string(dim) + "d-0");

    Dune::VTKWriter<IGridView> ivtkWriter( igridView );
    ivtkWriter.write("movement-" + std::to_string(dim) + "d-interface-0");

    // define some movement
    using GlobalCoordinate = Dune::FieldVector<double, dim>;
    auto movement = []( GlobalCoordinate x )
    {
      GlobalCoordinate m = x;
      m -= GlobalCoordinate( 0.5 );
      m *= 2e-2 * M_PI;

      GlobalCoordinate s( 0.0 );
      s[0] = m[1];
      s[1] = -m[0];

      return s;
    };

    // time loop
    for ( int t = 1; t <= 10; t++ )
    {
      std::cout << "t = " << t << std::endl;

      using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper< IGridView >;
      Mapper mapper( igridView, Dune::mcmgVertexLayout() );

      std::vector<GlobalCoordinate> shifts ( mapper.size() );

      for ( const auto& vertex : vertices( igridView ) )
        shifts[ mapper.index(vertex) ] = movement( vertex.geometry().center() );

      grid.preAdapt();

      // 1. ensure mesh conformity after movement by adapting the mesh
      grid.ensureInterfaceMovement( shifts );

      // 2. mark grid elements
      grid.markElements();

      // 3. adapt because of indicator value
      grid.adapt();

      // 4. transfer data
      for ( const auto& element : elements( gridView ) )
      {
        if ( element.isNew() )
        {
          const auto& component = element.impl().connectedComponent();

          // we just do a check here
          double sum = 0.0;
          for( const auto& e : component.children() )
            sum += e.intersectionVolume( element );

          assert( std::abs( sum - element.geometry().volume() ) < 1e-8 );
        }
      }

      // 5. update the shifts because the interface grid might have changed
      mapper.update(igridView);
      shifts.resize( mapper.size() );
      for ( const auto& vertex : vertices( igridView ) )
        shifts[ mapper.index(vertex) ] = movement( vertex.geometry().center() );

      // 6. actually move interface
      grid.moveInterface( shifts );

      grid.postAdapt();

      // write grids
      vtkWriter.write("movement-" + std::to_string(dim) + "d-" + std::to_string(t));
      ivtkWriter.write("movement-" + std::to_string(dim) + "d-interface-" + std::to_string(t));
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
