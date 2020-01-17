// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dune/grid/io/file/printgrid.hh>
#include <dune/mmesh/mmesh.hh>

using namespace Dune;

/** Test a FV scheme on the interface grid
 */
int main(int argc, char *argv[])
{
  try {
    std::cout << "-- MMesh interface grid fv test --" << std::endl;

    // Create MMesh
    // ------------
    static constexpr int dim = GRIDDIM;
    using Grid = Dune::MovingMesh<dim>;

    using GridFactory = Dune::GmshGridFactory< Grid >;
    GridFactory gridFactory( "grids/mifvmesh" + std::to_string(dim) + "d.msh" );

    Grid& grid = *gridFactory.grid();

    auto gridView = grid.leafGridView();

    // ---------------------------------------

    using InterfaceGrid = typename Grid::InterfaceGrid;
    const InterfaceGrid& igrid = grid.interfaceGrid();

    // ---------------------------------------

    using GlobalCoordinate = Dune::FieldVector<double, dim>;

    std::unordered_set<std::size_t> moveVertices;
    for( const auto& vertex : vertices( igrid.leafGridView() ) )
    {
      const GlobalCoordinate& x = vertex.geometry().center();
      if ( std::abs( x[0] - 0.1 ) < 1e-6 )
        moveVertices.insert( igrid.globalIdSet().id( vertex ) );
    }

    // define the movement
    auto movement = [moveVertices, &igrid]( const auto& vertex )
    {
      static const double speedx = 2.5;
      static const double speedy = 0.0;

      const auto& x = vertex.geometry().center();

      GlobalCoordinate m (0.0);
      if ( moveVertices.count( igrid.globalIdSet().id( vertex ) ) > 0 )
      {
        m[0] = speedx;
        m[1] = speedy * std::cos( 10. * M_PI * x[0] );
      }
      return m;
    };
    // ---------------------------------------

    VTKSequenceWriter<typename Grid::LeafGridView> vtkWriterBulk( gridView, "test-fvoninterfacegrid-2d-bulk", "", "" );
    vtkWriterBulk.write( 0.0 );

    std::cout << "Build LeafGridView and LeafIndexSet of interface..." << std::endl;
    auto igridView = igrid.leafGridView();
    auto& indexSet = igrid.leafIndexSet();
    auto& idSet = igrid.globalIdSet();

    std::vector<double> c( indexSet.size(0), 0.0 );

    VTKSequenceWriter<decltype(igridView)> vtkWriterInterface( igridView, "test-fvoninterfacegrid-2d-interface", "", "" );
    vtkWriterInterface.addCellData( c, "c" );
    vtkWriterInterface.write( 0.0 );

    const double dtInitial = 0.01;
    const double tEnd = 0.12;

    auto b = []( auto x ){ return x[0] < 1e-8; };
    FieldVector<double, 2> v ( {1.0, 0.0} );


    double dt = dtInitial;
    for( double t = 0.0; t <= tEnd + 1e-8; t += dt )
    {
      // set timestep
      dt = dtInitial;

      // printGrid( igrid,  MPIHelper::instance(argc,argv) );
      std::cout << "t = " << t << std::endl;

      // ==================
      // Move the interface
      // ==================

      // store interface data
      std::unordered_map<Dune::Impl::MultiId, double> data;
      for (const auto& e : elements( igridView ))
        data.insert( std::make_pair( idSet.id(e), c[ indexSet.index(e) ] * e.geometry().volume() ) );

      grid.preAdapt();
      std::vector<GlobalCoordinate> shifts( igridView.size(dim-1) );
      for( const auto& vertex : vertices( igridView ) )
      {
        // obtain the constant shift
        const std::size_t idx = igrid.leafIndexSet().index( vertex );
        shifts[ idx ] = movement(vertex);
        shifts[ idx ] *= dt;
      }

      grid.ensureInterfaceMovement( shifts );

      shifts.resize( igridView.size(dim-1) );
      for( const auto& vertex : vertices( igridView ) )
      {
        // obtain the constant shift
        const std::size_t idx = igrid.leafIndexSet().index( vertex );
        shifts[ idx ] = movement(vertex);
        shifts[ idx ] *= dt;
      }

      grid.moveInterface( shifts );

      c.resize( indexSet.size(0) );

      grid.postAdapt();

      // restore interface data
      for (const auto& e : elements( igridView ))
      {
        auto it = data.find( idSet.id(e) );
        if ( it != data.end() )
          c[ indexSet.index(e) ] = it->second / e.geometry().volume();
        else
        {
          assert( e.impl().hasConnectedComponent() );
          auto father = e.impl().connectedComponent();
          c[ indexSet.index(e) ] = data[ idSet.id(father) ] / ( 2.0 * e.geometry().volume() ); // a father has always to equally sized children
        }
      }

      // ==================
      //      FV step
      // ==================
      std::vector<double> update( indexSet.size(0) );
      for(const auto& e : elements(igridView))
      {
        const auto eIdx = indexSet.index(e);
        const auto eVolume = e.geometry().volume();

        for(const auto is : intersections(igridView, e))
        {
          if ( is.neighbor() )
          {
            const auto neighbor = is.impl().outside();
            const auto nbIdx = indexSet.index( neighbor );
            const auto n = is.centerUnitOuterNormal();
            const auto isVol = is.geometry().volume();

            update[ eIdx ] -= dt / eVolume * isVol * (v * n) * ( (v * n) > 0 ? c[ eIdx ] : c[ nbIdx ] ) / is.impl().numOutside();
          }
          else
          {
            const auto n = is.centerUnitOuterNormal();
            const auto isVol = is.geometry().volume();
            update[ eIdx ] -= dt / eVolume * isVol * (v * n) * ( (v * n) > 0 ? c[ eIdx ] : b( is.geometry().center() ) );
          }
        }
      }

      for ( std::size_t i = 0; i < update.size(); ++i )
        c[i] += update[i];

      // write data
      vtkWriterBulk.write( t + dt );
      vtkWriterInterface.write( t + dt );
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
  catch(std::exception& e) {
    std::cerr << "std exception thrown: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  catch(...) {
    std::cerr << "Unknown exception thrown!" << std::endl;
    return EXIT_FAILURE;
  }
}
