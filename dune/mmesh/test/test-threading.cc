// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <config.h>
#include <iostream>
#include <chrono>

#include <dune/mmesh/mmesh.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/threads/threaditerator.hh>
#include <dune/fem/misc/threads/threadpool.hh>

using namespace Dune;
using namespace std::chrono;

int main(int argc, char *argv[])
{
  try {
    Fem::MPIManager::initialize( argc, argv );

    std::cout << "-- MMesh threading test --" << std::endl;

    std::string gridfile = "grids/ellipse2d.msh";
    if (GRIDDIM == 3)
      gridfile = "grids/sphere3d.msh";

#if !INTERFACE
    using Grid = MovingMesh<GRIDDIM>;
    using GridFactory = GmshGridFactory< Grid >;
    GridFactory gridFactory( gridfile );
    Grid& grid = *gridFactory.grid();
#else
    using HGrid = MovingMesh<GRIDDIM>;
    using GridFactory = GmshGridFactory< HGrid >;
    GridFactory gridFactory( gridfile );
    HGrid& hgrid = *gridFactory.grid();
    using Grid = typename HGrid::InterfaceGrid;
    Grid& grid = hgrid.interfaceGrid();
#endif
    using GridPart = Fem::LeafGridPart< Grid >;
    GridPart gridPart( grid );

    Fem::ThreadManager::setNumThreads(2);

    Fem::ThreadIterator< GridPart > threadIterator( gridPart );
    threadIterator.update();

    auto execute = [&](int numThreads)
    {
      Fem::ThreadManager::setNumThreads(numThreads);

      Fem::ThreadIterator< GridPart > threadIterator( gridPart );
      threadIterator.update();

      std::cout << "Number of threads: " << Fem::ThreadManager::numThreads();

      double volume = 0.0;

      std::mutex mutex;
      auto doEval = [&threadIterator, &volume, &mutex] ()
      {
        double myVol = 0.0;
        const auto end = threadIterator.end();
        for( auto it = threadIterator.begin(); it != end; ++it )
        {
          const auto& entity = *it;
          myVol += entity.geometry().volume();
        }

        {
          std::lock_guard guard ( mutex );
          volume += myVol;
        }
      };

      auto start = high_resolution_clock::now();
      Fem::ThreadPool::run( doEval );
      auto stop = high_resolution_clock::now();

      auto duration = duration_cast<microseconds>(stop - start);
      std::cout << " - Took " << duration.count() << std::endl;

#if !INTERFACE
      double expVol = 1.0;
#else
      double expVol = (GRIDDIM == 2) ? 1.4495 : 0.2652;
#endif
      assert(std::abs(volume - expVol) < 1e-4);
    };

    execute(1);
    execute(2);
    execute(4);
    execute(8);

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
