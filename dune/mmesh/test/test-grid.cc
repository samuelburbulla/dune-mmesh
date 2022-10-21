// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <chrono>
#include <thread>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/timer.hh>

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checktwists.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

/** Test-template main program. Instantiate a single expression
 * template and evaluate it on a simple grid.
 */
int main(int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);

  // Create Grid
  // ------------
  static constexpr int dim = GRIDDIM;
  using Grid = Dune::MovingMesh<dim>;

  using GridFactory = Dune::GmshGridFactory< Grid >;
  GridFactory gridFactory( (dim == 2) ? "grids/junction2d.msh" : "grids/plane3d.msh" );

  Grid& grid = *gridFactory.grid();
  const auto& igrid = grid.interfaceGrid();
  grid.loadBalance();
  const auto& gv = grid.leafGridView();

  static constexpr bool verbose = false;
  if constexpr (verbose)
  {
    std::cout << "Rank " << grid.comm().rank() << ": Elements " << grid.size(0) << " (" << grid.ghostSize(0) << " ghost)" << std::endl;
    std::cout << "Rank " << grid.comm().rank() << ": Facets " << grid.size(1) << " (" << grid.ghostSize(1) << " ghost)" << std::endl;
    std::cout << "Rank " << grid.comm().rank() << ": Vertices " << grid.size(dim) << " (" << grid.ghostSize(dim) << " ghost)" << std::endl;

    std::cout << "Interface Rank " << grid.comm().rank() << ": Elements " << igrid.size(0) << " (" << igrid.ghostSize(0) << " ghost)" << std::endl;
    std::cout << "Interface Rank " << grid.comm().rank() << ": Vertices " << igrid.size(dim-1) << " (" << igrid.ghostSize(dim-1) << " ghost)" << std::endl;
  }


  auto computeVolume = [&](const auto& grid, bool interface = false)
  {
    const auto& gv = grid.leafGridView();
    Dune::Timer timer;
    timer.start();

    double vol = 0.0;
    for (const auto& e : elements(gv, Dune::Partitions::interior))
      vol += e.geometry().volume();
    double sumvol = grid.comm().sum(vol);

    double sumexact = !interface ? 1.0 : ((dim == 2) ? 1.29155 : 0.191342);
    assert(std::abs(sumvol - sumexact) < 1e-5);

    auto dt = timer.elapsed();
    if (grid.comm().rank() == 0)
      std::cout << " " << (!interface ? "Bulk" : "Interface") << " took " << dt << std::endl;
  };

  // Compute volume
  if (grid.comm().rank() == 0)
    std::cout << "Compute volume:" << std::endl;
  computeVolume(grid);
  computeVolume(igrid, true);

  // Call gridcheck from dune-grid
  gridcheck( grid );

  // Call grid check for interface grid
  gridcheck( igrid );

  if constexpr (verbose)
  {
    auto getElementPartition = [](const auto& grid)
    {
      const auto& indexSet = grid.leafIndexSet();
      std::vector<int> partition (grid.size(0));
      for (const auto& e : elements(grid.leafGridView(), Partitions::all))
        partition[indexSet.index(e)] = e.partitionType();
      return partition;
    };

    auto getConnectivity = [](const auto& grid)
    {
      const auto& indexSet = grid.leafIndexSet();
      std::vector<int> connectivity (grid.size(0));
      for (const auto& e : elements(grid.leafGridView(), Partitions::all))
        connectivity[indexSet.index(e)] = (grid.partitionHelper().connectivity(e).size() > 0);
      return connectivity;
    };

    auto getVertexPartition = [](const auto& grid, int dim)
    {
      const auto& indexSet = grid.leafIndexSet();
      std::vector<int> vpartition (grid.size(dim));
      for (const auto& v : vertices(grid.leafGridView(), Partitions::all))
        vpartition[indexSet.index(v)] = v.partitionType();
      return vpartition;
    };

    VTKWriter vtkWriter( grid.leafGridView() );
    auto ep = getElementPartition(grid);
    vtkWriter.addCellData( ep, "partition" );
    auto hc = getConnectivity(grid);
    vtkWriter.addCellData( hc, "connectivity" );
    auto vp = getVertexPartition(grid, dim);
    vtkWriter.addVertexData(vp, "partition");
    vtkWriter.write("test-grid");

    VTKWriter ivtkWriter( igrid.leafGridView() );
    auto iep = getElementPartition(igrid);
    ivtkWriter.addCellData( iep, "partition" );
    auto ihc = getConnectivity(igrid);
    ivtkWriter.addCellData( ihc, "connectivity" );
    auto ivp = getVertexPartition(igrid, dim-1);
    ivtkWriter.addVertexData( ivp, "partition");
    ivtkWriter.write("test-igrid");
  }

  return EXIT_SUCCESS;
}
