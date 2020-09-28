// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/test/gridcheck.hh>
#include <dune/grid/test/checktwists.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

/** Test-template main program. Instantiate a single expression
 * template and evaluate it on a simple grid.
 */
int main(int argc, char *argv[])
{
  std::cout << "-- Grid check --" << std::endl;

  // Create Grid
  // ------------
  static constexpr int dim = GRIDDIM;
  using Grid = Dune::MovingMesh<dim>;

  using GridFactory = Dune::GmshGridFactory< Grid >;
  GridFactory gridFactory( (dim == 2) ? "grids/line2d.msh" : "grids/plane3d.msh" );

  Grid& grid = *gridFactory.grid();

  // Call gridcheck from dune-grid
  gridcheck( grid );

  // Call grid check for interface grid
  gridcheck( grid.interfaceGrid() );

  return EXIT_SUCCESS;
}
