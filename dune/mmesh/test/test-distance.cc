// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/float_cmp.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

//! Write and check vertex distances
// tol: Acceptable tolerance difference ratio
template<class Grid>
void writeAndCheckDistance(const Grid& grid, double tol = 1e-8)
{
  static constexpr int dim = Grid::dimension;

  Dune::VTKWriter<typename Grid::LeafGridView> vtkWriter( grid.leafGridView() );
  vtkWriter.addVertexData(grid.distance(), "distance");
  vtkWriter.write("distance"+std::to_string(dim));

  const auto& distance = grid.distance();
  for (const auto& vertex : vertices(grid.leafGridView()))
  {
    double exactDistance = std::abs(vertex.geometry().center()[1] - 0.5);
    auto reportedDistance = distance(vertex);

    double err = std::abs(reportedDistance - exactDistance) / exactDistance;
    if (err > tol)
      DUNE_THROW(InvalidStateException,
        "Distance of vertex at " << vertex.geometry().center() << " is " << exactDistance << " but was reported as " << reportedDistance
        << ". This is " << err * 100 << "% off!" << std::endl
      );
  }
}

/** Test-template main program. Instantiate a single expression
 * template and evaluate it on a simple grid.
 */
int main(int argc, char *argv[])
{
  MPIHelper::instance(argc, argv);
  std::cout << "-- Distance test --" << std::endl;

  using Grid2D = Dune::MovingMesh<2>;
  using GridFactory2D = Dune::GmshGridFactory< Grid2D >;

  // First test grid without interface
  GridFactory2D gridFactory2dsimple( "grids/simple2d.msh" );
  Grid2D& grid2dsimple = *gridFactory2dsimple.grid();

  // Then, actually check the computed distance
  GridFactory2D gridFactory2d( "grids/line2d.msh" );
  Grid2D& grid2d = *gridFactory2d.grid();
  writeAndCheckDistance(grid2d);

  using Grid3D = Dune::MovingMesh<3>;
  using GridFactory3D = Dune::GmshGridFactory< Grid3D >;
  GridFactory3D gridFactory3d( "grids/flat3d.msh" );
  Grid3D& grid3d = *gridFactory3d.grid();
  writeAndCheckDistance(grid3d, 0.05);

  return EXIT_SUCCESS;
}
