/** \example gmsh.cc
 * This is an example of how to build an MMesh using the .msh file format.
 *
 * A simple `cube2d.geo` file for grid generation could be:
 * \include grids/cube2d.geo
 * We can generate the .msh file using Gmsh by calling: ``gmsh -2 -format msh2 cube2d.geo``.
 * For 3D use: ``gmsh -3 -format msh2 cube3d.geo``.
 *
 * Then, we can create an MMesh instance with the following code.
 */
#include <config.h>
#include <dune/mmesh/mmesh.hh>

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  static constexpr int dim = GRIDDIM;

  using Grid = Dune::MovingMesh< dim >;

  using GridFactory = Dune::GmshGridFactory< Grid >;

  GridFactory gridFactory( "grids/cube" + std::to_string(dim) + "d.msh" );

  Grid& grid = *gridFactory.grid();

  std::cout << grid.size(0) << std::endl;
}
