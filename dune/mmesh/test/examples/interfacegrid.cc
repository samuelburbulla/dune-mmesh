/** \example interfacegrid.cc
 * This is an example of how to build an MMesh using the .msh file format with a pre-defined interface
 * and obtain the interface grid after grid construction.
 *
 * A simple `horizontal2d.geo` file for a grid generation with interface could be:
 * \include grids/horizontal2d.geo
 * We can generate the .msh file using Gmsh by calling: ``gmsh -2 -format msh2 horizontal2d.geo``.
 *
 * Then, we can create an MMesh instance and obtain the interface grid with the following code.
 */
#include <config.h>

#include <dune/mmesh/mmesh.hh>

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  static constexpr int dim = GRIDDIM;

  using Grid = Dune::MovingMesh< dim >;

  using GridFactory = Dune::GmshGridFactory< Grid >;

  GridFactory gridFactory( "grids/horizontal" + std::to_string(dim) + "d.msh" );

  Grid& grid = *gridFactory.grid();

  std::cout << grid.size(0) << std::endl;

  using InterfaceGrid = typename Grid::InterfaceGrid;

  const InterfaceGrid& igrid = grid.interfaceGrid();

  std::cout << igrid.size(0) << std::endl;
}
