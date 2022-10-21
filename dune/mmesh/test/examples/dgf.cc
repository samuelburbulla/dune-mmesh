/** \example dgf.cc
 * This is an example of how to build an MMesh using a .dgf file.
 */
#include <config.h>
#include <dune/mmesh/mmesh.hh>

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  static constexpr int dim = 2;

  using Grid = Dune::MovingMesh< dim >;

  using GridFactory = Dune::DGFGridFactory< Grid >;

  GridFactory gridFactory( "grids/cube" + std::to_string(dim) + "d.dgf" );

  Grid& grid = *gridFactory.grid();

  std::cout << grid.size(0) << std::endl;
}
