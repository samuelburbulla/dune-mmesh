/** \example gmsh-implicit.cc
 * This is an example of how to build an MMesh using the .msh file format with the implicit grid factory.
 * The implicit grid factory only inserts the vertices to the triangulation, the connectivity is performed implicitly by CGAL.
 * Here, we use the Dune::DelaunayTriangulation< dim > type to create a DelaunayTriangulation of the given vertices.
 *
 * We use the same mesh file as in gmsh.cc.
 *
 * We can create an MMesh instance implicitly with the following code.
 */
#include <config.h>
#include <dune/mmesh/mmesh.hh>

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  static constexpr int dim = GRIDDIM;

  using Grid = Dune::DelaunayTriangulation< dim >;

  using GridFactory = Dune::GmshGridFactory< Grid, /*implicit=*/true >;

  GridFactory gridFactory( "grids/cube" + std::to_string(dim) + "d.msh" );

  Grid& grid = *gridFactory.grid();

  std::cout << grid.size(0) << std::endl;
}
