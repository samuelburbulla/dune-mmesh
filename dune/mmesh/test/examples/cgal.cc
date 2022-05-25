/** \example cgal.cc
 * This is an example of how to build an MMesh using a CGAL triangulation.
 */
#include <config.h>
#include <dune/mmesh/mmesh.hh>

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);
  using Triangulation = Dune::MMeshDefaults::Triangulation<2>::type;
  // alternative: using Triangulation = Dune::MMeshDefaults::Delaunay<2>::type;

  using Point = Triangulation::Point;

  std::vector<Point> points;
  points.push_back(Point(0,0));
  points.push_back(Point(1,0));
  points.push_back(Point(0,1));
  points.push_back(Point(1,1));

  Triangulation tr;
  tr.insert(points.begin(), points.end());

  using Grid = Dune::MMesh<Triangulation, 2>;
  Grid grid(tr);

  std::cout << grid.size(0) << std::endl;
}
