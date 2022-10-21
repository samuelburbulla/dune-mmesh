#include <dune/mmesh/mmesh.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <random>
#include <algorithm>
#include <dune/mmesh/grid/polygoncutting.hh>
#include <dune/common/fvector.hh>


int main(int argc, char** argv)
{
  try
  {
    MPIHelper::instance(argc, argv);
    //dimension of the grid
    static constexpr int dim = GRIDDIM;

    using Point = Dune::FieldVector<double, dim>;
    using PC = Dune::PolygonCutting<double, Point>;
    using PointList = std::vector<Point>;

    //create random number generator
    std::random_device rd;  //used to obtain a seed for the random number engine
    std::mt19937 rng(rd()); //mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> unif(0.0, 1.0);

    PointList polygon1;
    PointList polygon2;
    PointList polygonAfterCut;

    const int numberOfTests = 1e3;
    int numberOfSuccessfulTests = 0;

    for (; numberOfSuccessfulTests < numberOfTests; numberOfSuccessfulTests++)
    {
      polygon1 = PointList({Point({unif(rng), unif(rng)}),
        Point({unif(rng), unif(rng)}), Point({unif(rng), unif(rng)})});
      if (numberOfSuccessfulTests % 2 == 0)
      {
        polygon2 = PointList({Point({unif(rng), unif(rng)}),
          Point({unif(rng), unif(rng)}), Point({unif(rng), unif(rng)})});
      }
      else
      {
        polygon2 = PointList({polygon1[0], polygon1[1],
          Point({unif(rng), unif(rng)})});
      }

      //ensure COUNTERclockwise ordering (!!!)
      if (PC::polygonArea(polygon1) < 0)
      {
        std::reverse(polygon1.begin(), polygon1.end());
      }
      if (PC::polygonArea(polygon2) < 0)
      {
        std::reverse(polygon2.begin(), polygon2.end());
      }

      polygonAfterCut = PC::sutherlandHodgman(polygon1, polygon2);

      //compare intersection volume with CGAL algorithm
      const double intersectionArea = PC::polygonArea(polygonAfterCut);
      const double intersectionAreaCGAL =
        PC::intersectionVolumeCGAL(polygon1, polygon2);

      if (std::abs(intersectionArea - intersectionAreaCGAL) > 1e-10)
        DUNE_THROW( Dune::InvalidStateException, "The intersection volumes differ!" );
    }

    return 0;

  }
  catch (Dune::Exception& e)
  {
    std::cerr << "Dune reported error: " << e.what() << std::endl;
  }
  catch (std::exception& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
  }
  catch (...)
  {
    std::cerr << "Unknown exception thrown!" << std::endl;
  }
}
