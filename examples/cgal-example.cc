// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <CGAL/exceptions.h>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

/** Test-template main program. Instantiate a single expression
 * template and evaluate it on a simple grid.
 */
int main(int argc, char *argv[])
{
  try {
    // Use CGAL to create a triangulation
    // ==================================

    typedef MMeshDefaults::Delaunay<2>::Triangulation Triangulation;
    typedef Triangulation::Point Point;

    std::vector<Point> points;
    points.push_back(Point(0,0));
    points.push_back(Point(1,0));
    points.push_back(Point(0,1));
    points.push_back(Point(1,1));

    Triangulation tr;
    tr.insert(points.begin(), points.end());


    // Create Dune grid  wrapper
    // =========================
    MMesh<Triangulation, 2> grid(tr);


    // Test MMesh wrapper
    auto gridView = grid.leafGridView();
    auto indexSet = grid.leafIndexSet();

    std::cout << "Number of Cells: " << indexSet.size(0) << std::endl;

    for(auto e : elements(gridView)) {
      std::cout << "Entity Center: " << e.geometry().center() << std::endl;
      std::cout << "Volume of Element: " << e.geometry().volume() << std::endl;

      std::cout << "Index: " << indexSet.index( e ) << std::endl;

      const int codim = 2;
      for ( std::size_t i = 0; i < e.subEntities( codim ); ++i ) {
        std::cout << " Codim " << codim << " subindex " << i << ": " << indexSet.subIndex( e, i, codim );
        auto subEntity = e.subEntity<codim>( i );
        std::cout << " at " << subEntity.geometry().center() << std::endl;
      }

      for(auto is : intersections(gridView, e)) {
        std::cout << "  Intersection CenterUnitOuterNormal: " << is.centerUnitOuterNormal() << std::endl;
      }

    }

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
