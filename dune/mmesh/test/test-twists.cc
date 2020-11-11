// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

template< class GridView >
void checkTwists(const GridView& gridView)
{
  static constexpr int dim = GridView::dimension;

  for ( const auto& e : elements(gridView) )
  {
    for ( const auto& is : intersections(gridView, e) )
    {
      if (is.boundary())
        continue;

      auto eIn = is.inside();
      auto nIn = is.indexInInside();
      auto geoIn = eIn.geometry();

      auto eOut = is.outside();
      auto nOut = is.indexInOutside();
      auto geoOut = eOut.geometry();

      auto refIn = referenceElement( geoIn );
      auto refOut = referenceElement( geoOut );

      std::cout << "Inside" << std::endl;
      const auto& facetIn = eIn.template subEntity<1>(nIn).geometry();
      for ( int c = 0; c < dim; ++c )
        std::cout << facetIn.corner(c) << std::endl;

      const auto& isgeo = is.geometry();

      auto x0 = isgeo.corner( 0 );
      auto y0 = geoIn.corner( refIn.subEntity( nIn, 1, 0, dim ) );
      std::cout << x0 << "   -   " << y0 << std::endl;
      if( (x0 - y0).two_norm() > 1e-12 )
        DUNE_THROW(InvalidStateException, "twist in false");

      auto z0 = geoIn.global( is.geometryInInside().corner(0) );
      std::cout << x0 << "   -   " << z0 << std::endl;
      if( (x0 - z0).two_norm() > 1e-12 )
        DUNE_THROW(InvalidStateException, "geoInInside false");

      if constexpr (dim == 3)
      {
        auto x1 = isgeo.corner( 1 );
        auto y1 = geoIn.corner( refIn.subEntity( nIn, 1, 1, dim ) );
        std::cout << x1 << "   -   " << y1 <<std::endl;
        if( (x1 - y1).two_norm() > 1e-12 )
          DUNE_THROW(InvalidStateException, "twist sign in false");
      }


      std::cout << std::endl << "Outside" << std::endl;
      const auto& facetOut = eOut.template subEntity<1>(nOut).geometry();
      for ( int c = 0; c < dim; ++c )
        std::cout << facetOut.corner(c) << std::endl;

      auto c0 = isgeo.corner( 0 );
      auto d0 = geoOut.corner( refOut.subEntity( nOut, 1, 0, dim ) );
      std::cout << c0 << "   -   " << d0 << std::endl;
      if( (c0 - d0).two_norm() > 1e-12 )
        DUNE_THROW(InvalidStateException, "twist out false");

      auto e0 = geoOut.global( is.geometryInOutside().corner(0) );
      std::cout << c0 << "   -   " << e0 << std::endl;
      if( (c0 - e0).two_norm() > 1e-12 )
        DUNE_THROW(InvalidStateException, "geoInOutside false");

      if constexpr (dim == 3)
      {
        auto c1 = isgeo.corner( 1 );
        auto d1 = geoOut.corner( refOut.subEntity( nOut, 1, 1, dim ) );
        std::cout << c1 << "   -   " << d1 << std::endl;
        if( (c1 - d1).two_norm() > 1e-12 )
          DUNE_THROW(InvalidStateException, "twist sign out false");
      }

      std::cout << std::endl;
    }
  }
}
/** Test-template main program. Instantiate a single expression
 * template and evaluate it on a simple grid.
 */
int main(int argc, char *argv[])
{
  MPIHelper::instance( argc, argv );

  try {
    std::cout << "-- Reference element check --" << std::endl;
    static constexpr int dim = GRIDDIM;

    // Create Grid
    // ------------

    using Grid = Dune::MovingMesh<dim>;
    using GridFactory = Dune::GmshGridFactory< Grid >;
    GridFactory gridFactory( (dim == 2) ? "grids/line2d.msh" : "grids/plane3d.msh" );
    Grid& grid = *gridFactory.grid();
    const auto& gridView = grid.leafGridView();

    checkTwists(gridView);

    std::cout << std::endl << " = Check interface grid =" << std::endl << std::endl;

    const auto& igridView = grid.interfaceGrid().leafGridView();
    checkTwists(igridView);

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
