// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/test/gridcheck.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

template< class GridView >
void checkTwists(const GridView& gridView)
{
  using Grid = typename GridView::Grid;
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

      std::cout << "In Inside" << std::endl;
      auto lgeoIn = is.geometryInInside();
      for ( int c = 0; c < dim; ++c )
      {
        auto x = geoIn.global( lgeoIn.corner(c) );
        std::cout << x << std::endl;

        auto refx = geoIn.corner( refIn.subEntity( nIn, 1, c, dim ) );
        std::cout << refx << std::endl;
        if( (x - refx).two_norm() > 1e-12 )
          DUNE_THROW(InvalidStateException, "refElem inside is different");
      }

      int tin = Fem::TwistUtility<Grid>::twistInSelf(gridView.grid(), is);
      std::cout << "Twist: " << tin << std::endl;

      std::cout << "Outside" << std::endl;
      const auto& facetOut = eOut.template subEntity<1>(nOut).geometry();
      for ( int c = 0; c < dim; ++c )
        std::cout << facetOut.corner(c) << std::endl;

      std::cout << "In Outside" << std::endl;
      auto lgeoOut = is.geometryInOutside();
      for ( int c = 0; c < dim; ++c )
      {
        auto x = geoOut.global( lgeoOut.corner(c) );
        std::cout << x << std::endl;

        auto refx = geoOut.corner( refOut.subEntity( nOut, 1, c, dim ) );
        if( (x - refx).two_norm() > 1e-12 )
          DUNE_THROW(InvalidStateException, "refElem outside is different");
      }

      int tout = Fem::TwistUtility<Grid>::twistInNeighbor(gridView.grid(), is);
      std::cout << "Twist: " << tout << std::endl;

      const auto& igeo = is.geometry();

      int t0 = (tin < 0) ? -(tin+1) : tin;

      auto x0 = geoIn.corner( refIn.subEntity( nIn, 1, 0, dim ) );
      std::cout << x0 << std::endl;
      std::cout << igeo.corner(t0) << std::endl;
      if( (x0 - igeo.corner(t0)).two_norm() > 1e-12 )
        DUNE_THROW(InvalidStateException, "twist in false");

      if constexpr (dim == 3)
      {
        int t1 = ((tin < 0) ? -(tin+1)+(dim-1) : (tin+1))%dim;
        auto x1 = geoIn.corner( refIn.subEntity( nIn, 1, 1, dim ) );
        if( (x1 - igeo.corner(t1)).two_norm() > 1e-12 )
          DUNE_THROW(InvalidStateException, "twist sign in false");
      }

      t0 = (tout < 0) ? -(tout+1) : tout;

      auto c0 = geoOut.corner( refOut.subEntity( nOut, 1, 0, dim ) );
      if( (c0 - igeo.corner(t0)).two_norm() > 1e-12 )
        DUNE_THROW(InvalidStateException, "twist out false");

      if constexpr (dim == 3)
      {
        int t1 = ((tout < 0) ? -(tout+1)+(dim-1) : (tout+1))%dim;
        auto c1 = geoOut.corner( refOut.subEntity( nOut, 1, 1, dim ) );
        if( (c1 - igeo.corner(t1)).two_norm() > 1e-12 )
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
  try {
    std::cout << "-- Reference element check --" << std::endl;

    // Create Grid
    // ------------
    static constexpr int dim = GRIDDIM;
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
