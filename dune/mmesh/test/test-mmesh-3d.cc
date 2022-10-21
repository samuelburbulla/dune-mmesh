// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/mmesh/mmesh.hh>

using namespace Dune;

template<class T>
bool compare( const T& a, const T& b )
{
  return a == b;
}

template<>
bool compare( const double& a, const double& b )
{
  return std::abs(a - b) < std::numeric_limits<double>::epsilon();
}

template<>
bool compare( const Dune::FieldVector<double,3>& a, const Dune::FieldVector<double,3>& b )
{
  return (a - b).two_norm() < 1e-12;
}

//! Helper method for checking multiple properties given in a list with verbose output
template<class T>
void checkProperties( const std::string &name, std::initializer_list< std::initializer_list<T> > data )
{
  std::cout << "Check " << std::setfill(' ') << std::setw(30) << std::left << name;

  for( std::vector<T> datavector : data )
  {
    const T& resultValue = datavector[0];
    const T& testValue = datavector[1];

    bool correct = compare( resultValue, testValue );

    if ( !correct )
    {
      std::cout.precision(17);
      std::cout << " - failed: result [" << resultValue << "] was expected to be [" << testValue << "]" << std::endl;
    }
    assert( correct );
  }

  std::cout << " - correct" << std::endl;
}

//! Helper method for checking a single property
template<class T>
void checkProperty( const std::string &name, const T& resultValue, const T& testValue )
{
  return checkProperties( name, {{ resultValue, testValue }} );
}

/** Test-template main program. Instantiate a single expression
 * template and evaluate it on a simple grid.
 */
int main(int argc, char *argv[])
{
  try {
    MPIHelper::instance(argc, argv);
    std::cout << "-- MMesh implementation test for 3D --" << std::endl;

    std::cout << "Build simple grid..." << std::endl;

    // Create MMesh
    // ------------
    using Grid = Dune::MovingMesh<3>;

    using GridFactory = Dune::GmshGridFactory< Grid >;
    GridFactory gridFactory( "grids/twocell3d.msh" );

    Grid& mMesh = *gridFactory.grid();

    // ---------------------------------------

    std::cout << "Build LeafGridView and LeafIndexSet..." << std::endl;
    auto gridView = mMesh.leafGridView();
    auto indexSet = mMesh.leafIndexSet();

    // Test MMesh index set
    checkProperty( "size of element index set", indexSet.size(0), 2ul );
    checkProperty( "size of face index set", indexSet.size(1), 7ul );
    checkProperty( "size of edge index set", indexSet.size(2), 9ul );
    checkProperty( "size of vertex index set", indexSet.size(3), 5ul );

    // Test iterators
    std::size_t countCells = 0;
    for( const auto& cell : elements(gridView) )
    {
        cell.geometry();
        countCells++;
    }
    checkProperty( "number of iterated elements", countCells, indexSet.size(0) );

    std::size_t countFacets = 0;
    for( const auto& facet : facets(gridView) )
    {
        facet.geometry();
        countFacets++;
    }
    checkProperty( "number of iterated facets", countFacets, indexSet.size(1) );

    std::size_t countEdges = 0;
    for( const auto& edge : edges(gridView) )
    {
        edge.geometry();
        countEdges++;
    }
    checkProperty( "number of iterated edges", countEdges, indexSet.size(2) );

    std::size_t countVertices = 0;
    for( const auto& vertex : vertices(gridView) )
    {
        vertex.geometry().center();
        countVertices++;
    }
    checkProperty( "number of iterated vertices", countVertices, indexSet.size(3) );


    // Test MMesh mcmgmapper
    MultipleCodimMultipleGeomTypeMapper< decltype( gridView ) > vertexMapper ( gridView, mcmgVertexLayout() );
    checkProperty( "size of mcmg vertex mapper", vertexMapper.size(), 5ul );

    int elementCount = 0;
    for(auto e : elements(gridView))
    {
      const auto geo = e.geometry();

      if ( elementCount == 0)
      {
        std::cout << "- Check first element -" << std::endl;

        // check element
        checkProperty( "geometry center", geo.center(), { 0.375, 0.375, 0.25 } );
        checkProperty( "geometry volume", geo.volume(), 1.0/6.0 );
        checkProperty( "element index", indexSet.index( e ), 0ul );

        // check vertices
        checkProperties( "vertex sub indices",
          { { indexSet.subIndex( e, 0, 3 ), 0ul },
            { indexSet.subIndex( e, 1, 3 ), 1ul },
            { indexSet.subIndex( e, 2, 3 ), 3ul },
            { indexSet.subIndex( e, 3, 3 ), 4ul } }
        );
        checkProperties( "vertex positions",
          { { e.subEntity<3>( 0 ).geometry().center(), { 0, 0, 0 } },
            { e.subEntity<3>( 1 ).geometry().center(), { 1, 0, 0 } },
            { e.subEntity<3>( 2 ).geometry().center(), { 0, 1, 0 } },
            { e.subEntity<3>( 3 ).geometry().center(), { 0.5, 0.5, 1 } } }
        );

        // check faces
        checkProperties( "face centers",
          { { e.subEntity<1>( 0 ).geometry().center(), { 0.33333333333333331, 0.33333333333333331, 0 } },
            { e.subEntity<1>( 1 ).geometry().center(), { 0.5, 0.16666666666666671, 0.33333333333333343 } },
            { e.subEntity<1>( 2 ).geometry().center(), { 0.16666666666666666, 0.5, 0.33333333333333331 } },
            { e.subEntity<1>( 3 ).geometry().center(), { 0.50000000000000011, 0.5, 0.33333333333333331 } } }
        );

        // check intersections
        auto it = e.impl().ileafbegin();
        auto is0 = it.dereference();
        it.increment();
        auto is1 = it.dereference();
        it.increment();
        auto is2 = it.dereference();
        it.increment();
        auto is3 = it.dereference();

        checkProperties( "intersection normals",
          { { is0.centerUnitOuterNormal(), { 0, 0, -1 } },
            { is1.centerUnitOuterNormal(), { 0, -0.89442719099991586, 0.44721359549995793 } },
            { is2.centerUnitOuterNormal(), { -0.89442719099991586, 0, 0.44721359549995793 } },
            { is3.centerUnitOuterNormal(), { 0.70710678118654746, 0.70710678118654746, 0 } } }
        );

        checkProperties( "index in inside and outside",
          { { is1.indexInInside(), 1 },
            { is1.indexInOutside(), 0 } }
        );
      }

      if ( elementCount == 1)
      {
        std::cout << "- Check second element -" << std::endl;

        // check element
        checkProperty( "geometry center", geo.center(), { 0.625, 0.625, 0.25 } );
        checkProperty( "geometry volume", geo.volume(), 0.16666666666666666 );
        checkProperty( "element index", indexSet.index( e ), 1ul );

        // check intersections
        auto it = e.impl().ileafbegin();
        it.increment();
        it.increment();
        auto is2 = it.dereference();

        checkProperties( "index in inside and outside",
          { { is2.indexInInside(), 2 },
            { is2.indexInOutside(), 3 } }
        );
      }

      // Check entity seeds
      auto seed = e.seed();
      decltype(seed) seed2;
      seed2 = seed;
      const auto entity = mMesh.entity( seed2 );
      entity.geometry(); // dummy
      assert( e == entity );

      elementCount++;
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
