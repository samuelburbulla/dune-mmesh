// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/geometry/referenceelements.hh>

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
bool compare( const Dune::FieldVector<double,2>& a, const Dune::FieldVector<double,2>& b )
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
    std::cout << "-- MMesh implementation test --" << std::endl;

    std::cout << "Build simple grid..." << std::endl;

    // Create MMesh
    // ------------
    using Grid = Dune::MovingMesh<2>;

    using GridFactory = Dune::GmshGridFactory< Grid >;
    GridFactory gridFactory( "grids/simple2d.msh" );

    Grid& mMesh = *gridFactory.grid();

    // ---------------------------------------

    std::cout << "Build LeafGridView and LeafIndexSet..." << std::endl;
    auto gridView = mMesh.leafGridView();
    auto indexSet = mMesh.leafIndexSet();

    // Test MMesh index set
    checkProperty( "size of element index set", indexSet.size(0), 4ul );
    checkProperty( "size of edge index set", indexSet.size(1), 9ul );
    checkProperty( "size of vertex index set", indexSet.size(2), 6ul );

    // Test iterators
    std::size_t countCells = 0;
    for( const auto& cell : elements(gridView) )
    {
        cell.geometry();
        countCells++;
    }
    checkProperty( "number of iterated elements", countCells, indexSet.size(0) );

    std::size_t countEdges = 0;
    for( const auto& edge : edges(gridView) )
    {
        edge.geometry();
        countEdges++;
    }
    checkProperty( "number of iterated edges", countEdges, indexSet.size(1) );

    std::size_t countVertices = 0;
    for( const auto& vertex : vertices(gridView) )
    {
        vertex.geometry();
        countVertices++;
    }
    checkProperty( "number of iterated vertices", countVertices, indexSet.size(2) );

    // Test MMesh mcmgmapper for elements
    MultipleCodimMultipleGeomTypeMapper< decltype( gridView ) > elementMapper ( gridView, mcmgElementLayout() );
    checkProperty( "size of mcmg element mapper", elementMapper.size(), 4ul );
    unsigned int count = 0;
    for(auto e : elements(gridView))
      checkProperty( "index mapped by element mapper", elementMapper.index(e), count++ );

    // Test MMesh mcmgmapper for vertices
    MultipleCodimMultipleGeomTypeMapper< decltype( gridView ) > vertexMapper ( gridView, mcmgVertexLayout() );
    checkProperty( "size of mcmg vertex mapper", vertexMapper.size(), 6ul );
    unsigned int vertexCount = 0;
    for(auto v : vertices(gridView))
      checkProperty( "index mapped by vertex mapper", vertexMapper.index(v), vertexCount++ );

    int elementCount = 0;
    for(auto e : elements(gridView))
    {
      const auto geo = e.geometry();

      if (elementCount == 2)
      {
        std::cout << "- Check third element -" << std::endl;

        // check element
        checkProperty( "geometry center", geo.center(), { 0.66666666666666663, 0.083333333333333329 } );
        checkProperty( "geometry volume", geo.volume(), 0.125 );
        checkProperty( "element index", indexSet.index( e ), 2ul );

        // check vertices
        checkProperties( "vertex sub indices",
          { { indexSet.subIndex( e, 0, 2 ), 0ul },
            { indexSet.subIndex( e, 1, 2 ), 1ul },
            { indexSet.subIndex( e, 2, 2 ), 4ul } }
        );

        checkProperties( "vertex positions",
          { { e.subEntity<2>( 0 ).geometry().center(), { 0, 0 } },
            { e.subEntity<2>( 1 ).geometry().center(), { 1, 0 } },
            { e.subEntity<2>( 2 ).geometry().center(), { 1, 0.25 } } }
        );

        // check edges
        checkProperties( "edge centers",
          { { e.subEntity<1>( 0 ).geometry().center(), { 0.5, 0 } },
            { e.subEntity<1>( 1 ).geometry().center(), { 0.5, 0.125 } },
            { e.subEntity<1>( 2 ).geometry().center(), { 1, 0.125 } } }
        );

        // check intersections
        auto it = e.impl().ileafbegin();
        auto is0 = it.dereference();
        it.increment();
        auto is1 = it.dereference();
        it.increment();
        auto is2 = it.dereference();

        checkProperties( "boundary",
          { { is0.boundary(), true },
            { is1.boundary(), false },
            { is2.boundary(), true } }
        );

        checkProperties( "intersection normals",
          { { is0.centerUnitOuterNormal(), { 0, -1 } },
            { is1.centerUnitOuterNormal(), { -0.24253562503633297, 0.97014250014533188} },
            { is2.centerUnitOuterNormal(), { 1, 0 } } }
        );

        checkProperties( "index in inside and outside",
          { { is1.indexInInside(), 1 },
            { is1.indexInOutside(), 0 } }
        );

        // check global vertex indices obtained by intersection subentities
        const auto refElement = Dune::ReferenceElements<double, 2>::general(geo.type());

        const auto vIdxLocal1 = refElement.subEntity(0, 1, 0, 2);
        const auto vIdxGlobal1 = vertexMapper.subIndex(e, vIdxLocal1, 2);
        const auto vIdxLocal2 = refElement.subEntity(0, 1, 1, 2);
        const auto vIdxGlobal2 = vertexMapper.subIndex(e, vIdxLocal2, 2);

        checkProperties( "reference element mapping",
          { { vIdxGlobal1, 0u },
            { vIdxGlobal2, 1u } }
        );
      }

      if (elementCount == 1)
      {
        std::cout << "- Check second element -" << std::endl;

        // check element
        checkProperty( "geometry center", geo.center(), { 0.5, 0.41666666666666663 } );
        checkProperty( "geometry volume", geo.volume(), 0.4375 );
        checkProperty( "element index", indexSet.index( e ), 1ul );

        // check vertices
        checkProperties( "vertex sub indices",
          { { indexSet.subIndex( e, 0, 2 ), 0ul },
            { indexSet.subIndex( e, 1, 2 ), 4ul },
            { indexSet.subIndex( e, 2, 2 ), 5ul } }
        );
        checkProperties( "vertex positions",
          { { e.subEntity<2>( 0 ).geometry().center(), { 0, 0 } },
            { e.subEntity<2>( 1 ).geometry().center(), { 1, 0.25 } },
            { e.subEntity<2>( 2 ).geometry().center(), { 0.5, 1 } } }
        );

        // check edges
        checkProperties( "edge centers",
          { { e.subEntity<1>( 0 ).geometry().center(), { 0.5, 0.125  } },
            { e.subEntity<1>( 1 ).geometry().center(), { 0.25, 0.5 } },
            { e.subEntity<1>( 2 ).geometry().center(), { 0.75, 0.625 } } }
        );

        // check intersections
        auto it = e.impl().ileafbegin();
        auto is0 = it.dereference();
        it.increment();
        auto is1 = it.dereference();
        it.increment();
        auto is2 = it.dereference();

        checkProperties( "boundary",
          { { is0.boundary(), false },
            { is1.boundary(), false },
            { is2.boundary(), false } }
        );

        checkProperties( "isIntersection",
          { { mMesh.isInterface( is0 ), false },
            { mMesh.isInterface( is1 ), false },
            { mMesh.isInterface( is2 ), false } }
        );

        checkProperties( "intersection normals",
          { { is0.centerUnitOuterNormal(), { 0.24253562503633297, -0.97014250014533188 } },
            { is1.centerUnitOuterNormal(), { -0.89442719099991586, 0.44721359549995793 } },
            { is2.centerUnitOuterNormal(), { 0.83205029433784372, 0.55470019622522915 } } }
        );

        checkProperties( "index in inside and outside",
          { { is2.indexInInside(), 2 },
            { is2.indexInOutside(), 2 } }
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
