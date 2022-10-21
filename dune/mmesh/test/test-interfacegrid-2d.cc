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
#include <dune/grid/io/file/vtk/vtkwriter.hh>

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
  return std::abs(a - b) < 1e-10;
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
    std::cout << "-- MMesh interface implementation test --" << std::endl;

    std::cout << "Build two cell grid..." << std::endl;

    // Create MMesh
    // ------------
    static constexpr int dim = GRIDDIM;
    using Grid = Dune::MovingMesh<dim>;

    using GridFactory = Dune::GmshGridFactory< Grid >;
    GridFactory gridFactory( "grids/mimesh" + std::to_string(dim) + "d.msh" );

    Grid& mMesh = *gridFactory.grid();

    // ---------------------------------------

    using InterfaceGrid = typename Grid::InterfaceGrid;
    const InterfaceGrid& igrid = mMesh.interfaceGrid();

    // ---------------------------------------

    std::cout << "Build LeafGridView and LeafIndexSet..." << std::endl;
    auto gridView = igrid.leafGridView();
    auto indexSet = igrid.leafIndexSet();

    // Write grid
    VTKWriter<typename InterfaceGrid::LeafGridView> vtkWriter( gridView );
    vtkWriter.write("test-interfacegrid-2d-0");

    // Test MMesh index set
    checkProperty( "size of element index set", indexSet.size(0), 6ul );
    checkProperty( "size of vertex index set", indexSet.size(1), 7ul );

    // Test MMesh mcmgmapper for elements
    MultipleCodimMultipleGeomTypeMapper< decltype( gridView ) > elementMapper ( gridView, mcmgElementLayout() );
    checkProperty( "size of mcmg element mapper", elementMapper.size(), 6ul );
    unsigned int count = 0;
    for(auto e : elements(gridView))
      checkProperty( "index mapped by element mapper", elementMapper.index(e), count++ );

    // Test MMesh mcmgmapper for vertices
    MultipleCodimMultipleGeomTypeMapper< decltype( gridView ) > vertexMapper ( gridView, mcmgVertexLayout() );
    checkProperty( "size of mcmg vertex mapper", vertexMapper.size(), 7ul );
    unsigned int vertexCount = 0;
    std::vector< unsigned int > idxMap {{ 6, 0, 2, 5, 4, 1, 3 }};
    for(auto v : vertices(gridView))
      checkProperty( "index mapped by vertex mapper", vertexMapper.index(v), idxMap[ vertexCount++ ] );

    // Iterate over interface
    int numberOfInterfaceElements = 0;
    for ( const auto& segment : interfaceElements( mMesh.leafGridView() ) )
    {
      checkProperty( "isInterface", mMesh.isInterface( segment ), true );
      numberOfInterfaceElements++;
    }
    checkProperty( "number of interface elements", numberOfInterfaceElements, 6 );

    int elementCount = 0;
    for(auto e : elements(gridView))
    {
      const auto geo = e.geometry();

      if (elementCount == 1)
      {
        std::cout << "- Check second element -" << std::endl;

        // check element
        checkProperty( "geometry center", geo.center(), { 0.375, 0.5 } );
        checkProperty( "geometry volume", geo.volume(), 0.25 );
        checkProperty( "element index", indexSet.index( e ), 1ul );

        // check vertices
        checkProperties( "vertex sub indices",
          { { indexSet.subIndex( e, 0, 1 ), 2ul },
            { indexSet.subIndex( e, 1, 1 ), 1ul } }
        );

        checkProperties( "vertex positions",
          { { e.subEntity<1>( 0 ).geometry().center(), { 0.5, 0.5 } },
            { e.subEntity<1>( 1 ).geometry().center(), { 0.25000000000103079,  0.5 } } }
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
            { is1.boundary(), false } }
        );

        checkProperties( "intersection normals",
          { { is0.centerUnitOuterNormal(), {  1.0, 0.0 } },
            { is1.centerUnitOuterNormal(), {  1.0, 0.0 } },
            { is2.centerUnitOuterNormal(), { -1.0, 0.0 } } }
        );

        checkProperty( "index in inside", is1.indexInInside(), 0 );
        checkProperty( "index in outside", is1.indexInOutside(), 0 );

        const auto& neighborA = is0.impl().outside();
        const auto& neighborB = is1.impl().outside();
        const auto& neighborC = is2.impl().outside();

        // check global vertex indices obtained by intersection subentities
        const auto refElement = Dune::ReferenceElements<double, 1>::general(geo.type());

        const auto vIdxLocal1 = refElement.subEntity(0, 0, 0, 1);
        const auto vIdxGlobal1 = vertexMapper.subIndex(e, vIdxLocal1, 1);
        const auto vIdxLocal2 = refElement.subEntity(0, 0, 1, 1);
        const auto vIdxGlobal2 = vertexMapper.subIndex(e, vIdxLocal2, 1);

        checkProperties( "reference element mapping",
          { { vIdxGlobal1, 2u },
            { vIdxGlobal2, 1u } }
        );

        // try to obtain mmesh intersection
        const auto& bulkIs = mMesh.asIntersection( e );
        checkProperty( "asIntersection", bulkIs.geometry().center(), e.geometry().center() );
      }
      elementCount++;
    }

    int vi = 0;
    std::vector<int> numbOfInc = {{ 1, 1, 3, 1, 2, 2, 2 }};
    std::vector<int> numbOfIncV = {{ 1, 1, 3, 1, 3, 2, 3 }};
    for(auto v : vertices(gridView))
    {
      int c = 0;
      for(auto e : incidentInterfaceElements(v))
      {
        e.geometry();
        c++;
      }
      checkProperty( "number of incident elements", c, numbOfInc[vi] );

      c = 0;
      for(auto e : incidentInterfaceVertices(v))
      {
        e.geometry();
        c++;
      }
      checkProperty( "number of incident vertices", c, numbOfIncV[vi] );

      vi++;
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
