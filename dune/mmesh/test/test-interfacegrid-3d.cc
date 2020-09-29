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
    std::cout << "-- MMesh interface implementation test --" << std::endl;

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
    vtkWriter.write("test-interfacegrid-3d-0");

    // Test MMesh index set
    checkProperty( "size of element index set", indexSet.size(0), 4ul );
    checkProperty( "size of vertex index set", indexSet.size(2), 5ul );

    // Test MMesh mcmgmapper for elements
    MultipleCodimMultipleGeomTypeMapper< decltype( gridView ) > elementMapper ( gridView, mcmgElementLayout() );
    checkProperty( "size of mcmg element mapper", elementMapper.size(), 4u );
    unsigned int count = 0;
    for(auto e : elements(gridView))
      checkProperty( "index mapped by element mapper", elementMapper.index(e), count++ );

    // Test MMesh mcmgmapper for vertices
    MultipleCodimMultipleGeomTypeMapper< decltype( gridView ) > vertexMapper ( gridView, mcmgVertexLayout() );
    checkProperty( "size of mcmg vertex mapper", vertexMapper.size(), 5u );
    unsigned int vertexCount = 0;
    std::vector< unsigned int > idxMap {{ 4, 1, 0, 3, 2 }};
    for(auto v : vertices(gridView))
      checkProperty( "index mapped by vertex mapper", vertexMapper.index(v), idxMap[vertexCount++] );

    // Iterate over interface
    int numberOfInterfaceElements = 0;
    for ( const auto& segment : interfaceElements( mMesh.leafGridView() ) )
    {
      checkProperty( "isInterface", mMesh.isInterface( segment ), true );
      numberOfInterfaceElements++;
    }
    checkProperty( "number of interface elements", numberOfInterfaceElements, 4 );

    int elementCount = 0;
    for(auto e : elements(gridView))
    {
      const auto geo = e.geometry();

      if (elementCount == 0)
      {
        std::cout << "- Check first element -" << std::endl;

        // check element
        checkProperty( "geometry center", geo.center(), { 0.83333333333333337, 0.5, 0.5 } );
        checkProperty( "geometry volume", geo.volume(), 0.25 );
        checkProperty( "element index", indexSet.index( e ), 0ul );

        // check vertices
        checkProperties( "vertex sub indices",
          { { indexSet.subIndex( e, 0, 2 ), 1ul },
            { indexSet.subIndex( e, 1, 2 ), 0ul },
            { indexSet.subIndex( e, 2, 2 ), 2ul } }
        );

        checkProperties( "vertex positions",
          { { e.subEntity<2>( 0 ).geometry().center(), { 1.0, 0.5, 0.0 } },
            { e.subEntity<2>( 1 ).geometry().center(), { 1.0, 0.5, 1.0 } },
            { e.subEntity<2>( 2 ).geometry().center(), { 0.5, 0.5, 0.5 } } }
        );

        // check intersections
        auto it = e.impl().ileafbegin();
        auto is0 = it.dereference();
        it.increment();
        auto is1 = it.dereference();
        it.increment();
        auto is2 = it.dereference();

        checkProperties( "boundary",
          { { is0.boundary(), true  },
            { is1.boundary(), false },
            { is2.boundary(), false } }
        );

        checkProperties( "intersection centers",
          { { is0.geometry().center(), {  1.0, 0.5,  0.5 } },
            { is1.geometry().center(), { 0.75, 0.5, 0.25 } },
            { is2.geometry().center(), { 0.75, 0.5, 0.75 } } }
        );

        checkProperties( "intersection normals",
          { { is0.centerUnitOuterNormal(), { 1.0, 0.0, 0.0 } },
            { is1.centerUnitOuterNormal(), { -1.0 / std::sqrt(2), 0.0, -1.0 / std::sqrt(2) } },
            { is2.centerUnitOuterNormal(), { -1.0 / std::sqrt(2), 0.0,  1.0 / std::sqrt(2) } } }
        );

        checkProperties( "index in inside",
          { { is0.indexInInside(), 0 },
            { is1.indexInInside(), 1 },
            { is2.indexInInside(), 2 } }
        );

        const auto& neighborA = is1.impl().outside();
        const auto& neighborB = is2.impl().outside();

        checkProperties( "neighbors centers",
          { { neighborA.geometry().center(), { 0.5, 0.5, 0.16666666666666666 } },
            { neighborB.geometry().center(), { 0.50000000000000011, 0.5, 0.83333333333333337 } } }
        );

        // check global vertex indices obtained by intersection subentities
        const auto refElement = Dune::ReferenceElements<double, 2>::general(geo.type());

        const auto vIdxLocal1 = refElement.subEntity(0, 0, 0, 2);
        const auto vIdxGlobal1 = vertexMapper.subIndex(e, vIdxLocal1, 2);
        const auto vIdxLocal2 = refElement.subEntity(0, 0, 1, 2);
        const auto vIdxGlobal2 = vertexMapper.subIndex(e, vIdxLocal2, 2);
        const auto vIdxLocal3 = refElement.subEntity(0, 0, 2, 2);
        const auto vIdxGlobal3 = vertexMapper.subIndex(e, vIdxLocal3, 2);

        checkProperties( "reference element mapping",
          { { vIdxGlobal1, 1u },
            { vIdxGlobal2, 0u },
            { vIdxGlobal3, 2u } }
        );
      }
      elementCount++;
    }

    int vi = 0;
    std::vector<int> numbOfInc = {{ 2, 2, 2, 2, 4 }};
    std::vector<int> numbOfIncV = {{ 3, 3, 3, 3, 4 }};
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
