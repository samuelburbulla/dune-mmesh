#ifndef DUNE_MMESH_GRID_DECLARATION_HH
#define DUNE_MMESH_GRID_DECLARATION_HH

namespace Dune
{

  namespace MMeshDefaults
  {
    template< int dim >
    class Triangulation;

    template< int dim >
    class Delaunay;
  }

  // Forward declarations
  template<int dim, class HostGrid>
  struct MMeshFamily;

  // MMesh
  template<class HostGrid, int dim>
  class MMesh;

  // Type of wrapper triangulation
  template< int dim >
  class TriangulationWrapper;

  // Type shortcut with default triangulation
  template< int dim >
  using MovingMesh = MMesh< TriangulationWrapper<dim>, dim >;

  // Type of Delaunay wrapper triangulation
  template< int dim >
  class DelaunayTriangulationWrapper;

  // Type shortcut with delaunay triangulation
  template<int dim>
  using DelaunayTriangulation = MMesh< DelaunayTriangulationWrapper<dim>, dim >;
}
#endif // #ifndef DUNE_MMESH_GRID_DECLARATION_HH
