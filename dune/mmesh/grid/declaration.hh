#ifndef DUNE_MMESH_GRID_DECLARATION_HH
#define DUNE_MMESH_GRID_DECLARATION_HH

namespace Dune
{
  // Forward declarations
  template<int dim, class HostGrid>
  struct MMeshFamily;

  // MMesh
  template<class HostGrid, int dim, class GridFamily = MMeshFamily<dim, HostGrid>>
  class MMesh;
}
#endif // #ifndef DUNE_MMESH_GRID_DECLARATION_HH
