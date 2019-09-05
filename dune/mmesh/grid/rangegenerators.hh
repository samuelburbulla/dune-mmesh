// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_RANGEGENERATORS_HH
#define DUNE_MMESH_RANGEGENERATORS_HH

#include <dune/common/iteratorrange.hh>

namespace Dune
{
  /**
   * Incident elements range implementation without PartitionSet parameter. The default implementation obtains the
   * iterators by calling incidentBegin() and incidentEnd() without specifying a partition type.
   *
   */
  template<typename Vertex>
  inline auto incidentElements(const Vertex& e)
    -> IteratorRange<decltype(e.impl().incidentBegin())>
  {
    static_assert(Vertex::mydimension == 0, "Incident element range iterator is only available for vertices!");
    return IteratorRange<decltype(e.impl().incidentBegin())>(e.impl().incidentBegin(),e.impl().incidentEnd());
  }

  /**
   * Incident facets range implementation without PartitionSet parameter. The default implementation obtains the
   * iterators by calling incidentFacetsBegin() and incidentFacetsEnd() without specifying a partition type.
   *
   */
  template<typename Vertex>
  inline auto incidentFacets(const Vertex& e)
    -> IteratorRange<decltype(e.impl().incidentFacetsBegin())>
  {
    static_assert(Vertex::mydimension == 0, "Incident facet range iterator is only available for vertices!");
    return IteratorRange<decltype(e.impl().incidentFacetsBegin())>(e.impl().incidentFacetsBegin(),e.impl().incidentFacetsEnd());
  }

  /**
   * Interface elements range implementation without PartitionSet parameter. The default implementation obtains the
   * iterators by calling interfaceBegin() and interfaceEnd() without specifying a partition type.
   *
   */
  template<typename GridView, int codim = 1>
  inline auto interfaceElements(const GridView& gv, bool includeBoundary = false)
    -> IteratorRange<decltype(gv.grid().template interfaceBegin<codim>( includeBoundary ))>
  {
    return IteratorRange<decltype(gv.grid().template interfaceBegin<codim>( includeBoundary ))>(
      gv.grid().template interfaceBegin<codim>( includeBoundary ),
      gv.grid().template interfaceEnd<codim>( includeBoundary )
    );
  }

  /**
   * Interface vertices range implementation without PartitionSet parameter. The default implementation obtains the
   * iterators by calling interfaceVerticesBegin() and interfaceVerticesEnd() without specifying a partition type.
   *
   */
  template<typename GridView>
  inline auto interfaceVertices(const GridView& gv, bool includeBoundary = false)
    -> IteratorRange<decltype(gv.grid().interfaceVerticesBegin( includeBoundary ))>
  {
    return IteratorRange<decltype(gv.grid().interfaceVerticesBegin( includeBoundary ))>(
      gv.grid().interfaceVerticesBegin( includeBoundary ),
      gv.grid().interfaceVerticesEnd( includeBoundary )
    );
  }

} // namespace Dune

#endif // DUNE_MMESH_RANGEGENERATORS_HH
