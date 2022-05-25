// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_RANGEGENERATORS_HH
#define DUNE_MMESH_RANGEGENERATORS_HH

#include <dune/common/iteratorrange.hh>

namespace Dune
{
  /**
   * \brief Elements incident to a given entity.
   */
  template<typename Entity>
  inline auto incidentElements(const Entity& entity)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    -> IteratorRange<decltype(entity.impl().incidentBegin())>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    static_assert(Entity::mydimension <= 1, "Incident element range iterator is only available for vertices and edges!");
    return IteratorRange<decltype(entity.impl().incidentBegin())>(entity.impl().incidentBegin(), entity.impl().incidentEnd());
  }

  /**
   * \brief Facets incident to a given vertex.
   */
  template<typename Vertex>
  inline auto incidentFacets(const Vertex& vertex)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    -> IteratorRange<decltype(vertex.impl().incidentFacetsBegin())>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    static_assert(Vertex::mydimension == 0, "Incident facet range iterator is only available for vertices!");
    return IteratorRange<decltype(vertex.impl().incidentFacetsBegin())>(vertex.impl().incidentFacetsBegin(),vertex.impl().incidentFacetsEnd());
  }

  /**
   * \brief Vertices incident to a given vertex.
   */
  template<typename Vertex>
  inline auto incidentVertices(const Vertex& vertex, bool includeInfinite = false)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    -> IteratorRange<decltype(vertex.impl().incidentVerticesBegin( includeInfinite ))>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    static_assert(Vertex::mydimension == 0, "Incident vertices range iterator is only available for vertices!");
    return IteratorRange<decltype(vertex.impl().incidentVerticesBegin( includeInfinite ))>(
      vertex.impl().incidentVerticesBegin( includeInfinite ),
      vertex.impl().incidentVerticesEnd( includeInfinite )
    );
  }

  /**
   * \brief All interface elements.
   */
  template<typename GridView, int codim = 1>
  inline auto interfaceElements(const GridView& gv, bool includeBoundary = false)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    -> IteratorRange<decltype(gv.grid().template interfaceBegin<codim>( includeBoundary ))>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    return IteratorRange<decltype(gv.grid().template interfaceBegin<codim>( includeBoundary ))>(
      gv.grid().template interfaceBegin<codim>( includeBoundary ),
      gv.grid().template interfaceEnd<codim>( includeBoundary )
    );
  }

  /**
   * \brief All interface vertices.
   */
  template<typename GridView>
  inline auto interfaceVertices(const GridView& gv, bool includeBoundary = false)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    -> IteratorRange<decltype(gv.grid().interfaceVerticesBegin( includeBoundary ))>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    return IteratorRange<decltype(gv.grid().interfaceVerticesBegin( includeBoundary ))>(
      gv.grid().interfaceVerticesBegin( includeBoundary ),
      gv.grid().interfaceVerticesEnd( includeBoundary )
    );
  }

  /**
   * \brief Incident interface vertices.
   */
  template<typename Vertex>
  inline auto incidentInterfaceVertices(const Vertex& vertex)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    -> IteratorRange<decltype(vertex.impl().incidentInterfaceVerticesBegin())>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    static_assert(Vertex::mydimension == 0, "Incident interface vertices range iterator is only available for interface vertices!");
    return IteratorRange<decltype(vertex.impl().incidentInterfaceVerticesBegin())>(
      vertex.impl().incidentInterfaceVerticesBegin(),
      vertex.impl().incidentInterfaceVerticesEnd()
    );
  }

  /**
   * \brief Incident interface elements.
   */
  template<typename Entity>
  inline auto incidentInterfaceElements(const Entity& entity)
#ifndef DOXYGEN_SHOULD_SKIP_THIS
    -> IteratorRange<decltype(entity.impl().incidentInterfaceElementsBegin())>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    static_assert(Entity::mydimension <= 1, "Incident interface vertices range iterator is only available for interface vertices and edges!");
    return IteratorRange<decltype(entity.impl().incidentInterfaceElementsBegin())>(
      entity.impl().incidentInterfaceElementsBegin(),
      entity.impl().incidentInterfaceElementsEnd()
    );
  }

} // namespace Dune

#endif // DUNE_MMESH_RANGEGENERATORS_HH
