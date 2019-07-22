// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_RANGEGENERATORS_HH
#define DUNE_MMESH_RANGEGENERATORS_HH

#include <dune/common/iteratorrange.hh>

namespace Dune
{
  /**
   * Incident elements range implementation without PartitionSet parameter. The default implementation obtains the
   * iterators by calling beginIncident() and endIncident() without specifying a partition type.
   *
   */
  template<typename Entity>
  inline auto incidentElements(const Entity& e)
    -> IteratorRange<decltype(e.impl().incidentBegin())>
  {
    static_assert(Entity::mydimension == 0, "Incident element range iterator is only available for vertices!");
    return IteratorRange<decltype(e.impl().incidentBegin())>(e.impl().incidentBegin(),e.impl().incidentEnd());
  }

} // namespace Dune

#endif // DUNE_MMESH_RANGEGENERATORS_HH
