// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_MISC_PARTITIONHELPER_HH
#define DUNE_MMESH_MISC_PARTITIONHELPER_HH

#include <dune/grid/common/partitionset.hh>
#include <dune/mmesh/grid/rangegenerators.hh>

namespace Dune
{

// PartitionHelper
struct PartitionHelper
{

  template <class Entity>
  static bool contains(PartitionIteratorType pitype, const Entity& e)
  {
    if (!isInteriorOrGhost(e))
      return false;

    auto partitionType = e.partitionType();
    switch( pitype )
    {
      case All_Partition:
        return true;

      case Interior_Partition:
        return partitionType == InteriorEntity;

      case InteriorBorder_Partition:
      case Overlap_Partition:
      case OverlapFront_Partition:
        return partitionType != GhostEntity;

      case Ghost_Partition:
        return partitionType == GhostEntity;
    }
    return false;
  }

  template <class Entity>
  static bool isInteriorOrGhost(const Entity& e)
  {
    static constexpr int codim = Entity::codimension;
    static constexpr int dim = Entity::dimension;
    static constexpr int dimworld = Entity::Geometry::coorddimension;

    const auto& grid = e.impl().grid();

    if constexpr (codim == 0)
    {
      if (e.partitionType() == InteriorEntity)
        return true;

      // In the bulk, entities are ghost if a neighbor is interior
      if constexpr (dim == dimworld)
      {
        for (const auto& is : intersections(e.impl().grid().leafGridView(), e))
          if (is.neighbor())
            if (is.outside().partitionType() == InteriorEntity)
              return true;
      }
      // On the interface, entities are ghost if all adjacent bulk elements are interior or ghost
      else
      {
        for (std::size_t i = 0; i < e.subEntities(dim); ++i)
        {
          const auto v = e.template subEntity<dim>(i);

          for (const auto incident : incidentInterfaceElements(v) )
            if (incident.partitionType() == InteriorEntity)
              return true;
        }
      }
    }

    else if constexpr (codim == dim)
    {
      if constexpr (codim == dimworld)
      {
        for (const auto& incident : incidentElements( e ))
          if (isInteriorOrGhost(incident))
            return true;
      }
      else
      {
        for (const auto& incident : incidentInterfaceElements( e ))
          if (isInteriorOrGhost(incident))
            return true;
      }
    }

    else if constexpr (codim == 1)
    {
      const auto is = grid.asIntersection( e );

      if (is.neighbor())
      {
        if (isInteriorOrGhost(is.inside()) || isInteriorOrGhost(is.outside()))
          return true;
      }
      else
        return isInteriorOrGhost(is.inside());
    }

    return false;
  }

  template <class Entity>
  static int rank(const Entity& e)
  {
    int size = e.impl().grid().comm().size();
    double xmin = 0.0, xmax = 1.0;
    double dx = (xmax - xmin) / size;

    for (int i = 0; i < size; ++i)
      if (e.geometry().center()[0] <= xmin + (i+1) * dx)
        return i;
    return 0;
  }
};

} // end namespace Dune

#endif
