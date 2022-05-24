// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_MISC_PARTITIONHELPER_HH
#define DUNE_MMESH_MISC_PARTITIONHELPER_HH

#include <dune/grid/common/partitionset.hh>
#include <dune/mmesh/grid/rangegenerators.hh>

namespace Dune
{

// PartitionHelper
template <class Grid>
struct PartitionHelper
{
  static constexpr int dim = Grid::dimension;
  using IdType = typename Grid::IdType;

  PartitionHelper(const Grid& grid) : grid_(grid) {}

  void distribute()
  {
    if (grid().comm().size() == 1)
      return;

    setRanks();
    computePartitions();
    computeInterfacePartitions();
  }

  template <class Entity>
  bool contains(PartitionIteratorType pitype, const Entity& e) const
  {
    int partitionType = partition(e);

    if (partitionType == -1)
      return false;

    switch( pitype )
    {
      case All_Partition:
        return true;

      case Interior_Partition:
        return partitionType == 0;

      case InteriorBorder_Partition:
      case Overlap_Partition:
      case OverlapFront_Partition:
        return partitionType != 2;

      case Ghost_Partition:
        return partitionType == 2;
    }
    return false;
  }

  template <class Entity>
  PartitionType partitionType(const Entity& e) const
  {
    int partitionType = partition(e);

    if (partitionType == 0)
      return InteriorEntity;

    if (partitionType == 1)
      return BorderEntity;

    if (partitionType == 2)
      return GhostEntity;

    // we return GhostEntity even if it's not actually ghost (partitionType == -1)
    return GhostEntity;
  }

private:

  //! Get partition marker
  template <class Entity>
  int partition(const Entity& e) const
  {
    if (grid().comm().size() == 1)
      return 0; // interior

    if constexpr (Entity::dimension == dim)
    {
      auto entry = partition_[Entity::codimension].find(e.impl().id());
      if (entry != partition_[Entity::codimension].end())
        return entry->second;
    }
    else
    {
      auto entry = interfacePartition_[Entity::codimension].find(e.impl().id());
      if (entry != interfacePartition_[Entity::codimension].end())
        return entry->second;
    }

    return 0;
  }

  //! Set partition marker
  template <class Entity>
  int& partition(const Entity& e)
  {
    if constexpr (Entity::dimension == dim)
      return partition_[Entity::codimension][e.impl().id()];
    else
      return interfacePartition_[Entity::codimension][e.impl().id()];
  }

  //! Compute partition type for every entity
  void computePartitions()
  {
    for (int i = 0; i < dim+1; ++i)
      partition_[i].clear();

    // Set elements
    forEntityDim<dim>([this](const auto& fc){
      const auto e = grid().entity(fc);
      if (e.impl().hostEntity()->info().rank == grid().comm().rank())
        partition(e) = 0;
      else
        partition(e) = -1;
    });

    // Set interior and boundary facets and find ghost elements
    forEntityDim<dim-1>([this](const auto& fc){
      const auto e = grid().entity(fc);
      const auto is = grid().asIntersection( e );

      int pIn = partition( is.inside() );

      if (is.neighbor())
      {
        int pOut = partition( is.outside() );

        // interior
        if (pIn == 0 && pOut == 0)
          partition(e) = 0;

        // border
        else if (pIn == 0 && pOut != 0)
        {
          partition(e) = 1;
          partition(is.outside()) = 2; // set outside as ghost
        }

        else if (pIn != 0 && pOut == 0)
        {
          partition(e) = 1;
          partition(is.inside()) = 2; // set inside as ghost
        }

        // none
        else
          partition(e) = -1;
      }
      else
      {
        // interior
        if (pIn == 0)
          partition(e) = 0;
      }
    });

    // Set for ghost facets
    forEntityDim<dim-1>([this](const auto& fc){
      const auto e = grid().entity(fc);
      const auto is = grid().asIntersection( e );

      if (partition(e) == 1)
        return;

      if (partition( is.inside() ) == 2)
        partition(e) = 2;

      if (is.neighbor())
        if (partition( is.outside() ) == 2)
          partition(e) = 2;
    });

    // Codim 2 in 3D
    if constexpr (dim == 3)
    {
      forEntityDim<1>([this](const auto& fc){
        const auto edge = grid().entity(fc);
        std::size_t count = 0, interior = 0, ghost = 0;
        for (const auto e : incidentElements(edge))
        {
          count++;
          if (partition(e) == 0) interior++;
          if (partition(e) == 2) ghost++;
        }

        // interior
        if (interior == count)
          partition(edge) = 0;

        // ghost
        else if (interior > 0 and ghost > 0)
          partition(edge) = 1;

        // ghost
        else if (ghost > 0)
          partition(edge) = 2;

        // none
        else
          partition(edge) = -1;
      });
    }

    // Set vertices
    forEntityDim<0>([this](const auto& fc){
      const auto v = grid().entity(fc);
      std::size_t count = 0, interior = 0, ghost = 0;
      for (const auto e : incidentElements(v))
      {
        count++;
        if (partition(e) == 0) interior++;
        if (partition(e) == 2) ghost++;
      }

      // interior
      if (interior == count)
        partition(v) = 0;

      // ghost
      else if (interior > 0 and ghost > 0)
        partition(v) = 1;

      // ghost
      else if (ghost > 0)
        partition(v) = 2;

      // none
      else
        partition(v) = -1;
    });
  }

  //! Compute partition type for every interface entity
  void computeInterfacePartitions()
  {
//    for (int i = 0; i < dim; ++i)
//      interfacePartition_[i].clear();
//
//    // Set elements
//    forEntityDim<dim-1>([this](const auto& fc){
//      const auto e = grid().interfaceGrid().entity(fc);
//      partition(e) = 0;
//    });
  }

  //! Set rank for every entity. We use a naiv partition slicing the x-axis here.
  void setRanks()
  {
    xbounds_[0] = std::numeric_limits<double>::max();
    xbounds_[1] = std::numeric_limits<double>::lowest();

    forEntityDim<0>([this](const auto& fc){
      auto v = grid().entity(fc);
      auto x = v.geometry().center()[0];
      xbounds_[0] = std::min(xbounds_[0], x);
      xbounds_[1] = std::max(xbounds_[1], x);
    });

    auto rank = [this](const auto& e)
    {
      int size = grid().comm().size();
      double dx = (xbounds_[1] - xbounds_[0]) / size;

      auto x = e.geometry().center()[0];
      for (int i = 0; i < size; ++i)
        if (x <= xbounds_[0] + (i+1) * dx)
          return i;
      return 0;
    };

    forEntityDim<dim>([this, rank](const auto& fc){
      auto e = grid().entity(fc);
      e.impl().hostEntity()->info().rank = rank( e );
    });
  }

  template <int edim, class F>
  void forEntityDim(const F& f)
  {
    if constexpr (edim == dim)
    {
      if constexpr (dim == 2)
        for (auto fc = grid().getHostGrid().finite_faces_begin(); fc != grid().getHostGrid().finite_faces_end(); ++fc)
          f( fc );

      if constexpr (dim == 3)
        for (auto fc = grid().getHostGrid().finite_cells_begin(); fc != grid().getHostGrid().finite_cells_end(); ++fc)
          f( fc );
    }

    if constexpr (edim == 0)
      for (auto fc = grid().getHostGrid().finite_vertices_begin(); fc != grid().getHostGrid().finite_vertices_end(); ++fc)
        f( fc );

    if constexpr (edim == 1)
      for (auto fc = grid().getHostGrid().finite_edges_begin(); fc != grid().getHostGrid().finite_edges_end(); ++fc)
        f( *fc );

    if constexpr (dim == 3 && edim == 2)
      for (auto fc = grid().getHostGrid().finite_facets_begin(); fc != grid().getHostGrid().finite_facets_end(); ++fc)
        f( *fc );
  }

  const Grid& grid() const { return grid_; }

  std::array<double, 2> xbounds_;
  std::array<std::unordered_map<IdType, int>, dim+1> partition_;
  std::array<std::unordered_map<IdType, int>, dim> interfacePartition_;
  const Grid& grid_;
};

} // end namespace Dune

#endif




//
//static constexpr int codim = Entity::codimension;
//static constexpr int edim = Entity::dimension;
//static constexpr int dimworld = Entity::Geometry::coorddimension;
//
//if constexpr (codim == 0)
//{
//  if (e.impl().hostEntity()->info().rank == grid().comm().rank())
//    return InteriorEntity;
//  else
//    return GhostEntity;
//}
//
//if constexpr (edim == dimworld)
//{
//  if constexpr (codim == dim)
//  {
//    if (grid().comm().size() == 1)
//      return InteriorEntity;
//
//    std::size_t interior = 0, count = 0;
//    for (const auto& e : incidentElements( e ))
//    {
//      count++;
//      if (e.partitionType() == InteriorEntity)
//        interior++;
//    }
//
//    if (interior == count)
//      return InteriorEntity;
//    else if (interior == 0)
//      return GhostEntity;
//    else
//      return BorderEntity;
//  }
//  else if constexpr (codim == 1)
//  {
//    const auto is = grid().asIntersection( e );
//
//    auto pIn = is.inside().partitionType();
//    if (is.neighbor())
//    {
//      auto pOut = is.inside().partitionType();
//      if (pIn == InteriorEntity && pOut == InteriorEntity)
//        return InteriorEntity;
//
//      if ((pIn == InteriorEntity && pOut == GhostEntity)
//          || (pIn == GhostEntity && pOut == InteriorEntity))
//        return BorderEntity;
//
//      return GhostEntity;
//    }
//    else
//      return (pIn == InteriorEntity) ? InteriorEntity : GhostEntity;
//  }
//  else if constexpr (codim == 2)
//  {
//    return InteriorEntity; // TODO
//  }
//  else
//    return InteriorEntity;
//}
//else // (edim != dimworld)
//{
//  if constexpr (codim == 0)
//  {
//    const auto is = grid().getMMesh().asIntersection( *this );
//
//    auto pIn = is.inside().partitionType();
//    if (is.neighbor())
//    {
//      auto pOut = is.outside().partitionType();
//      if (pIn == InteriorEntity && pOut == InteriorEntity)
//        return InteriorEntity;
//
//      if (pIn == GhostEntity && pOut == GhostEntity)
//        return GhostEntity;
//
//      if (is.inside().impl().hostEntity()->info().rank == grid().comm().rank())
//        return InteriorEntity;
//      else
//        return GhostEntity;
//    }
//    else
//      return (pIn == InteriorEntity) ? InteriorEntity : GhostEntity;
//  }
//  else
//  {
//    if constexpr (codim == dim)
//    {
//      int interior = 0, count = 0;
//      for (const auto& incident : incidentInterfaceElements( e ))
//      {
//        count++;
//        if (incident.partitionType() == InteriorEntity)
//          interior++;
//      }
//
//      if (interior == count)
//        return InteriorEntity;
//      else if (interior == 0)
//        return GhostEntity;
//      else
//        return BorderEntity;
//    }
//    if constexpr (codim == 2)
//    {
//      return InteriorEntity; // TODO
//    }
//    else
//      return InteriorEntity;
//  }
//}
//}
//
//
//
//
//// Codim 0 (interface)
//template<class Entity, int dim>
//int pt<Entity, 0, dim, true> (const Entity& e) const
//{
//  const auto is = grid().getMMesh().asIntersection( e );
//
//  int pIn = mmeshPartitionType( is.inside() );
//  if (is.neighbor())
//  {
//    int pOut = mmeshPartitionType( is.outside() );
//
//    if (pIn == 0 && pOut == 0)
//      return 0;
//
//    if (pIn == 2 && pOut == 2)
//      return GhostEntity;
//
//    if (is.inside().impl().hostEntity()->info().rank == grid().comm().rank())
//      return 0;
//    else
//      return 2;
//  }
//  else
//    return (pIn == InteriorEntity) ? InteriorEntity : GhostEntity;
//}
//
//
//if constexpr (codim == 0)
//{
//  if (e.impl().hostEntity()->info().rank == grid().comm().rank())
//    return true;
//
//  // In the bulk, entities are ghost if a neighbor is interior
//  if constexpr (dim == dimworld)
//  {
//    for (const auto& is : intersections(grid().leafGridView(), e))
//      if (is.neighbor())
//        if (is.outside().hostEntity()->info().rank == grid().comm().rank())
//          return true;
//  }
//  // On the interface, entities are ghost if some adjacent bulk element is interior
//  else
//  {
//    for (std::size_t i = 0; i < e.subEntities(dim); ++i)
//    {
//      const auto f = e.template subEntity<1>(i);
//
//      for (const auto incident : incidentInterfaceElements(f) )
//        if (incident.partitionType() == InteriorEntity)
//          return true;
//    }
//  }
//}
//
//else if constexpr (codim == dim)
//{
//  if constexpr (codim == dimworld)
//  {
//    for (const auto& incident : incidentElements( e ))
//      if (isPartOfThisRank(incident))
//        return true;
//  }
//  else
//  {
//    for (const auto& incident : incidentInterfaceElements( e ))
//      if (isPartOfThisRank(incident))
//        return true;
//  }
//}
//
//else if constexpr (codim == 1)
//{
//  if constexpr (codim == dimworld)
//  {
//    const auto is = grid().asIntersection( e );
//
//    if (is.neighbor())
//    {
//      if (isPartOfThisRank(is.inside()) || isPartOfThisRank(is.outside()))
//        return true;
//    }
//    else
//      return isPartOfThisRank(is.inside());
//  }
//  else
//    return true; // TODO
//}
//
//else if constexpr (codim == 2)
//{
//  return true; // TODO
//}
//
//return false;
//}
