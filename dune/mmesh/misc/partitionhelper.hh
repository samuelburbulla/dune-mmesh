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
    for (int i = 0; i <= dim; ++i)
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

      // border
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
    static constexpr int idim = dim-1;

    for (int i = 0; i <= idim; ++i)
      interfacePartition_[i].clear();

    // Set interior elements
    forEntityDim<idim>([this](const auto& fc){
      if (!grid().isInterface( grid().entity(fc) ))
        return;

      const auto e = grid().interfaceGrid().entity(fc);
      const auto is = grid().asIntersection(e);

      if (is.inside().impl().hostEntity()->info().rank == grid().comm().rank())
        partition(e) = 0;
      else
        partition(e) = -1;
    });

    // Set ghost elements
    forEntityDim<idim-1>([this](const auto& fc){
      if (!grid().isInterface( grid().entity(fc) ))
        return;

      const auto e = grid().interfaceGrid().entity(fc);

      bool haveIncidentInterior = [this, &e](){
        for (const auto& incident : incidentInterfaceElements(e))
          if (partition(incident) == 0)
            return true;
        return false;
      }();

      if (haveIncidentInterior)
        for (const auto& incident : incidentInterfaceElements(e))
          if (partition(incident) == -1)
            partition(incident) = 2;
    });

    // Set facets
    forEntityDim<idim-1>([this](const auto& fc){
      if (!grid().isInterface( grid().entity(fc) ))
        return;

      const auto e = grid().interfaceGrid().entity(fc);

      std::size_t count = 0, interior = 0, ghost = 0;
      for (const auto& incident : incidentInterfaceElements(e))
      {
        count++;
        if (partition(incident) == 0) interior++;
        if (partition(incident) == 2) ghost++;
      }

      // interior
      if (interior == count)
        partition(e) = 0;

      // border
      else if (interior > 0 and ghost > 0)
        partition(e) = 1;

      // ghost
      else if (ghost > 0)
        partition(e) = 2;

      // none
      else
        partition(e) = -1;
    });

    // Set vertices in 3d
    if constexpr (dim == 3)
    {
      forEntityDim<idim-2>([this](const auto& fc){
        if (!grid().isInterface( grid().entity(fc) ))
          return;

        const auto v = grid().interfaceGrid().entity(fc);

        std::size_t count = 0, interior = 0, ghost = 0;
        for (const auto& incident : incidentInterfaceElements(v))
        {
          count++;
          if (partition(incident) == 0) interior++;
          if (partition(incident) == 2) ghost++;
        }

        // interior
        if (interior == count)
          partition(v) = 0;

        // border
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
