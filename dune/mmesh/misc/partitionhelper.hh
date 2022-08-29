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
    return contains(pitype, partition(e));
  }

  bool contains(PartitionIteratorType pitype, int partitionType) const
  {
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
    return partitionType( partition(e) );
  }

  PartitionType partitionType(int partition) const
  {
    if (partition == 0)
      return InteriorEntity;

    if (partition == 1)
      return BorderEntity;

    else
      return GhostEntity;
  }

  void updatePartitions()
  {
    computePartitions();
    computeInterfacePartitions();
  }

  //! List of connected ranks
  std::vector< int > links () const
  {
    const int rank = comm().rank();
    const int size = comm().size();
    std::vector< int > links;
    if (rank > 0)
      links.push_back( rank - 1 );
    if (rank < size-1)
      links.push_back( rank + 1 );

    // TODO: is this sufficient? Improve efficiency.
    return links;
  }

  auto& comm() const
  {
    return grid().comm();
  }

private:

  //! Get partition marker
  template <class Entity>
  int partition(const Entity& e) const
  {
    static constexpr int cd = Entity::codimension;
    if (grid().comm().size() == 1)
      return 0; // interior

    if constexpr (Entity::dimension == dim)
    {
      if constexpr (Entity::codimension == 0 || Entity::codimension == dim)
        return e.impl().hostEntity()->info().partition;
      else
      {
        auto entry = partition_[cd-1].find(e.impl().id());
        if (entry != partition_[cd-1].end())
          return entry->second;
      }
    }
    else
    {
      auto entry = interfacePartition_[cd].find(e.impl().id());
      if (entry != interfacePartition_[cd].end())
        return entry->second;
    }

    DUNE_THROW(InvalidStateException, "Partition not set yet!");
    return 0;
  }

  //! Set partition marker
  template <class Entity>
  void setPartition(const Entity& e, int partition)
  {
    if constexpr (Entity::dimension == dim)
    {
      if constexpr (Entity::codimension == 0 || Entity::codimension == dim)
        e.impl().hostEntity()->info().partition = partition;
      else
        partition_[Entity::codimension-1][e.impl().id()] = partition;
    }
    else
      interfacePartition_[Entity::codimension][e.impl().id()] = partition;
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
        setPartition(e, 0); // interior
      else
        setPartition(e, 2); // ghost
    });

    // Set facets
    forEntityDim<dim-1>([this](const auto& fc){
      const auto e = grid().entity(fc);
      const auto is = grid().asIntersection( e );

      int pIn = partition( is.inside() );

      if (is.neighbor())
      {
        int pOut = partition( is.outside() );

        // interior
        if (pIn == 0 && pOut == 0)
          setPartition(e, 0);

        // border
        else if ((pIn == 0 && pOut != 0) or (pIn != 0 && pOut == 0))
          setPartition(e, 1);

        // ghost
        else
          setPartition(e, 2);
      }
      else
      {
        // interior
        if (pIn == 0)
          setPartition(e, 0);

        // ghost
        else
          setPartition(e, 2);
      }
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
          setPartition(edge, 0);

        // border
        else if (interior > 0 and ghost > 0)
          setPartition(edge, 1);

        // ghost
        else
          setPartition(edge, 2);
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
        setPartition(v, 0);

      // border
      else if (interior > 0 and ghost > 0)
        setPartition(v, 1);

      // ghost
      else
        setPartition(v, 2);
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
        // interior
        setPartition(e, 0);
      else
        // ghost
        setPartition(e, 2);
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
        setPartition(e, 0);

      // border
      else if (interior > 0 and ghost > 0)
        setPartition(e, 1);

      // ghost
      else
        setPartition(e, 2);
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
          setPartition(v, 0);

        // border
        else if (interior > 0 and ghost > 0)
          setPartition(v, 1);

        // ghost
        else
          setPartition(v, 2);
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
        if ((x > xbounds_[0] + i * dx) and (x <= xbounds_[0] + (i+1) * dx))
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
  std::array<std::unordered_map<IdType, int>, dim-1> partition_;
  std::array<std::unordered_map<IdType, int>, dim> interfacePartition_;
  const Grid& grid_;
};

} // end namespace Dune

#endif
