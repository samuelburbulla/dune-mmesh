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
  using ConnectivityType = std::unordered_set<int>;
  using LinksType = std::vector<int>;
  using LeafIterator = typename Grid::LeafIterator::Implementation::HostGridLeafIterator;

  PartitionHelper(const Grid& grid) : grid_(grid) {}

  void distribute()
  {
    // Initialize leafBegin_ and leafEnd_
    if constexpr (dim == 2)
    {
      leafBegin_ = grid().getHostGrid().finite_faces_begin();
      leafEnd_ = grid().getHostGrid().finite_faces_end();
    }
    else if constexpr (dim == 3)
    {
      leafBegin_ = grid().getHostGrid().finite_cells_begin();
      leafEnd_ = grid().getHostGrid().finite_cells_end();
    }

    if (grid().comm().size() == 1)
      return;

    setRanks();
    computePartitions();
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
  }

  //! List of connected ranks
  const LinksType& links () const
  {
    return links_;
  }

  auto& comm() const
  {
    return grid().comm();
  }

  //! Get connectivity (list of ranks)
  template <class Entity>
  ConnectivityType connectivity(const Entity &e) const
  {
    if constexpr (Entity::dimension == dim)
      return e.impl().hostEntity()->info().connectivity;
    else
      try {
        return interfaceConnectivity_.at(e.impl().id());
      } catch (std::out_of_range&) {
        return ConnectivityType();
      }
  }

  //! Get rank of an entity
  template <class Entity>
  int rank(const Entity &e) const
  {
    if constexpr (Entity::dimension == dim)
      return e.impl().hostEntity()->info().rank;
    else
      return grid().asIntersection(e).inside().impl().hostEntity()->info().rank;
  }

  const LeafIterator& leafInteriorBegin() const
  {
    return leafBegin_;
  }

  const LeafIterator& leafInteriorEnd() const
  {
    return leafEnd_;
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

      // handle empty interface
      return -1;
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

  //! Add connectivity
  template <class Entity>
  void addConnectivity(const Entity &e, int rank)
  {
    if constexpr (Entity::dimension == dim)
      e.impl().hostEntity()->info().connectivity.insert(rank);
    else
      interfaceConnectivity_[e.impl().id()].insert(rank);

    if (partition(e) == 0)
      if (std::find(links_.begin(), links_.end(), rank) == links_.end())
        links_.push_back(rank);
  }

  //! Clear connectivity
  template <class Entity>
  void clearConnectivity(const Entity &e)
  {
    if constexpr (Entity::dimension == dim)
      e.impl().hostEntity()->info().connectivity.clear();
    else
      interfaceConnectivity_[e.impl().id()].clear();
  }

  //! Set rank for every entity. We use a naiv partitioning using the entity iterator.
  void setRanks()
  {
    int rank = grid().comm().rank();
    int size = grid().comm().size();

    std::size_t N;
    if constexpr (dim == 2)
      N = grid().getHostGrid().number_of_faces();
    else
      N = grid().getHostGrid().number_of_cells();

    std::size_t i = 0;
    forEntityDim<dim>([this, &i, &rank, &size, &N](const auto& fc)
    {
      auto getRank = [&size, &N](std::size_t i){ return i * size / N; };
      int r = getRank(i);

      // store begin iterator
      if (r == rank && getRank(i-1) == rank-1)
        this->leafBegin_ = fc;

      // store end iterator
      if (r == rank+1 && getRank(i-1) == rank)
        this->leafEnd_ = fc;

      auto e = grid().entity(fc);
      e.impl().hostEntity()->info().rank = r;
      i++;
    });
  }

  //! Compute partition type for every entity
  void computePartitions()
  {
    links_.clear();

    for (int i = 0; i <= dim; ++i)
      partition_[i].clear();

    // Set interior elements
    forEntityDim<dim>([this](const auto& fc){
      const auto e = grid().entity(fc);
      if (rank(e) == grid().comm().rank())
        setPartition(e, 0); // interior
      else
        setPartition(e, -1); // none

      clearConnectivity(e);
    });

    // Set ghosts
    forEntityDim<dim-1>([this](const auto& fc){
      const auto e = grid().entity(fc);
      const auto is = grid().asIntersection( e );
      if (is.neighbor())
      {
        const auto& inside = is.inside();
        const auto& outside = is.outside();
        int pIn = partition( inside );
        int pOut = partition( outside );

        // border
        if ((pIn == 0 && pOut != 0) or (pIn != 0 && pOut == 0))
        {
          setPartition(pIn == 0 ? outside : inside, 2); // ghost
          addConnectivity(inside, rank(outside));
          addConnectivity(outside, rank(inside));
        }
      }
    });

    // Compute interface partitions based on Codim 0 entities, this might add further ghost elements to the bulk
    computeInterfacePartitions();

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
        else if ((pIn == 0 && pOut == 2) or (pIn == 2 && pOut == 0))
          setPartition(e, 1);

        // ghost
        else if (pIn == 2 || pOut == 2)
          setPartition(e, 2);

        // none
        else
          setPartition(e, -1);
      }
      else
      {
        // interior
        if (pIn == 0)
          setPartition(e, 0);

        // ghost
        else if (pIn == 2)
          setPartition(e, 2);

        // none
        else
          setPartition(e, -1);
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
        else if (interior == 0 and ghost > 0)
          setPartition(edge, 2);

        // none
        else
          setPartition(edge, -1);
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
      else if (interior == 0 and ghost > 0)
        setPartition(v, 2);

      // none
      else
        setPartition(v, -1);
    });
  }

public:
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

      clearConnectivity(e);

      if (rank(e) == grid().comm().rank())
      {
        setPartition(e, 0); // interior

        // add connectivity to this entity being ghost on outside rank
        if (is.neighbor())
          if (rank(is.outside()) != grid().comm().rank())
            addConnectivity(e, rank(is.outside()));
      }
      else
        setPartition(e, -1); // none

      // if outside bulk entity is interior, the interface element is at least ghost
      if (is.neighbor())
        if (rank(is.outside()) == grid().comm().rank())
          if (partition(e) == -1)
          {
            setPartition(e, 2); // ghost
            addConnectivity(e, rank(is.inside()));
          }
    });

    // Set ghosts
    forEntityDim<idim-1>([this](const auto& fc){
      if (!grid().isInterface( grid().entity(fc) ))
        return;

      // convert to interface vertex
      const auto e = grid().interfaceGrid().entity(fc);

      std::size_t count = 0, interior = 0, other = 0;
      for (const auto& incident : incidentInterfaceElements(e))
      {
        count++;
        if (partition(incident) == 0) interior++;
        if (partition(incident) != 0) other++;
      }

      // if we find both interior and other entities we set none to ghost
      if (interior > 0 and other > 0)
        for (const auto& incident : incidentInterfaceElements(e))
          if (partition(incident) == -1)
          {
            setPartition(incident, 2); // ghost

            // make sure that both adjacent bulk entities are at least ghost
            auto intersection = grid().asIntersection(incident);
            if (partition(intersection.inside()) == -1)
              setPartition(intersection.inside(), 2); // ghost
            if (intersection.neighbor())
              if (partition(intersection.outside()) == -1)
                setPartition(intersection.outside(), 2); // ghost
          }
    });

    // Set facets
    forEntityDim<idim-1>([this](const auto& fc){
      if (!grid().isInterface( grid().entity(fc) ))
        return;

      // convert to interface vertex
      const auto e = grid().interfaceGrid().entity(fc);

      std::size_t count = 0, interior = 0, ghost = 0;
      std::unordered_set<int> connectivity;
      for (const auto& incident : incidentInterfaceElements(e))
      {
        count++;
        if (partition(incident) == 0) interior++;
        if (partition(incident) == 2) ghost++;
        connectivity.insert( rank(incident) );
      }

      // interior
      if (interior == count)
        setPartition(e, 0);

      // border
      else if (interior > 0 and ghost > 0)
      {
        setPartition(e, 1);

        for (const auto& incident : incidentInterfaceElements(e))
          for (auto r : connectivity)
            if (r != rank(incident))
              addConnectivity(incident, r);
      }

      // ghost
      else if (ghost > 0)
        setPartition(e, 2);

      // none
      else
        setPartition(e, -1);
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
        else if (ghost > 0)
          setPartition(v, 2);

        // none
        else
          setPartition(v, -1);
      });
    }
  }

private:
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
  std::unordered_map<IdType, ConnectivityType> interfaceConnectivity_;
  LeafIterator leafBegin_, leafEnd_;
  LinksType links_;
  const Grid& grid_;
};

} // end namespace Dune

#endif
