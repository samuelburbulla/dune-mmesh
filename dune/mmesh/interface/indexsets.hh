// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_INDEXSETS_HH
#define DUNE_MMESH_INTERFACE_INDEXSETS_HH

/** \file
 * \brief The index and id sets for the MMesh class
 */

// Dune includes
#include <dune/grid/common/indexidset.hh>
#include <dune/mmesh/grid/multiid.hh>

namespace Dune
{

  template<class GridImp>
  class MMeshInterfaceGridLeafIndexSet :
    public IndexSet<GridImp, MMeshInterfaceGridLeafIndexSet<GridImp>>
  {
    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;

    template<int codim>
    using HostGridEntity = typename GridImp::template MMeshInterfaceEntity<codim>;

  public:
    typedef std::size_t IndexType;
    typedef const std::vector< GeometryType > Types;

    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    enum {dimensionworld = std::remove_const<GridImp>::type::dimensionworld};
    enum {dimension = dimensionworld-1};

    using LocalIndexMap = std::unordered_map< std::size_t, std::size_t >;
    using IndexMap = std::unordered_map< std::array< std::size_t, dimension >, LocalIndexMap, HashUIntArray >;
    using EdgeIndexMap = std::unordered_map< std::array< std::size_t, dimension >, std::size_t, HashUIntArray >;
    using VertexIndexMap = std::unordered_map< std::size_t, std::size_t >;

    //! constructor stores reference to a grid and level
    MMeshInterfaceGridLeafIndexSet (const GridImp* grid)
      : grid_(grid)
    {}

    MMeshInterfaceGridLeafIndexSet (const MMeshInterfaceGridLeafIndexSet& leafIndexSet)
       : grid_(leafIndexSet.grid_),
         sizeOfCodim_(leafIndexSet.sizeOfCodim_),
         indexMap_(leafIndexSet.indexMap_),
         edgeIndexMap_(leafIndexSet.edgeIndexMap_),
         vertexIndices_(leafIndexSet.vertexIndices_)
    {}

    //! get index of an codim 0 entity
    template<int codim>
    std::enable_if_t< codim == 0, IndexType >
    index (const Entity< codim, dimension, GridImp, MMeshInterfaceGridEntity>& e) const
    {
      auto hostEntity = e.impl().hostEntity();

      // handle invalid entity (occurs with empty interface)
      using HostGridHandle = typename GridImp::MMeshType::template HostGridEntity<0>;
      if (!grid_->canBeMirrored(hostEntity))
      {
        if constexpr (dimension == 1)
          return (hostEntity.first->vertex(0) == grid_->getMMesh().getHostGrid().finite_faces_end()->vertex(0)) ? 1 : 0;
        else // dimension == 2
          return (hostEntity.first->vertex(0) == grid_->getMMesh().getHostGrid().finite_cells_end()->vertex(0)) ? 1 : 0;
      }

      std::array<std::size_t, dimensionworld> ids;
      for( int i = 0; i < dimensionworld; ++i )
        try {
          ids[i] = vertexIndices_.at( hostEntity.first->vertex((hostEntity.second+i+1)%(dimensionworld+1))->info().id );
        } catch (std::exception &e) {
          DUNE_THROW(InvalidStateException, e.what());
        }
      std::sort(ids.begin(), ids.end());

      std::array<std::size_t, dimension> edgeIds;
      for( int i = 0; i < dimension; ++i )
        edgeIds[i] = ids[i];

      try {
        IndexType index = indexMap_.at( edgeIds ).at( ids[dimension] );
        assert( index <= size(codim) );
        return index;
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }
    }

    //! get index of an codim 1 entity (3D)
    template<int codim>
    std::enable_if_t< codim == 1 && dimension == 2, IndexType >
    index (const Entity< codim, dimension, GridImp, MMeshInterfaceGridEntity>& e) const
    {
      auto hostEntity = e.impl().hostEntity();

      std::array<std::size_t, 2> ids;
      try {
          ids[0] = vertexIndices_.at( hostEntity.first->vertex(hostEntity.second)->info().id );
          ids[1] = vertexIndices_.at( hostEntity.first->vertex(hostEntity.third)->info().id );
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }
      std::sort(ids.begin(), ids.end());

      try {
        IndexType index = edgeIndexMap_.at( ids );
        assert( index <= size(codim) );
        return index;
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }
    }

    //! get index of an codim dimension entity
    template<int codim>
    std::enable_if_t< codim == dimension, IndexType >
    index (const Entity< codim, dimension, GridImp, MMeshInterfaceGridEntity>& e) const
    {
      auto hostEntity = e.impl().hostEntity();
      IndexType index;
      try {
        index = vertexIndices_.at( hostEntity->info().id );
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }
      assert( index <= size(codim) );
      return index;
    }

    //! get subIndex of subEntity i with given codim of an entity
    template<class Entity>
    IndexType subIndex (const Entity& e, int i, int codim)
    {
     return subIndex< Entity::codimension >( e, i, codim );
    }

    //! get subIndex of a codim dimension entity
    template<int cc>
    std::enable_if_t< cc == dimension, IndexType >
    subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      assert( i == 0 && codim == dimension );
      const HostGridEntity<dimension> hostEntity = e.impl().hostEntity();
      try {
        return vertexIndices_.at( hostEntity->info().id );
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }
    };

    //! get subIndex of a codim 0 entity
    template<int cc>
    std::enable_if_t< cc == 0, IndexType > subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      assert ( codim >= 0 && codim <= dimension );

      if ( codim == 0 )
        return index( e );
      else if ( codim == dimension )
        return index( e.template subEntity<dimension>( i ) );
      else if ( codim == 1 && dimension == 2)
        return index( e.template subEntity<1>( i ) );
      else
          DUNE_THROW( NotImplemented, "SubIndices for codim == " << codim << " and dimension == " << dimension );

      return 0;
    }

    //! provide member function subIndex for other codims but disable the usage
    template<int cc>
    std::enable_if_t< cc != 0 && cc != dimension, IndexType > subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      DUNE_THROW( NotImplemented, "SubIndex for cc != 0." );
    };

    //! get number of entities of given type
    std::size_t size (GeometryType type) const
    {
      if( type == GeometryTypes::vertex )
          return size(dimension);
      else if( type == GeometryTypes::line )
          return size(dimension-1);
      else if( type == GeometryTypes::triangle )
          return size(dimension-2);
      else
          return 0;
    }

    //! get number of entities of given codim
    std::size_t size (int codim) const
    {
      assert( (0 <= codim) && (codim <= dimension) );
      return sizeOfCodim_[codim];
    }

    /** \brief Deliver all geometry types used in this grid */
    const Types geomTypes (int codim) const
    {
      return types( codim );
    }

    /** \brief Deliver all geometry types used in this grid */
    Types types (int codim) const
    {
      switch ( dimension - codim ) {
        case 0:
          return {{ GeometryTypes::vertex }};
        case 1:
          return {{ GeometryTypes::line }};
        case 2:
          return {{ GeometryTypes::triangle }};
        default:
          DUNE_THROW(InvalidStateException, "Codim is not within 0 <= codim <= dimension.");
      }
    }

    /** \brief Return true if the given entity is contained in the index set in 2d */
    template< class EntityType >
    bool contains (const EntityType& e) const
    {
      return grid_->isInterface( e );
    }

    //! update index set in 2d
    template< int d = dimensionworld >
    std::enable_if_t< d == 2, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      vertexIndices_.clear();
      indexMap_.clear();

      // Count the finite edges and build index map
      std::size_t vertexCount = 0;
      std::size_t elementCount = 0;
      for (const auto& element : elements(grid_->leafGridView(), Partitions::all))
      {
          auto eh = &element.impl().hostEntity();
          auto vh0 = eh->first->vertex((eh->second+1)%3);
          auto vh1 = eh->first->vertex((eh->second+2)%3);

          std::size_t idx0 = vh0->info().id;
          std::size_t idx1 = vh1->info().id;

          addVertexIndex( idx0, vertexCount );
          addVertexIndex( idx1, vertexCount );

          try {
            std::size_t id0 = vertexIndices_.at( idx0 );
            std::size_t id1 = vertexIndices_.at( idx1 );

            // we store the indices for each vertex
            indexMap_[{id0}].insert( { id1, elementCount } );
            indexMap_[{id1}].insert( { id0, elementCount++ } );
          } catch (std::exception &e) {
            DUNE_THROW(InvalidStateException, e.what());
          }
      }

      // Cache sizes since it is expensive to compute them
      sizeOfCodim_[0] = elementCount;
      sizeOfCodim_[1] = vertexCount;
    }

  private:
    void addVertexIndex( std::size_t index, std::size_t& vertexCount )
    {
      auto it = vertexIndices_.find( index );
      if ( it == vertexIndices_.end() )
        vertexIndices_.insert( { index, vertexCount++ } );
    }

  public:
    //! update index set in 3d
    template< int d = dimensionworld >
    std::enable_if_t< d == 3, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      vertexIndices_.clear();
      indexMap_.clear();
      edgeIndexMap_.clear();

      // Count the finite edges and build index map
      std::size_t vertexCount = 0;
      std::size_t edgeCount = 0;
      std::size_t elementCount = 0;
      for (const auto& element : elements(grid_->leafGridView(), Partitions::all))
      {
          auto eh = &element.impl().hostEntity();
          std::array<std::size_t, dimensionworld> ids;
          for( int i = 0; i < dimensionworld; ++i )
          {
            std::size_t idx = eh->first->vertex((eh->second+i+1)%4)->info().id;
            addVertexIndex( idx, vertexCount );
            try {
              ids[i] = vertexIndices_.at( idx );
            } catch (std::exception &e) {
              DUNE_THROW(InvalidStateException, e.what());
            }
          }
          std::sort(ids.begin(), ids.end());

          // we store the indices for each codim 1 edge
          indexMap_[{ids[0], ids[1]}].insert( { ids[2], elementCount } );
          indexMap_[{ids[1], ids[2]}].insert( { ids[0], elementCount } );
          indexMap_[{ids[0], ids[2]}].insert( { ids[1], elementCount++ } );

          // store the index for edges
          for( int i = 0; i < 3; ++ i )
          {
            std::array<std::size_t, dimension> edgeIds;
            edgeIds[0] = ids[i];
            edgeIds[1] = ids[(i+1)%3];
            std::sort(edgeIds.begin(), edgeIds.end());

            if ( edgeIndexMap_.count( edgeIds ) == 0 )
              edgeIndexMap_.insert( { edgeIds, edgeCount++ } );
          }
      }

      // Cache sizes since it is expensive to compute them
      sizeOfCodim_[0] = elementCount;
      sizeOfCodim_[1] = edgeCount;
      sizeOfCodim_[2] = vertexCount;
    }

    const IndexMap& indexMap() const
    {
      return indexMap_;
    }

    const VertexIndexMap& vertexIndexMap() const
    {
      return vertexIndices_;
    }

  private:
    GridImp* grid_;
    std::array<std::size_t, dimension+1> sizeOfCodim_;
    IndexMap indexMap_;
    EdgeIndexMap edgeIndexMap_;
    VertexIndexMap vertexIndices_;
  };

  template <class GridImp>
  class MMeshInterfaceGridGlobalIdSet :
    public IdSet<GridImp, MMeshInterfaceGridGlobalIdSet<GridImp>, MMeshImpl::MultiId>
  {
    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;

    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    enum {dimensionworld = std::remove_const<GridImp>::type::dimensionworld};
    enum {dimension = dimensionworld-1};

    template<int codim>
    using HostGridEntity = typename GridImp::template MMeshInterfaceEntity<codim>;

  public:
    //! constructor stores reference to a grid
    MMeshInterfaceGridGlobalIdSet (const GridImp* g)
    {}

    //! define the type used for persistent indices
    using IdType = MMeshImpl::MultiId;

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template< int cd >
    IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return id( e );
    }

    IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e) const
    {
      return IdType( e.impl().id() );
    }

    template< int d = dimension >
    std::enable_if_t< d == 2, IdType > id (const typename std::remove_const<GridImp>::type::Traits::template Codim<1>::Entity& e) const
    {
      const auto& host = e.impl().hostEntity();
      return id( host );
    }

    //! Helper function to obtian id of MMesh codim 1 entity
    template< int cd >
    IdType id (const typename std::remove_const<GridImp>::type::MMeshType::Traits::template Codim<cd>::Entity& e) const
    {
      static_assert( cd == dimension );
      const auto& host = e.impl().hostEntity();
      return id( host );
    }

    template< int d = dimensionworld >
    std::enable_if_t< d == 3, IdType >
    id (const typename std::remove_const<GridImp>::type::template MMeshInterfaceEntity<1>& host) const
    {
      std::vector< std::size_t > ids ( 2 );
      ids[0] = host.first->vertex(host.second)->info().id;
      ids[1] = host.first->vertex(host.third)->info().id;
      std::sort( ids.begin(), ids.end() );
      return IdType( ids );
    }

    IdType id (const typename std::remove_const<GridImp>::type::Traits::template Codim<dimension>::Entity& e) const
    {
      return { e.impl().hostEntity()->info().id };
    }

    //! get id of subEntity
    /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    IdType subId (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const
    {
      assert( 0 <= codim && codim <= dimension );
      IdType dummyId ( { std::size_t(-3), std::size_t(-2) } );

      if( codim == 0 )
        return e.impl().id();

      if constexpr (dimension == 1)
      {
        // ( codim == 1 )
        {
          if (e.impl().id() != dummyId )
            return e.impl().id().vt()[i];
          else
            return IdType( std::size_t(-3 + i) );
        }
      }
      else // (dimension == 2)
      {
        if ( codim == 1 )
          return id<1>( e.impl().template subEntity<1>( i ) );

        if( codim == 2 )
          return id<2>( e.impl().template subEntity<2>( i ) );
      }

      DUNE_THROW( NotImplemented, "InterfaceGrid: subId of codim != 0 or codim != 1 or codim != 2" );
      return IdType();
    }

    //! update id set
    void update(const GridImp* grid) {}
  };

}  // end namespace Dune

#endif
