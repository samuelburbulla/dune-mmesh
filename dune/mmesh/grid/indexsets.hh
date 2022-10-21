// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_INDEXSETS_HH
#define DUNE_MMESH_GRID_INDEXSETS_HH

/** \file
 * \brief The index and id sets for the MMesh class
 */

// Dune includes
#include <dune/grid/common/indexidset.hh>

#include <dune/mmesh/grid/multiid.hh>

namespace Dune
{

  template<class GridImp>
  class MMeshLeafIndexSet :
    public IndexSet<GridImp, MMeshLeafIndexSet<GridImp>>
  {
    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;

    template<int codim>
    using HostGridEntity = typename GridImp::template HostGridEntity<codim>;

  public:
    typedef std::size_t IndexType;
    typedef const std::vector< GeometryType > Types;

    typedef std::unordered_map< MMeshImpl::MultiId, std::size_t > CodimIndexMap;

    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    enum {dim = std::remove_const<GridImp>::type::dimension};

    //! constructor stores reference to a grid and level
    MMeshLeafIndexSet (const GridImp* grid)
      : grid_(grid)
    {}

    MMeshLeafIndexSet (const MMeshLeafIndexSet& leafIndexSet)
       : grid_(leafIndexSet.grid_), sizeOfCodim_(leafIndexSet.sizeOfCodim_)
    {}

    //! get index of an codim 0 or codim 2 entity
    template<int codim>
    std::enable_if_t< codim == 0 || codim == dim, IndexType >
    index (const Entity< codim, dim, GridImp, MMeshEntity>& e) const
    {
      auto hostEntity = e.impl().hostEntity();
      IndexType index = hostEntity->info().index;
      assert( index <= size(codim) );
      return index;
    }

    //! get index of an codim 1 entity
    template<int codim>
    std::enable_if_t< codim == 1, IndexType >
    index (const Entity< codim, dim, GridImp, MMeshEntity>& e) const
    {
      try {
        return codimIndexMap_[0].at( grid_->globalIdSet().id( e ) );
      } catch (...) {
        DUNE_THROW(InvalidStateException, "Id of codim 1 entity not found.");
      }
    }

    //! get index of an codim 2 entity
    template<int codim>
    std::enable_if_t< codim == 2 && dim == 3, IndexType >
    index (const Entity< codim, dim, GridImp, MMeshEntity>& e) const
    {
      try {
        return codimIndexMap_[1].at( grid_->globalIdSet().id( e ) );
      } catch (...) {
        DUNE_THROW(InvalidStateException, "Id of codim 2 entity not found.");
      }
    }

    //! get subIndex of subEntity i with given codim of an entity
    template<class Entity>
    IndexType subIndex (const Entity& e, int i, int codim)
    {
     return subIndex< Entity::codimension >( e, i, codim );
    }

    //! get subIndex of a codim dim entity
    template<int cc>
    std::enable_if_t< cc == dim, IndexType >
    subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      const HostGridEntity<dim> hostEntity = e.impl().hostEntity();
      return hostEntity->info().index;
    };

    //! get subIndex of a codim 0 entity
    template<int cc>
    std::enable_if_t< cc == 0, IndexType > subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      assert ( codim >= 0 && codim <= dim );

      if ( codim == 0 )
          return index( e );
      if ( codim == dim )
          return e.impl().template subEntity<dim>( i ).impl().hostEntity()->info().index;
      else if ( codim == 1 )
          return codimIndexMap_[0].at( grid_->globalIdSet().id( e.impl().template subEntity<1>( i ) ) );
      else if ( codim == 2 )
          return codimIndexMap_[1].at( grid_->globalIdSet().id( e.impl().template subEntity<2>( i ) ) );
      else
          DUNE_THROW( InvalidStateException, "subIndex() was called for codim " << codim );

      return 0;
    }

    //! provide member function subIndex for other codims but disable the usage
    template<int cc>
    std::enable_if_t< cc != 0 && cc != dim, IndexType > subIndex (const typename std::remove_const<GridImp>::type::Traits::template Codim<cc>::Entity& e, int i, int codim) const
    {
      DUNE_THROW( NotImplemented, "SubIndices for codim != 0 || codim != dim." );
    };

    //! get number of entities of given type
    std::size_t size (GeometryType type) const
    {
      if( type == GeometryTypes::vertex )
          return size(dim);
      else if( type == GeometryTypes::line )
          return size(dim-1);
      else if( type == GeometryTypes::triangle )
          return size(dim-2);
      else if( type == GeometryTypes::tetrahedron )
          return size(dim-3);
      else
          return 0;
    }

    //! get number of entities of given codim
    std::size_t size (int codim) const
    {
      assert( (0 <= codim) && (codim <= dim) );
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
      switch ( dim - codim ) {
        case 0:
          return {{ GeometryTypes::vertex }};
        case 1:
          return {{ GeometryTypes::line }};
        case 2:
          return {{ GeometryTypes::triangle }};
        case 3:
          return {{ GeometryTypes::tetrahedron }};
        default:
          DUNE_THROW(InvalidStateException, "Codim is not within 0 <= codim <= dim.");
      }
    }

    /** \brief Return true if the given entity is contained in the index set in 2d */
    template< class EntityType, int d = dim >
    std::enable_if_t< d == 2 && EntityType::codimension == 0, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().is_face( hostEntity->vertex(0), hostEntity->vertex(1), hostEntity->vertex(2) );
    }

    template< class EntityType, int d = dim >
    std::enable_if_t< d == 2 && EntityType::codimension == 1, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().is_edge(
        hostEntity.first->vertex((hostEntity.second+1)%3),
        hostEntity.first->vertex((hostEntity.second+2)%3)
      );
    }

    template< class EntityType, int d = dim >
    std::enable_if_t< d == 2 && EntityType::codimension == 2, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().tds().is_vertex( hostEntity );
    }

    /** \brief Return true if the given entity is contained in the index set in 3d */
    template< class EntityType, int d = dim >
    std::enable_if_t< d == 3 && EntityType::codimension == 0, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().is_cell( hostEntity );
    }

    /** \brief Return true if the given entity is contained in the index set in 3d */
    template< class EntityType, int d = dim >
    std::enable_if_t< d == 3 && EntityType::codimension == 1, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().tds().is_facet(
        hostEntity.first,
        hostEntity.second
      );
    }

    /** \brief Return true if the given entity is contained in the index set in 3d */
    template< class EntityType, int d = dim >
    std::enable_if_t< d == 3 && EntityType::codimension == 2, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().tds().is_edge(
        hostEntity.first,
        hostEntity.second,
        hostEntity.third
      );
    }

    /** \brief Return true if the given entity is contained in the index set in 3d */
    template< class EntityType, int d = dim >
    std::enable_if_t< d == 3 && EntityType::codimension == 3, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().tds().is_vertex( hostEntity );
    }

    //! update index set in 2d
    template< int d = dim >
    std::enable_if_t< d == 2, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      // Store face indices within face infos
      std::size_t elementCount = 0;
      for (const auto& element : elements(grid_->leafGridView(), Partitions::all))
        element.impl().hostEntity()->info().index = elementCount++;

      // Store vertex indices within vertex infos
      std::size_t vertexCount = 0;
      for (const auto& vertex : vertices(grid_->leafGridView(), Partitions::all))
        vertex.impl().hostEntity()->info().index = vertexCount++;

      // Store the finite edge indices in a map
      codimIndexMap_[0].clear();
      std::size_t edgeCount = 0;
      for (const auto& edge : edges(grid_->leafGridView(), Partitions::all))
        codimIndexMap_[0][ grid_->globalIdSet().id( edge ) ] = edgeCount++;

      // Cache sizes since it is expensive to compute them
      sizeOfCodim_[0] = elementCount;
      sizeOfCodim_[1] = edgeCount;
      sizeOfCodim_[2] = vertexCount;
    }

    //! update index set in 3d
    template< int d = dim >
    std::enable_if_t< d == 3, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      // Store cell indices within cell infos
      std::size_t elementCount = 0;
      for (const auto& element : elements(grid_->leafGridView(), Partitions::all))
        element.impl().hostEntity()->info().index = elementCount++;

      // Store vertex indices within vertex infos
      std::size_t vertexCount = 0;
      for (const auto& vertex : vertices(grid_->leafGridView(), Partitions::all))
        vertex.impl().hostEntity()->info().index = vertexCount++;

      // Store the finite facet indices in a map
      codimIndexMap_[0].clear();
      std::size_t facetCount = 0;
      for (const auto& facet : facets(grid_->leafGridView(), Partitions::all))
        codimIndexMap_[0][ grid_->globalIdSet().id( facet ) ] = facetCount++;

      // Store the finite edge indices in a map
      codimIndexMap_[1].clear();
      std::size_t edgeCount = 0;
      for (const auto& edge : edges(grid_->leafGridView(), Partitions::all))
        codimIndexMap_[1][ grid_->globalIdSet().id( edge ) ] = edgeCount++;

      // Cache sizes since it is expensive to compute them
      sizeOfCodim_[0] = elementCount;
      sizeOfCodim_[1] = facetCount;
      sizeOfCodim_[2] = edgeCount;
      sizeOfCodim_[3] = vertexCount;
    }

    GridImp* grid_;
    std::array<std::size_t, dim+1> sizeOfCodim_;
    std::array<CodimIndexMap, dim-1> codimIndexMap_;
  };


  template <class GridImp>
  class MMeshGlobalIdSet :
    public IdSet<GridImp, MMeshGlobalIdSet<GridImp>, MMeshImpl::MultiId>
  {
    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;

    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    enum {dim = std::remove_const<GridImp>::type::dimension};

    template<int codim>
    using HostGridEntity = typename GridImp::template HostGridEntity<codim>;

  public:
    //! define the type used for persistent indices
    using IdType = MMeshImpl::MultiId;

    //! constructor stores reference to a grid
    MMeshGlobalIdSet (const GridImp* g) : grid_(g), nextVertexId_(0)
    {
      init();
    }

    //! store element and vertex id count
    void init ()
    {
      const auto& hostgrid = grid_->getHostGrid();

      // Determine nextVertexId_
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( vh->info().idWasSet && vh->info().id >= nextVertexId_ )
          nextVertexId_ = vh->info().id+1;
    }

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    std::enable_if_t< cd == dim, IdType >
    id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return IdType( e.impl().hostEntity()->info().id );
    }

    template<int cd>
    std::enable_if_t< cd != dim, IdType >
    id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return IdType( e.impl().id() );
    }

    //! get id of subEntity
    template< int d = dim >
    std::enable_if_t< d == 2, IdType >
    subId (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const
    {
      assert( 0 <= codim && codim <= dim );
      IdType dummyId ( { std::size_t(-4), std::size_t(-3), std::size_t(-2) } );
      switch( codim )
      {
        case 0:
          return e.impl().id();
        case 1:
          if (e.impl().id() != dummyId )
          {
            auto id0 = e.impl().id().vt()[i < 2 ? 0 : 1];
            auto id1 = e.impl().id().vt()[i == 0 ? 1 : 2];
            return IdType( { std::min(id0, id1), std::max(id0, id1) } );
          }
          else
          {
            std::size_t id0 = (i < 2 ? -4 : -3);
            std::size_t id1 = (i == 0 ? -3 : -2);
            return IdType( { id0, id1 } );
          }
        case 2:
          if (e.impl().id() != dummyId )
            return e.impl().id().vt()[i];
          else
            return IdType( std::size_t(-4 + i) );
      };
      return IdType();
    }

    template< int d = dim >
    std::enable_if_t< d == 3, IdType >
    subId (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const
    {
      assert( 0 <= codim && codim <= dim );
      switch( codim )
      {
        case 0: return id<0>( e.impl().template subEntity<0>( i ) );
        case 1: return id<1>( e.impl().template subEntity<1>( i ) );
        case 2: return id<2>( e.impl().template subEntity<2>( i ) );
        case 3: return id<3>( e.impl().template subEntity<3>( i ) );
      };
      return IdType();
    }

    //! update id set in 2d
    template< int d = dim >
    std::enable_if_t< d == 2, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      // Store vertex ids within vertex infos
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( !vh->info().idWasSet )
        {
          vh->info().id = nextVertexId_++;
          vh->info().idWasSet = true;
        }

      // Compute mapping DUNE vertex index to CGAL vertex index
      for ( auto fh = hostgrid.all_faces_begin(); fh != hostgrid.all_faces_end(); ++fh)
        fh->info().cgalIndex = MMeshImpl::computeCGALIndices<decltype(fh), 2>( fh );
    }

    //! update id set in 3d
    template< int d = dim >
    std::enable_if_t< d == 3, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      // Store vertex ids within vertex infos
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( !vh->info().idWasSet )
        {
          vh->info().id = nextVertexId_++;
          vh->info().idWasSet = true;
        }

      // Compute mapping DUNE vertex index to CGAL vertex index
      for ( auto ch = hostgrid.all_cells_begin(); ch != hostgrid.all_cells_end(); ++ch)
        ch->info().cgalIndex = MMeshImpl::computeCGALIndices<decltype(ch), 3>( ch );
    }

    //! advanced method to set the id of a vertex manually
    std::size_t setNextId( HostGridEntity<dim> vh ) const
    {
      assert( !vh->info().idWasSet );
      vh->info().id = nextVertexId_++;
      vh->info().idWasSet = true;
      return vh->info().id;
    }

    GridImp* grid_;
    mutable std::size_t nextVertexId_;
  };

}  // end namespace Dune

#endif
