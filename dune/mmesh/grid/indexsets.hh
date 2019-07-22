// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_INDEXSETS_HH
#define DUNE_MMESH_GRID_INDEXSETS_HH

/** \file
 * \brief The index and id sets for the MMesh class
 */

// Dune includes
#include <dune/grid/common/indexidset.hh>

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
      DUNE_THROW( NotImplemented, "Index for codim 1 entities." );
      return 0;
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
      assert ( codim > 0 && codim <= dim );
      const HostGridEntity<0> hostEntity = e.impl().hostEntity();

      if ( codim == dim )
          return hostEntity->vertex( i )->info().index;
      else
          DUNE_THROW( NotImplemented, "SubIndices for codim != 0 || codim != dim." );

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
    const Types& geomTypes (int codim) const
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
    std::enable_if_t< d == 2, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().is_face( hostEntity->vertex(0), hostEntity->vertex(1), hostEntity->vertex(2) );
    }

    /** \brief Return true if the given entity is contained in the index set in 3d */
    template< class EntityType, int d = dim >
    std::enable_if_t< d == 3, bool >
    contains (const EntityType& e) const
    {
      const auto hostEntity = e.impl().hostEntity();
      return grid_->getHostGrid().is_cell( hostEntity );
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
      for ( auto fc = hostgrid.finite_faces_begin(); fc != hostgrid.finite_faces_end(); ++fc)
        fc->info().index = elementCount++;

      // Store vertex indices within vertex infos
      std::size_t vertexCount = 0;
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        vh->info().index = vertexCount++;

      // Count the finite edges
      std::size_t edgeCount = 0;
      for ( auto eh = hostgrid.finite_edges_begin(); eh != hostgrid.finite_edges_end(); ++eh)
        edgeCount++;

      // Cache sizes since it is expensive to compute them
      sizeOfCodim_[0] = hostgrid.number_of_faces();
      sizeOfCodim_[1] = edgeCount;
      sizeOfCodim_[2] = hostgrid.number_of_vertices();
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
      for ( auto fc = hostgrid.finite_cells_begin(); fc != hostgrid.finite_cells_end(); ++fc)
          fc->info().index = elementCount++;

      // Store vertex indices within vertex infos
      std::size_t vertexCount = 0;
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
          vh->info().index = vertexCount++;

      // Cache sizes since it is expensive to compute them
      sizeOfCodim_[0] = hostgrid.number_of_finite_cells();
      sizeOfCodim_[1] = hostgrid.number_of_finite_facets();
      sizeOfCodim_[2] = hostgrid.number_of_finite_edges();
      sizeOfCodim_[3] = hostgrid.number_of_vertices();
    }

    GridImp* grid_;
    std::array<std::size_t, dim+1> sizeOfCodim_;
  };

  template <class GridImp>
  class MMeshGlobalIdSet :
    public IdSet<GridImp,MMeshGlobalIdSet<GridImp>,
                 /*IdType=*/std::size_t>
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
    //! constructor stores reference to a grid
    MMeshGlobalIdSet (const GridImp* g) : grid_(g), nextElementId_(0), nextVertexId_(0)
    {
      init();
    }

    //! store element and vertex id count
    template< int d = dim >
    std::enable_if_t< d == 2, void >
    init ()
    {
      const auto& hostgrid = grid_->getHostGrid();

      // Determine nextElementId_
      for ( auto fc = hostgrid.finite_faces_begin(); fc != hostgrid.finite_faces_end(); ++fc)
        if( fc->info().idWasSet && fc->info().id >= nextElementId_ )
          nextElementId_ = fc->info().id+1;

      // Determine nextVertexId_
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( vh->info().idWasSet && vh->info().id >= nextVertexId_ )
          nextVertexId_ = vh->info().id+1;
    }

    //! store element and vertex id count
    template< int d = dim >
    std::enable_if_t< d == 3, void >
    init ()
    {
      const auto& hostgrid = grid_->getHostGrid();

      // Determine nextElementId_
      for ( auto fc = hostgrid.finite_cells_begin(); fc != hostgrid.finite_cells_end(); ++fc)
        if( fc->info().idWasSet && fc->info().id >= nextElementId_ )
          nextElementId_ = fc->info().id+1;

      // Determine nextVertexId_
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( vh->info().idWasSet && vh->info().id >= nextVertexId_ )
          nextVertexId_ = vh->info().id+1;
    }

    //! define the type used for persistent indices
    typedef std::size_t IdType;

    //! get id of an entity
    /*
       We use the remove_const to extract the Type from the mutable class,
       because the const class is not instantiated yet.
     */
    template<int cd>
    std::enable_if_t< cd == 0 || cd == dim, IdType >
    id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return e.impl().hostEntity()->info().id;
    }

    template<int cd>
    std::enable_if_t< cd == 1, IdType >
    id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      DUNE_THROW( NotImplemented, "Id for codim 1 entities." );
      return 0;
    }

    //! get id of subEntity
    /*
        We use the remove_const to extract the Type from the mutable class,
        because the const class is not instantiated yet.
     */
    IdType subId (const typename std::remove_const<GridImp>::type::Traits::template Codim<0>::Entity& e, int i, int codim) const
    {
      assert( codim == 0 || codim == dim );
      if( codim == dim )
        return e.impl().hostEntity()->vertex( i )->info().id;
      else
        return e.impl().hostEntity()->info().id;
    }

    //! update id set in 2d
    template< int d = dim >
    std::enable_if_t< d == 2, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      // Store cell ids within cell infos
      for ( auto fc = hostgrid.finite_faces_begin(); fc != hostgrid.finite_faces_end(); ++fc)
        if( !fc->info().idWasSet )
        {
          fc->info().id = nextElementId_++;
          fc->info().idWasSet = true;
        }

      // Store vertex ids within vertex infos
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( !vh->info().idWasSet )
        {
          vh->info().id = nextVertexId_++;
          vh->info().idWasSet = true;
        }
    }

    //! update id set in 3d
    template< int d = dim >
    std::enable_if_t< d == 3, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      // Store cell ids within cell infos
      for ( auto fc = hostgrid.finite_cells_begin(); fc != hostgrid.finite_cells_end(); ++fc)
        if( !fc->info().idWasSet )
        {
          fc->info().id = nextElementId_++;
          fc->info().idWasSet = true;
        }

      // Store vertex ids within vertex infos
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( !vh->info().idWasSet )
        {
          vh->info().id = nextVertexId_++;
          vh->info().idWasSet = true;
        }
    }

    GridImp* grid_;
    std::size_t nextElementId_, nextVertexId_;
  };


  template<class GridImp>
  class MMeshLocalIdSet :
    public IdSet<GridImp,MMeshLocalIdSet<GridImp>,
                 /*IdType=*/std::size_t>
  {
  private:
    typedef typename std::remove_const<GridImp>::type::HostGridType HostGrid;

    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    enum {dim = std::remove_const<GridImp>::type::dimension};

    template<int codim>
    using HostGridEntity = typename GridImp::template HostGridEntity<codim>;

  public:
    //! define the type used for persistent local ids
    typedef std::size_t IdType;

    //! constructor stores reference to a grid
    MMeshLocalIdSet (const GridImp* g) : grid_(g), nextElementId_(0), nextVertexId_(0)
    {
      init();
    }

    //! store element and vertex id count
    template< int d = dim >
    std::enable_if_t< d == 2, void >
    init ()
    {
      const auto& hostgrid = grid_->getHostGrid();

      // Determine nextElementId_
      for ( auto fc = hostgrid.finite_faces_begin(); fc != hostgrid.finite_faces_end(); ++fc)
        if( fc->info().idWasSet && fc->info().id >= nextElementId_ )
          nextElementId_ = fc->info().id+1;

      // Determine nextVertexId_
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( vh->info().idWasSet && vh->info().id >= nextVertexId_ )
          nextVertexId_ = vh->info().id+1;
    }

    //! store element and vertex id count
    template< int d = dim >
    std::enable_if_t< d == 3, void >
    init ()
    {
      const auto& hostgrid = grid_->getHostGrid();

      // Determine nextElementId_
      for ( auto fc = hostgrid.finite_cells_begin(); fc != hostgrid.finite_cells_end(); ++fc)
        if( fc->info().idWasSet && fc->info().id >= nextElementId_ )
          nextElementId_ = fc->info().id+1;

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
    std::enable_if_t< cd == 0 || cd == dim, IdType >
    id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      return e.impl().hostEntity()->info().id;
    }

    template<int cd>
    std::enable_if_t< cd == 1, IdType >
    id (const typename std::remove_const<GridImp>::type::Traits::template Codim<cd>::Entity& e) const
    {
      DUNE_THROW( NotImplemented, "Id for codim 1 entities." );
      return 0;
    }

    //! get id of subEntity
    /*
     * We use the remove_const to extract the Type from the mutable class,
     * because the const class is not instantiated yet.
     */
    IdType subId (const typename std::remove_const<GridImp>::type::template Codim<0>::Entity& e, int i, int codim) const
    {
      assert( codim == 0 || codim == dim );
      if( codim == dim )
        return e.impl().hostEntity()->vertex( i )->info().id;
      else
        return e.impl().hostEntity()->info().id;
    }

    //! update id set in 2d
    template< int d = dim >
    std::enable_if_t< d == 2, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto& hostgrid = grid_->getHostGrid();

      // Store cell ids within cell infos
      for ( auto fc = hostgrid.finite_faces_begin(); fc != hostgrid.finite_faces_end(); ++fc)
        if( !fc->info().idWasSet )
        {
          fc->info().id = nextElementId_++;
          fc->info().idWasSet = true;
        }

      // Store vertex ids within vertex infos
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( !vh->info().idWasSet )
        {
          vh->info().id = nextVertexId_++;
          vh->info().idWasSet = true;
        }
    }

    //! update id set in 3d
    template< int d = dim >
    std::enable_if_t< d == 3, void >
    update(const GridImp* grid)
    {
      grid_ = grid;
      const auto hostgrid = grid_->getHostGrid();

      // Store cell ids within cell infos
      for ( auto fc = hostgrid.finite_cells_begin(); fc != hostgrid.finite_cells_end(); ++fc)
        if( !fc->info().idWasSet )
        {
          fc->info().id = nextElementId_++;
          fc->info().idWasSet = true;
        }

      // Store vertex ids within vertex infos
      for ( auto vh = hostgrid.finite_vertices_begin(); vh != hostgrid.finite_vertices_end(); ++vh)
        if( !vh->info().idWasSet )
        {
          vh->info().id = nextVertexId_++;
          vh->info().idWasSet = true;
        }
    }

    GridImp* grid_;
    std::size_t nextElementId_, nextVertexId_;
  };

}  // end namespace Dune

#endif
