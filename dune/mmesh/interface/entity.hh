// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_ENTITY_HH
#define DUNE_MMESH_INTERFACE_ENTITY_HH

/** \file
 * \brief The MMeshInterfaceGridEntity class
 */

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/mmesh/interface/common.hh>
#include <dune/mmesh/interface/incidentiterator.hh>

namespace Dune
{
  // External forward declarations
  template<class Grid>
  struct HostGridAccess;


  //**********************************************************************
  //
  // --MMeshInterfaceGridEntity
  // --Entity
  //
  /** \brief The implementation of entities in a MMesh interface grid
   *   \ingroup MMeshInterfaceGrid
   */
  template<int codim, int dim, class GridImp>
  class MMeshInterfaceGridEntity
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   : public EntityDefaultImplementation <codim,dim,GridImp,MMeshInterfaceGridEntity>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    template <class GridImp_>
    friend class MMeshInterfaceGridLeafIndexSet;

    template <class GridImp_>
    friend class MMeshInterfaceGridLocalIdSet;

    template <class GridImp_>
    friend class MMeshInterfaceGridGlobalIdSet;

    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

  private:
    typedef typename GridImp::ctype ctype;

    // equivalent entity in the host grid
    typedef typename GridImp::template MMeshInterfaceEntity<codim> MMeshInterfaceEntity;

    //! define the type used for persistent indices
    using IdType = MMeshImpl::MultiId;

  public:
    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<codim>::EntitySeed EntitySeed;

    MMeshInterfaceGridEntity()
      : grid_(nullptr)
    {}

    MMeshInterfaceGridEntity(const GridImp* grid, const MMeshInterfaceEntity& hostEntity)
      : hostEntity_(hostEntity)
      , grid_(grid)
    {}

    MMeshInterfaceGridEntity(const GridImp* grid, MMeshInterfaceEntity&& hostEntity)
      : hostEntity_(std::move(hostEntity))
      , grid_(grid)
    {}

    MMeshInterfaceGridEntity(const MMeshInterfaceGridEntity& original)
      : hostEntity_(original.hostEntity_)
      , grid_(original.grid_)
    {}

    MMeshInterfaceGridEntity(MMeshInterfaceGridEntity&& original)
      : hostEntity_(std::move(original.hostEntity_))
      , grid_(original.grid_)
    {}

    MMeshInterfaceGridEntity& operator=(const MMeshInterfaceGridEntity& original)
    {
      if (this != &original)
      {
        grid_ = original.grid_;
        hostEntity_ = original.hostEntity_;
      }
      return *this;
    }

    MMeshInterfaceGridEntity& operator=(MMeshInterfaceGridEntity&& original)
    {
      if (this != &original)
      {
        grid_ = original.grid_;
        hostEntity_ = std::move(original.hostEntity_);
      }
      return *this;
    }

    // Comparator for vertices
    template <int cc = codim>
    std::enable_if_t< cc == dim, bool >
    equals(const MMeshInterfaceGridEntity& other) const
    {
      return hostEntity_ == other.hostEntity_;
    }

    // Comparator for edges
    template <int cc = codim>
    std::enable_if_t< cc == 1 && dim == 2, bool >
    equals(const MMeshInterfaceGridEntity& other) const
    {
      const auto& vh1 = hostEntity_.first->vertex(hostEntity_.second);
      const auto& vh2 = hostEntity_.first->vertex(hostEntity_.third);
      const auto& hostOther_ = other.hostEntity_;
      const auto& vo1 = hostOther_.first->vertex(hostOther_.second);
      const auto& vo2 = hostOther_.first->vertex(hostOther_.third);

      return ( (vh1 == vo1) && ( vh2 == vo2 ) ) || ( (vh1 == vo2) && ( vh2 == vo1 ) );
    }

    //! returns true if connected component entity exists
    bool hasConnectedComponent () const {
      return false;
    }

    //! Return entity seed
    EntitySeed seed () const
    {
      return EntitySeed( hostEntity_ );
    }

    //! level of this element
    int level () const {
      // we only have one level
      return 0;
    }

    //! The partition type for parallel computing
    PartitionType partitionType () const
    {
      return grid().getMMesh().partitionHelper().partitionType( grid().entity(hostEntity_) );
    }

    //! Return the number of subEntities of codimension codim
    unsigned int subEntities (unsigned int cc) const
    {
      if( dim == 1 )
        return (cc == 0) ? 0 : 2;

      if( dim == 2 )
        return (cc == 0) ? 0 : 3;
    }

    /** \brief Provide access to sub entity i for cc == dim
     */
    template <int cc>
    std::enable_if_t< cc == codim, typename GridImp::template Codim<cc>::Entity >
    subEntity (std::size_t i) const {
      return *this;
    }

    /** \brief Provide access to sub entity i for cc == dim-1
     */
    template <int cc>
    std::enable_if_t< cc == codim+1, typename GridImp::template Codim<cc>::Entity >
    subEntity (std::size_t i) const {
      DUNE_THROW(NotImplemented, "subEntity<1> for codim 1 entity");
      return typename GridImp::template Codim<cc>::Entity();
    }

    //! geometry of this entity
    Geometry geometry () const
    {
      return Geometry( hostEntity_ );
    }

    //! returns the geometry type
    GeometryType type () const
    {
      return GeometryTypes::simplex(dim-codim);
    }

    //! returns that entity is part of the interface
    bool isInterface() const
    {
      return true;
    }

    //! Return if this vertex is a tip
    bool isTip() const
    {
      if constexpr ( codim == dim )
      {
        if constexpr ( dim == 2 )
          return false;
        else
        {
          const auto& hostgrid = grid().getHostGrid();

          int count = 0;
          auto circulator = hostgrid.incident_edges( hostEntity_ );
          for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
          {
            // at boundary
            if (hostgrid.is_infinite(circulator))
              return false;

             if (grid().isInterface(*circulator))
              count++;
          }

          return (count == 1);
        }
      }
      else
        DUNE_THROW( NotImplemented, "isTip() for codim != dim" );
    };

    //! Return boundary flag (-1 = not set, 0 = can be removed, 1 = important for domain boundary)
    int boundaryFlag() const
    {
      if constexpr ( codim == dim )
        return hostEntity_->info().boundaryFlag;
      else
        DUNE_THROW( NotImplemented, "boundaryFlag() for codim != dim" );
    }

    //! Return the insertion level of the vertex
    std::size_t insertionLevel() const
    {
      if constexpr ( codim == dim )
        return hostEntity_->info().insertionLevel;
      else
        DUNE_THROW( NotImplemented, "boundaryFlag() for codim != dim" );
    }

    //! First incident vertex
    auto incidentVerticesBegin ( bool includeInfinite ) const
    {
      if constexpr ( codim == dim )
      {
        using Impl = typename MMeshIncidentVerticesIterator<GridImp>::Implementation;
        return MMeshIncidentVerticesIterator<GridImp>( Impl( grid_, hostEntity_, includeInfinite) );
      }
      else
        DUNE_THROW( NotImplemented, "incidentVerticesBegin() for codim != dim" );
    }

    //! Last incident vertex
    auto incidentVerticesEnd ( bool includeInfinite ) const
    {
      if constexpr ( codim == dim )
      {
        using Impl = typename MMeshIncidentVerticesIterator<GridImp>::Implementation;
        return MMeshIncidentVerticesIterator<GridImp>( Impl( grid_, hostEntity_, includeInfinite, true ) );
      }
      else
        DUNE_THROW( NotImplemented, "incidentVerticesEnd() for codim != dim" );
    }

    //! First incident vertex
    auto incidentInterfaceVerticesBegin () const
    {
      if constexpr ( codim == dim )
      {
        using Impl = typename MMeshIncidentInterfaceVerticesIterator<GridImp>::Implementation;
        return MMeshIncidentInterfaceVerticesIterator<GridImp>( Impl( grid_, hostEntity_ ) );
      }
      else
        DUNE_THROW( NotImplemented, "incidentInterfaceVerticesBegin() for codim != dim" );
    }

    //! Last incident vertex
    auto incidentInterfaceVerticesEnd () const
    {
      if constexpr ( codim == dim )
      {
        using Impl = typename MMeshIncidentInterfaceVerticesIterator<GridImp>::Implementation;
        return MMeshIncidentInterfaceVerticesIterator<GridImp>( Impl( grid_, hostEntity_, true ) );
      }
      else
        DUNE_THROW( NotImplemented, "incidentInterfaceVerticesEnd() for codim != dim" );
    }

    //! First incident element
    auto incidentInterfaceElementsBegin () const
    {
      if constexpr ( codim == dim )
      {
        using Impl = typename MMeshIncidentInterfaceElementsIterator<GridImp>::Implementation;
        return MMeshIncidentInterfaceElementsIterator<GridImp>( Impl( grid_, hostEntity_ ) );
      }
      else if constexpr ( codim == dim-1 )
      {
        using Impl = typename MMeshEdgeIncidentInterfaceElementsIterator<GridImp>::Implementation;
        return MMeshEdgeIncidentInterfaceElementsIterator<GridImp>( Impl( grid_, hostEntity_ ) );
      }
      else
        DUNE_THROW( NotImplemented, "incidentInterfaceElementsBegin() for codim <= dim-1" );
    }

    //! Last incident element
    auto incidentInterfaceElementsEnd () const
    {
      if constexpr ( codim == dim )
      {
        using Impl = typename MMeshIncidentInterfaceElementsIterator<GridImp>::Implementation;
        return MMeshIncidentInterfaceElementsIterator<GridImp>( Impl( grid_, hostEntity_, true ) );
      }
      else if constexpr ( codim == dim-1 )
      {
        using Impl = typename MMeshEdgeIncidentInterfaceElementsIterator<GridImp>::Implementation;
        return MMeshEdgeIncidentInterfaceElementsIterator<GridImp>( Impl( grid_, hostEntity_, true ) );
      }
      else
        DUNE_THROW( NotImplemented, "incidentInterfaceElementsEnd() for codim != dim" );
    }

    //! Return reference to the host entity
    const MMeshInterfaceEntity& hostEntity () const
    {
      return hostEntity_;
    }

    //! Return id
    IdType id() const
    {
      if constexpr ( codim == dim )
        return hostEntity_->info().id;
      else
        return grid().globalIdSet().id( *this );
    }

    //! returns the grid
    const GridImp& grid () const
    {
      return *grid_;
    }

  private:
    //! the host entity of this entity
    MMeshInterfaceEntity hostEntity_;

    //! the grid implementation
    const GridImp* grid_;
  };


  //***********************
  //
  //  --MMeshInterfaceGridEntity
  //
  //***********************
  /** \brief Specialization for codim-0-entities.
   * \ingroup MMesh
   *
   * This class embodies the topological parts of elements of the grid.
   * It has an extended interface compared to the general entity class.
   * For example, Entities of codimension 0 allow to visit all neighbors.
   */
  template<int dim, class GridImp>
  class MMeshInterfaceGridEntity<0,dim,GridImp>
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   : public EntityDefaultImplementation<0,dim,GridImp, MMeshInterfaceGridEntity>
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

    typedef Entity< 0, dim, GridImp, MMeshInterfaceGridEntity > EntityType;

  public:
    // equivalent entity in the host grid
    typedef typename GridImp::template MMeshInterfaceEntity<0> MMeshInterfaceEntity;

    typedef typename GridImp::template MMeshInterfaceEntity<dim> MMeshInterfaceVertex;

    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    //! The Iterator over intersections on the leaf level
    typedef MMeshInterfaceGridLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef MMeshHierarchicIterator<GridImp> HierarchicIterator;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

    //! The type of the connected component
    typedef typename GridImp::ConnectedComponent ConnectedComponent;

    //! define the type used for persistent indices
    using IdType = MMeshImpl::MultiId;

    // type of global coordinate
    typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

    // type of local coordinate
    typedef typename Geometry::LocalCoordinate LocalCoordinate;

    //! define the type used for storage the vertices of a caching entity
    using VertexStorage = std::array<GlobalCoordinate, dim+1>;

    MMeshInterfaceGridEntity()
      : grid_(nullptr), isLeaf_(true)
    {}

    MMeshInterfaceGridEntity(const GridImp* grid, const MMeshInterfaceEntity& hostEntity)
      : hostEntity_(hostEntity), grid_(grid), isLeaf_(true)
    {
      if (grid->canBeMirrored(hostEntity_))
      {
        const auto mirrored = grid_->mirrorHostEntity( hostEntity_ );
        if ( hostEntity_.first->info().insertionIndex > mirrored.first->info().insertionIndex )
          hostEntity_ = mirrored;
      }
    }

    MMeshInterfaceGridEntity(const GridImp* grid, MMeshInterfaceEntity&& hostEntity)
      : hostEntity_(std::move(hostEntity)), grid_(grid), isLeaf_(true)
    {
      if (grid->canBeMirrored(hostEntity_))
      {
        const auto mirrored = grid_->mirrorHostEntity( hostEntity_ );
        if ( hostEntity_.first->info().insertionIndex > mirrored.first->info().insertionIndex )
          hostEntity_ = mirrored;
      }
    }

    MMeshInterfaceGridEntity(const GridImp* grid, const MMeshInterfaceEntity& hostEntity, const IdType& id)
      : hostEntity_(hostEntity), id_(id), grid_(grid), isLeaf_(false)
    {}

    MMeshInterfaceGridEntity(const GridImp* grid, const VertexStorage& vertex)
      : id_( /*caching id*/ IdType({ std::size_t(-3), std::size_t(-2) }) ), grid_(grid), isLeaf_(false), vertex_(vertex)
    {}

    MMeshInterfaceGridEntity(const MMeshInterfaceGridEntity& original)
      : hostEntity_(original.hostEntity_), id_(original.id_)
      , grid_(original.grid_), isLeaf_(original.isLeaf_), vertex_(original.vertex_)
    {}

    MMeshInterfaceGridEntity(MMeshInterfaceGridEntity&& original)
      : hostEntity_(std::move(original.hostEntity_)), id_(original.id_)
      , grid_(original.grid_), isLeaf_(original.isLeaf_), vertex_(original.vertex_)
    {}

    MMeshInterfaceGridEntity& operator=(const MMeshInterfaceGridEntity& original)
    {
      if (this != &original)
      {
        grid_ = original.grid_;
        hostEntity_ = original.hostEntity_;
        isLeaf_ = original.isLeaf_;
        id_ = original.id_;
        vertex_ = original.vertex_;
      }
      return *this;
    }

    MMeshInterfaceGridEntity& operator=(MMeshInterfaceGridEntity&& original)
    {
      if (this != &original)
      {
        grid_ = original.grid_;
        hostEntity_ = std::move(original.hostEntity_);
        isLeaf_ = original.isLeaf_;
        id_ = original.id_;
        vertex_ = original.vertex_;
      }
      return *this;
    }

    // Comparator for edges
    template <int d = dim>
    std::enable_if_t< d == 1, bool >
    equals(const MMeshInterfaceGridEntity& other) const
    {
      return hostEntity_ == other.hostEntity_;
    }

    // Comparator for facets
    template <int d = dim>
    std::enable_if_t< d == 2, bool >
    equals(const MMeshInterfaceGridEntity& other) const
    {
      return hostEntity_ == other.hostEntity_;
    }

    //! returns true if host entities are equal
    bool operator==(const MMeshInterfaceGridEntity& other) const
    {
      return this->equals(other);
    }

    //! returns true if host entities are equal
    bool operator<(const MMeshInterfaceGridEntity& other) const
    {
      return hostEntity_ < other.hostEntity_;
    }

    //! returns the connected component
    const ConnectedComponent& connectedComponent () const
    {
      return grid_->getConnectedComponent( *this );
    }

    //! returns true if a connected component exists
    bool hasConnectedComponent () const
    {
      return grid_->hasConnectedComponent( *this );
    }

    //! returns the father entity
    auto father () const
    {
      DUNE_THROW( InvalidStateException, "MMesh entities do no have a father, but a connectedComponent instead!" );
      return *this;
    }

    //! returns true if father entity exists
    bool hasFather () const {
      return false;
    }

    void bindFather( const EntityType& father ) const
    {
      father_ = &father;
    }

    //! Geometry of this entity in bounded father entity ( assumption: this \subset father )
    LocalGeometry geometryInFather() const
    {
      if constexpr( dim != 1 )
        DUNE_THROW(NotImplemented, "geometryInFather() for dim != 1");
      else
      {
        assert( father_ != nullptr );

        auto thisPoints = this->vertex_;

        if( isLeaf_ )
          for ( int i = 0; i < 2; ++i )
            thisPoints[i] = geometry().corner(i);

        std::array< LocalCoordinate, 2 > local;
        for ( int i = 0; i < 2; ++i )
          local[i] = father_->impl().geometry().local( thisPoints[i] );

        return LocalGeometry( local );
      }
    }

    //! returns true if this entity is new after adaptation
    const bool isNew () const
    {
      return hasConnectedComponent();
    }

    //! set if this entity is new after adaptation
    void setIsNew ( bool isNew ) const
    {}

    //! returns true if this entity will vanish after adaptation
    const bool mightVanish () const
    {
      return hostEntity_.first->info().mightVanish;
    }

    //! set if this entity will vanish after adaptation
    void setWillVanish ( bool mightVanish ) const
    {
      hostEntity_.first->info().mightVanish = mightVanish;
    }

    //! mark entity for refine or coarse
    void mark ( int refCount ) const
    {
      hostEntity_.first->info().mark = refCount;
    }

    //! get mark of entity
    int getMark () const
    {
      return hostEntity_.first->info().mark;
    }

    //! Create EntitySeed
    EntitySeed seed () const
    {
      return EntitySeed(hostEntity_);
    }

    //! Level of this element
    int level () const
    {
      // we only have one level
      return 0;
    }

    //! The partition type for parallel computing
    PartitionType partitionType () const
    {
      return grid().getMMesh().partitionHelper().partitionType( grid().entity(hostEntity_) );
    }

    //! Geometry of this entity
    Geometry geometry () const
    {
      if (isLeaf_)
        return Geometry( hostEntity_ );
      else
        return Geometry( this->vertex_ );
    }

    //! Return the number of subEntities of codimension cc
    std::size_t subEntities (std::size_t cc) const
    {
      if( dim == 1 )
        return (cc == 0) ? 1 : 2;

      if( dim == 2 )
        return (cc == 0) ? 1 : 3;
    }

    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template <int cc>
    std::enable_if_t< cc == 0, typename GridImp::template Codim<cc>::Entity >
    subEntity (std::size_t i) const {
      return *this;
    }

    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template <int cc>
    std::enable_if_t< cc == dim, typename GridImp::template Codim<cc>::Entity >
    subEntity (std::size_t i) const {
      assert( i < subEntities( cc ) );

      auto cgalIndex = MMeshInterfaceImpl::computeCGALIndices<MMeshInterfaceEntity, dim>( hostEntity_ );
      return MMeshInterfaceGridEntity<cc, dim, GridImp>(
        grid_, hostEntity_.first->vertex( cgalIndex[i] )
      );
    }

    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template <int cc>
    std::enable_if_t< cc == 1 && dim == 2, typename GridImp::template Codim<cc>::Entity >
    subEntity (std::size_t i) const {
      assert( i < subEntities( cc ) );

      auto cgalIndex = MMeshInterfaceImpl::computeCGALIndices<MMeshInterfaceEntity, dim>( hostEntity_ );
      auto cell = hostEntity_.first;

      int v1 = cgalIndex[i==2 ? 1 : 0];
      int v2 = cgalIndex[i==0 ? 1 : 2];

      return MMeshInterfaceGridEntity<cc, dim, GridImp>(
        grid_, CGAL::Triple<decltype(cell), int, int>( cell, v1, v2 )
      );
    }

    //! First leaf intersection
    MMeshInterfaceGridLeafIntersectionIterator<GridImp> ileafbegin () const {
      return MMeshInterfaceGridLeafIntersectionIterator<GridImp>(
      grid_,
      hostEntity_);
    }

    //! Reference to one past the last leaf intersection
    MMeshInterfaceGridLeafIntersectionIterator<GridImp> ileafend () const {
      return MMeshInterfaceGridLeafIntersectionIterator<GridImp>(
         grid_,
         hostEntity_,
         true);
    }

    //! We only have one level
    MMeshInterfaceGridLeafIntersectionIterator<GridImp> ilevelbegin () const {
      return ileafbegin();
    }

    MMeshInterfaceGridLeafIntersectionIterator<GridImp> ilevelend () const {
      return ileafend();
    }

    //! First hierarchic entity, i.e. this entity, because we only have one level
    MMeshInterfaceGridHierarchicIterator<GridImp> hbegin (int maxlevel) const {
      return MMeshInterfaceGridHierarchicIterator<GridImp>(
      grid_,
      *this,
      maxlevel);
    }

    //! Reference to one past the last hierarchic entity
    MMeshInterfaceGridHierarchicIterator<GridImp> hend (int maxlevel) const {
      return MMeshInterfaceGridHierarchicIterator<GridImp>(
        grid_,
         *this,
         maxlevel,
         true);
    }

    //! returns true if Entity has no children
    bool isLeaf() const {
      return isLeaf_ && !isNew();
    }

    //! returns if grid was refined
    bool wasRefined () const
    {
      return false;
    }

    //! returns if grid might be coarsened
    bool mightBeCoarsened () const
    {
      return false;
    }

    //! returns the geometry type
    GeometryType type () const
    {
      return GeometryTypes::simplex(dim);
    }

    //! returns the host entity
    const MMeshInterfaceEntity& hostEntity () const
    {
      return hostEntity_;
    }

    //! returns the host entity
    const GridImp& grid () const
    {
      return *grid_;
    }

    //! return cached id
    IdType id() const
    {
      // cache id
      if (id_ == IdType())
      {
        typename IdType::VT idlist( dim+1 );
        for( std::size_t i = 0; i < this->subEntities(dim); ++i )
          if (grid_->canBeMirrored(hostEntity_))
            idlist[i] = this->subEntity<dim>(i).impl().hostEntity()->info().id;
          else
            idlist[i] = -i;
        std::sort( idlist.begin(), idlist.end() );
        id_ = IdType( idlist );
      }

      return id_;
    }

  private:
    //! the host entity of this entity
    MMeshInterfaceEntity hostEntity_;

    //! the cached id
    mutable IdType id_;

    //! the grid implementation
    const GridImp* grid_;

    //! return if leaf
    bool isLeaf_;

  protected:
    //! the vertices of the host entity object of this entity (for caching entity)
    VertexStorage vertex_;

    mutable const EntityType* father_;

  }; // end of MMeshInterfaceGridEntity codim = 0

} // namespace Dune

#endif
