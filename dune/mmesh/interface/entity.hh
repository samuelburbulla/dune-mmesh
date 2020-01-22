// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_ENTITY_HH
#define DUNE_MMESH_INTERFACE_ENTITY_HH

/** \file
 * \brief The MMeshInterfaceGridEntity class
 */

// Dune includes
#include <dune/grid/common/grid.hh>
#include <dune/mmesh/interface/incidentiterator.hh>

// CGAL includes
#include <CGAL/utility.h>

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
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   *
   */
  template<int codim, int dim, class GridImp>
  class MMeshInterfaceGridEntity :
    public EntityDefaultImplementation <codim,dim,GridImp,MMeshInterfaceGridEntity>
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
    PartitionType partitionType () const {
      return PartitionType( 0 );
    }

    //! Return the number of subEntities of codimension codim
    unsigned int subEntities (unsigned int cc) const
    {
      if( dim == 1 )
        return (cc == 0) ? 0 : 2;

      if( dim == 2 )
        return (cc == 0) ? 0 : 3;
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

    //! Return the insertion level of the vertex
    std::size_t insertionLevel() const
    {
      static_assert( codim == dim );
      return hostEntity_->info().insertionLevel;
    }

    //! First incident vertex
    auto incidentVerticesBegin ( bool includeInfinite ) const
    {
      return grid_->getMMesh().entity( hostEntity_ ).impl().incidentVerticesBegin( includeInfinite );
    }

    //! Last incident vertex
    auto incidentVerticesEnd ( bool includeInfinite ) const
    {
      return grid_->getMMesh().entity( hostEntity_ ).impl().incidentVerticesEnd( includeInfinite );
    }

    //! First incident vertex
    auto incidentInterfaceVerticesBegin () const
    {
      using Impl = typename MMeshIncidentInterfaceVerticesIterator<GridImp>::Implementation;
      return MMeshIncidentInterfaceVerticesIterator<GridImp>( Impl( grid_, hostEntity_) );
    }

    //! Last incident vertex
    auto incidentInterfaceVerticesEnd () const
    {
      using Impl = typename MMeshIncidentInterfaceVerticesIterator<GridImp>::Implementation;
      return MMeshIncidentInterfaceVerticesIterator<GridImp>( Impl( grid_, hostEntity_, true ) );
    }

    //! returns the host entity
    const MMeshInterfaceEntity& hostEntity () const
    {
      return hostEntity_;
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
  class MMeshInterfaceGridEntity<0,dim,GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp, MMeshInterfaceGridEntity>
  {
    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

  public:
    // equivalent entity in the host grid
    typedef typename GridImp::template MMeshInterfaceEntity<0> MMeshInterfaceEntity;

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

    MMeshInterfaceGridEntity()
      : grid_(nullptr)
    {}

    MMeshInterfaceGridEntity(const GridImp* grid, const MMeshInterfaceEntity& hostEntity)
      : hostEntity_(hostEntity)
      , grid_(grid)
    {
      const auto mirrored = grid_->mirrorHostEntity( hostEntity_ );
      if ( hostEntity_.first->info().index > mirrored.first->info().index )
        hostEntity_ = mirrored;
    }

    MMeshInterfaceGridEntity(const GridImp* grid, MMeshInterfaceEntity&& hostEntity)
      : hostEntity_(std::move(hostEntity))
      , grid_(grid)
    {
      const auto mirrored = grid_->mirrorHostEntity( hostEntity_ );
      if ( hostEntity_.first->info().index > mirrored.first->info().index )
        hostEntity_ = mirrored;
    }

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

    LocalGeometry geometryInFather() const {
      DUNE_THROW( InvalidStateException, "MMesh entities do no implement a geometry in father!" );
      // return LocalGeometry( std::array<Vertex, 3>() ); // dummy
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
      return hostEntity_->info().mark;
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
    PartitionType partitionType () const {
      return PartitionType::InteriorEntity; /* dummy */
    }

    //! Geometry of this entity
    Geometry geometry () const
    {
      return Geometry( hostEntity_ );
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
      return MMeshInterfaceGridEntity<cc, dim, GridImp>(
        grid_, hostEntity_.first->vertex( (hostEntity_.second+i+1)%(dim+2) )
      );
    }

    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template <int cc>
    std::enable_if_t< cc == 1 && dim == 2, typename GridImp::template Codim<cc>::Entity >
    subEntity (std::size_t i) const {
      assert( i < subEntities( cc ) );
      const int j = hostEntity_.second;
      const int v1 = (j+i+1)%4;
      const int v2 = (i==2) ? ((j+1)%4) : ((j+i+2)%4);

      return MMeshInterfaceGridEntity<cc, dim, GridImp>(
        grid_, CGAL::Triple<decltype(hostEntity_.first), int, int>( hostEntity_.first, v1, v2 )
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
      return true;
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

  private:
    //! the host entity of this entity
    MMeshInterfaceEntity hostEntity_;

    //! the grid implementation
    const GridImp* grid_;

  }; // end of MMeshInterfaceGridEntity codim = 0

} // namespace Dune

#endif
