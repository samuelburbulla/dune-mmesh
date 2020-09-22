// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_ENTITY_HH
#define DUNE_MMESH_GRID_ENTITY_HH

/** \file
 * \brief The MMeshEntity class
 */

// Dune includes
#include <dune/grid/common/grid.hh>

namespace Dune
{
  // Forward declarations
  template<int codim, int dim, class GridImp>
  class MMeshEntity;
}


// MMesh includes
#include <dune/mmesh/grid/connectedcomponent.hh>
#include <dune/mmesh/grid/polygoncutting.hh>

namespace Dune
{
  template<class GridImp>
  class MMeshLeafIntersectionIterator;

  template<class GridImp>
  class MMeshHierarchicIterator;

  // External forward declarations
  template< class Grid >
  struct HostGridAccess;


  //**********************************************************************
  //
  // --MMeshEntity
  // --Entity
  //
  /** \brief The implementation of entities in a MMesh
   *   \ingroup MMesh
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   *
   */
  template<int codim, int dim, class GridImp>
  class MMeshEntity :
    public EntityDefaultImplementation <codim,dim,GridImp,MMeshEntity>
  {
    template <class GridImp_>
    friend class MMeshLeafIndexSet;

    template <class GridImp_>
    friend class MMeshLocalIdSet;

    template <class GridImp_>
    friend class MMeshGlobalIdSet;

    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

    //! define the type used for persistent indices
    using IdType = MMeshImpl::MultiId;

  private:
    typedef typename GridImp::ctype ctype;

    // equivalent entity in the host grid
    typedef typename GridImp::template HostGridEntity<codim> HostGridEntity;

  public:
    typedef typename GridImp::template Codim<codim>::Geometry Geometry;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<codim>::EntitySeed EntitySeed;

    MMeshEntity()
      : mMesh_(nullptr)
    {}

    MMeshEntity(const GridImp* mMesh, const HostGridEntity& hostEntity)
      : hostEntity_(hostEntity)
      , mMesh_(mMesh)
    {}

    MMeshEntity(const GridImp* mMesh, HostGridEntity&& hostEntity)
      : hostEntity_(std::move(hostEntity))
      , mMesh_(mMesh)
    {}

    MMeshEntity(const MMeshEntity& original)
      : hostEntity_(original.hostEntity_)
      , mMesh_(original.mMesh_)
    {}

    MMeshEntity(MMeshEntity&& original)
      : hostEntity_(std::move(original.hostEntity_))
      , mMesh_(original.mMesh_)
    {}

    MMeshEntity& operator=(const MMeshEntity& original)
    {
      if (this != &original)
      {
        mMesh_ = original.mMesh_;
        hostEntity_ = original.hostEntity_;
      }
      return *this;
    }

    MMeshEntity& operator=(MMeshEntity&& original)
    {
      if (this != &original)
      {
        mMesh_ = original.mMesh_;
        hostEntity_ = std::move(original.hostEntity_);
      }
      return *this;
    }

    // Comparator for vertices
    template <int cc = codim>
    std::enable_if_t< cc == dim, bool >
    equals(const MMeshEntity& other) const
    {
      return hostEntity_ == other.hostEntity_;
    }

    // Comparator for edges
    template <int cc = codim>
    std::enable_if_t< cc == 1 && dim == 2, bool >
    equals(const MMeshEntity& other) const
    {
      return (hostEntity_ == other.hostEntity_)
       || ( mMesh_->getHostGrid().mirror_edge( hostEntity_ ) == other.hostEntity_ );
    }

    // Comparator for facets
    template <int cc = codim>
    std::enable_if_t< cc == 1 && dim == 3, bool >
    equals(const MMeshEntity& other) const
    {
      return (hostEntity_ == other.hostEntity_)
       || ( hostEntity_ == mMesh_->getHostGrid().mirror_facet( other.hostEntity_ ) );
    }

    // Comparator for edges in 3d
    template <int cc = codim>
    std::enable_if_t< cc == 2 && dim == 3, bool >
    equals(const MMeshEntity& other) const
    {
      const auto& v0 = hostEntity_.first->vertex( hostEntity_.second );
      const auto& v1 = hostEntity_.first->vertex( hostEntity_.third );
      const auto& w0 = other.hostEntity_.first->vertex( other.hostEntity_.second );
      const auto& w1 = other.hostEntity_.first->vertex( other.hostEntity_.third );
      return (v0 == w0 && v1 == w1) || (v0 == w1 && v1 == w0);
    }

    //! returns true if father entity exists
    bool hasFather () const {
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
      // we have a simplex grid
      int n = dim-codim+1;
      int k = dim-cc+1;

      // binomial: n over k
      int binomial=1;
      for (int i=n-k+1; i<=n; i++)
        binomial *= i;
      for (long i=2; i<=k; i++)
        binomial /= i;

      return binomial;
    }

    //! Obtain a cc dim subEntity of a codim 1 entity
    template <int cc>
    std::enable_if_t< codim == 1 && cc == dim, typename GridImp::template Codim<dim>::Entity >
    subEntity (unsigned int i) const
    {
      assert( i < subEntities( cc ) );
      // remark: the i-th edge in CGAL corresponds to the (dim-i)-th edge in DUNE,
      // but the mapping should do the right thing here
      auto& edgeIdx = hostEntity_.second;
      return MMeshEntity<cc, dim, GridImp> (
        mMesh_,
        typename GridImp::template HostGridEntity<dim> (
          hostEntity_.first->vertex( (edgeIdx+1+i)%(dim+1) )
        )
      );
    }

    //! Obtain a cc 3 subEntity of a codim 2 entity (only for 3d)
    template <int cc>
    std::enable_if_t< codim == 2 && cc == 3, typename GridImp::template Codim<3>::Entity >
    subEntity (unsigned int i) const
    {
      assert( i < subEntities( cc ) );
      auto& v0Idx = hostEntity_.second;
      auto& v1Idx = hostEntity_.third;
      // return the two vertices but do not care about any ordering here
      return MMeshEntity<cc, dim, GridImp> (
        mMesh_,
        typename GridImp::template HostGridEntity<dim> (
          hostEntity_.first->vertex( (i==0) ? v0Idx : v1Idx )
        )
      );
    }

    //! Obtain a cc dim subEntity of a codim dim entity
    template <int cc>
    std::enable_if_t< codim == dim && cc == dim, typename GridImp::template Codim<dim>::Entity >
    subEntity (unsigned int i) const
    {
      assert( i < subEntities( cc ) );
      return MMeshEntity<cc, dim, GridImp> (
        mMesh_,
        hostEntity_
      );
    }

    //! First incident element
    MMeshIncidentIterator<GridImp> incidentBegin () const {
      using Impl = typename MMeshIncidentIterator<GridImp>::Implementation;
      return MMeshIncidentIterator<GridImp>( Impl( mMesh_, hostEntity_) );
    }

    //! Last incident element
    MMeshIncidentIterator<GridImp> incidentEnd () const {
      using Impl = typename MMeshIncidentIterator<GridImp>::Implementation;
      return MMeshIncidentIterator<GridImp>( Impl( mMesh_, hostEntity_, true ) );
    }

    //! First incident facet
    MMeshIncidentFacetsIterator<GridImp> incidentFacetsBegin () const {
      using Impl = typename MMeshIncidentFacetsIterator<GridImp>::Implementation;
      return MMeshIncidentFacetsIterator<GridImp>( Impl( mMesh_, hostEntity_) );
    }

    //! Last incident facet
    MMeshIncidentFacetsIterator<GridImp> incidentFacetsEnd () const {
      using Impl = typename MMeshIncidentFacetsIterator<GridImp>::Implementation;
      return MMeshIncidentFacetsIterator<GridImp>( Impl( mMesh_, hostEntity_, true ) );
    }

    //! First incident vertex
    MMeshIncidentVerticesIterator<GridImp> incidentVerticesBegin ( bool includeInfinite ) const {
      using Impl = typename MMeshIncidentVerticesIterator<GridImp>::Implementation;
      return MMeshIncidentVerticesIterator<GridImp>( Impl( mMesh_, hostEntity_, includeInfinite) );
    }

    //! Last incident vertex
    MMeshIncidentVerticesIterator<GridImp> incidentVerticesEnd ( bool includeInfinite ) const {
      using Impl = typename MMeshIncidentVerticesIterator<GridImp>::Implementation;
      return MMeshIncidentVerticesIterator<GridImp>( Impl( mMesh_, hostEntity_, includeInfinite, true ) );
    }

    //! Return insertion level of vertex
    template <int cd = codim>
    std::enable_if_t< cd == dim, std::size_t >
    insertionLevel() const
    {
      if ( codim == dim )
        return hostEntity_->info().insertionLevel;
    }

    //! Return insertion level (maximal insertionLevel of the corresponding vertices)
    template <int cd = codim>
    std::enable_if_t< cd != dim, std::size_t >
    insertionLevel() const
    {
      std::size_t insertionLevel = 0;
      for( std::size_t i = 0; i < subEntities(dim); ++i )
        insertionLevel = std::max( insertionLevel, this->template subEntity<dim>(i).impl().insertionLevel() );
      return insertionLevel;
    }

    //! Return if vertex is part of the interface
    bool isInterface() const
    {
      static_assert( codim == dim );
      return hostEntity_->info().isInterface;
    }

    //! Return boundary flag (-1 = not set, 0 = can be removed, 1 = important for domain boundary)
    int boundaryFlag() const
    {
      static_assert( codim == dim );
      return hostEntity_->info().boundaryFlag;
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

    //! returns the host entity
    const HostGridEntity& hostEntity () const
    {
      return hostEntity_;
    }

    //! returns the host entity
    HostGridEntity& hostEntity ()
    {
      return hostEntity_;
    }

    //! returns the grid
    const GridImp& grid () const
    {
      return *mMesh_;
    }

    //! returns id computed by vertex ids
    IdType id() const
    {
      typename IdType::VT idlist( dim+1-codim );
      for( std::size_t i = 0; i < this->subEntities(dim); ++i )
        idlist[i] = this->template subEntity<dim>(i).impl().hostEntity()->info().id;
      std::sort( idlist.begin(), idlist.end() );

      return IdType( idlist );
    }

  private:
    //! the host entity of this entity
    HostGridEntity hostEntity_;

    //! the grid implementation
    const GridImp* mMesh_;
  };


  //***********************
  //
  //  --MMeshEntity
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
  class MMeshEntity<0,dim,GridImp> :
    public EntityDefaultImplementation<0,dim,GridImp, MMeshEntity>
  {
    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

    typedef Entity< 0, dim, GridImp, MMeshEntity > EntityType;

  public:
    // equivalent entity in the host grid
    typedef typename GridImp::template HostGridEntity<0> HostGridEntity;

    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    // type of global coordinate
    typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

    //! The Iterator over intersections on the leaf level
    typedef MMeshLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef MMeshHierarchicIterator<GridImp> HierarchicIterator;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

    //! The type of a ConnectedComponent
    using ConnectedComponent = MMeshConnectedComponent< GridImp >;

    //! define the type used for persistent indices
    using IdType = MMeshImpl::MultiId;

    //! define the type used for storage the vertices of a caching entity
    using VertexStorage = std::array<GlobalCoordinate, dim+1>;

    MMeshEntity()
      : isLeaf_(true), mMesh_(nullptr)
    {}

    MMeshEntity(const GridImp* mMesh, const HostGridEntity& hostEntity)
      : hostEntity_(hostEntity), isLeaf_(true)
      , mMesh_(mMesh)
    {}

    MMeshEntity(const GridImp* mMesh, HostGridEntity&& hostEntity)
      : hostEntity_(std::move(hostEntity)), isLeaf_(true)
      , mMesh_(mMesh)
    {}

    MMeshEntity(const GridImp* mMesh, const HostGridEntity& hostEntity, const IdType& id)
      : hostEntity_(hostEntity)
      , id_(id), isLeaf_(false)
      , mMesh_(mMesh)
    {}

    MMeshEntity(const GridImp* mMesh, const VertexStorage& vertex)
      : id_( /*caching id*/ IdType({42,42,42}) ), isLeaf_(false), mMesh_(mMesh), vertex_(vertex)
    {}

    MMeshEntity(const MMeshEntity& original)
      : hostEntity_(original.hostEntity_), id_(original.id_), isLeaf_(original.isLeaf_)
      , mMesh_(original.mMesh_), vertex_(original.vertex_)
    {}

    MMeshEntity(MMeshEntity&& original)
      : hostEntity_(std::move(original.hostEntity_)), id_(original.id_), isLeaf_(original.isLeaf_)
      , mMesh_(original.mMesh_), vertex_(original.vertex_)
    {}

    MMeshEntity& operator=(const MMeshEntity& original)
    {
      if (this != &original)
      {
        mMesh_ = original.mMesh_;
        hostEntity_ = original.hostEntity_;
        isLeaf_ = original.isLeaf_;
        id_ = original.id_;
        vertex_ = original.vertex_;
      }
      return *this;
    }

    MMeshEntity& operator=(MMeshEntity&& original)
    {
      if (this != &original)
      {
        mMesh_ = original.mMesh_;
        hostEntity_ = std::move(original.hostEntity_);
        isLeaf_ = original.isLeaf_;
        id_ = original.id_;
        vertex_ = original.vertex_;
      }
      return *this;
    }

    //! returns true if host entities are equal
    bool equals(const MMeshEntity& other) const
    {
      return hostEntity_ == other.hostEntity_;
    }

    //! returns true if host entities are equal
    bool operator==(const MMeshEntity& other) const
    {
      return this->equals(other);
    }

    //! returns true if host entities are equal
    bool operator<(const MMeshEntity& other) const
    {
      return hostEntity_ < other.hostEntity_;
    }

    //! returns the father entity
    MMeshEntity father () const
    {
      DUNE_THROW( InvalidStateException, "MMesh entities do no have a father, but a connectedComponent instead!" );
      return *this;
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return false;
    }

    //! returns the connected component as father entity
    const ConnectedComponent& connectedComponent () const
    {
      assert( isNew() == true );
      return mMesh_->getConnectedComponent(*this);
    }

    //! returns true if this entity is new after adaptation
    const bool isNew () const
    {
      return hostEntity_->info().isNew;
    }

    //! set if this entity is new after adaptation
    void setIsNew ( bool isNew ) const
    {
      hostEntity_->info().isNew = isNew;
    }

    //! returns true if this entity will vanish after adaptation
    const bool mightVanish () const
    {
      return hostEntity_->info().mightVanish;
    }

    //! set if this entity will vanish after adaptation
    void setWillVanish ( bool mightVanish ) const
    {
      hostEntity_->info().mightVanish = mightVanish;
    }

    //! mark entity for refine or coarse
    void mark ( int refCount ) const
    {
      hostEntity_->info().mark = refCount;
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
      if (isLeaf_)
        return Geometry( hostEntity_ );
      else
        return Geometry( this->vertex_ );
    }

    void bindFather( const EntityType& father )
    {
      father_ = &father;
    }

    //! Geometry of this entity in bounded father entity ( assumption: this \subset father )
    LocalGeometry geometryInFather() const
    {
      if constexpr( dim != 2 )
        DUNE_THROW(NotImplemented, "geometryInFather() for dim != 2");
      else
      {
        assert( father_ != nullptr );

        auto thisPoints = this->vertex_;

        if( isLeaf_ )
          for ( int i = 0; i < 3; ++i )
            thisPoints[i] = geometry().corner(i);

        std::array< GlobalCoordinate, 3 > local;
        for ( int i = 0; i < 3; ++i )
          local[i] = father_->impl().geometry().local( thisPoints[i] );

        return LocalGeometry( local );
      }
    }

    //! Return the number of subEntities of codimension cc
    unsigned int subEntities (unsigned int cc) const
    {
      // we have a simplex grid
      int n = dim+1;
      int k = dim-cc+1;

      // binomial: n over k
      int binomial=1;
      for (int i=n-k+1; i<=n; i++)
        binomial *= i;
      for (long i=2; i<=k; i++)
        binomial /= i;

      return binomial;
    }

    /** \brief Provide access to sub entity i of given codimension. Entities
     *  are numbered 0 ... subEntities(cc)-1
     */
    template <int cc>
    std::enable_if_t< cc == 0, typename GridImp::template Codim<cc>::Entity >
    subEntity (unsigned int i) const {
      assert( i < subEntities( cc ) );
      return *this;
    }

    template <int cc>
    std::enable_if_t< cc == dim, typename GridImp::template Codim<cc>::Entity >
    subEntity (unsigned int i) const {
      assert( i < subEntities( cc ) );
      if (hostEntity_ != HostGridEntity()) {
        return MMeshEntity<cc, dim, GridImp>( mMesh_, hostEntity_->vertex( i ) );
      }
      else {
        auto vh = mMesh_->getHostGrid().infinite_vertex();
        // TODO should distinguish for different i
        return MMeshEntity<cc, dim, GridImp>( mMesh_, vh );
      }
    }

    template <int cc>
    std::enable_if_t< cc == 1, typename GridImp::template Codim<cc>::Entity >
    subEntity (unsigned int i) const {
      assert( i < subEntities( cc ) );
      // remark: the i-th edge in CGAL corresponds to the (dim-i)-th edge in DUNE
      return MMeshEntity<cc, dim, GridImp> ( mMesh_, typename GridImp::template HostGridEntity<1> ( hostEntity_, dim-i ) );
    }

    template <int cc>
    std::enable_if_t< cc == 2 && dim == 3, typename GridImp::template Codim<cc>::Entity >
    subEntity (unsigned int i) const {
      assert( i < subEntities( cc ) );
      // permutate i to match the DUNE edge numbering
      static constexpr std::array<unsigned int, 6> perm = {{ 0, 4, 1, 3, 5, 2 }};
      i = perm[i];
      return MMeshEntity<cc, dim, GridImp> ( mMesh_, typename GridImp::template HostGridEntity<2> ( hostEntity_, i%4, (i+1)%4+(i>3)) );
    }

    //! First leaf intersection
    MMeshLeafIntersectionIterator<GridImp> ileafbegin () const {
      return MMeshLeafIntersectionIterator<GridImp>(
        mMesh_,
        hostEntity_);
    }

    //! Reference to one past the last leaf intersection
    MMeshLeafIntersectionIterator<GridImp> ileafend () const {
      return MMeshLeafIntersectionIterator<GridImp>(
         mMesh_,
         hostEntity_,
         true);
    }

    //! We only have one level
    MMeshLeafIntersectionIterator<GridImp> ilevelbegin () const {
      return ileafbegin();
    }

    MMeshLeafIntersectionIterator<GridImp> ilevelend () const {
      return ileafend();
    }

    //! First hierarchic entity, i.e. this entity, because we only have one level
    MMeshHierarchicIterator<GridImp> hbegin (int maxlevel) const {
      return MMeshHierarchicIterator<GridImp>(
      mMesh_,
      *this,
      maxlevel);
    }

    //! Reference to one past the last hierarchic entity
    MMeshHierarchicIterator<GridImp> hend (int maxlevel) const {
      return MMeshHierarchicIterator<GridImp>(
        mMesh_,
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

    //! return domain marker of entity
    std::size_t domainMarker() const
    {
      return hostEntity_->info().domainMarker;
    }

    //! returns the host entity
    const HostGridEntity& hostEntity () const
    {
      return hostEntity_;
    }

    //! returns the host entity
    HostGridEntity& hostEntity ()
    {
      return hostEntity_;
    }

    //! returns the host grid
    const GridImp& grid () const
    {
      return *mMesh_;
    }

    //! returns id computed by vertex ids
    IdType id() const
    {
      // cache id
      if (id_ == IdType())
      {
        typename IdType::VT idlist( dim+1 );
        for( std::size_t i = 0; i < this->subEntities(dim); ++i )
          idlist[i] = this->template subEntity<dim>(i).impl().hostEntity()->info().id;
        std::sort( idlist.begin(), idlist.end() );
        id_ = IdType( idlist );
      }

      return id_;
    }

  private:
    //! the host entity of this entity
    HostGridEntity hostEntity_;

    //! the cached id
    mutable IdType id_;

    //! return if leaf
    bool isLeaf_;

    //! the grid implementation
    const GridImp* mMesh_;

  protected:
    //! the vertices of the host entity object of this entity (for caching entity)
    VertexStorage vertex_;

    const EntityType* father_;

  }; // end of MMeshEntity codim = 0

} // namespace Dune

#include <dune/mmesh/grid/cachingentity.hh>

#endif
