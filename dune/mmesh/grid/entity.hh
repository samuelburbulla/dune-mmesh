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
       || ( hostEntity_ == mMesh_->getHostGrid().mirror_edge( other.hostEntity_ ) );
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

    //! Obtain a cc 2 subEntity of a codim 1 entity
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

  public:
    // equivalent entity in the host grid
    typedef typename GridImp::template HostGridEntity<0> HostGridEntity;

    typedef typename GridImp::template Codim<0>::Geometry Geometry;

    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    //! The Iterator over intersections on the leaf level
    typedef MMeshLeafIntersectionIterator<GridImp> LeafIntersectionIterator;

    //! Iterator over descendants of the entity
    typedef MMeshHierarchicIterator<GridImp> HierarchicIterator;

    //! The type of the EntitySeed interface class
    typedef typename GridImp::template Codim<0>::EntitySeed EntitySeed;

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
      DUNE_THROW( InvalidStateException, "MMesh entities do no have a father!" );
      return *this;
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return false;
    }

    //! returns true if this entity is new after adaptation
    const bool isNew () const
    {
      return false;
    }

    //! returns true if this entity will vanish after adaptation
    const bool mightVanish () const
    {
      return false;
    }

    //! set if this entity will vanish after adaptation
    void setWillVanish ( bool mightVanish ) const {}

    //! mark entity for refine or coarse
    void mark ( int refCount ) const {}

    //! get mark of entity
    int getMark () const
    {
      return 0;
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
    std::enable_if_t< cc == dim, typename GridImp::template Codim<cc>::Entity >
    subEntity (unsigned int i) const {
      assert( i < subEntities( cc ) );
      return MMeshEntity<cc, dim, GridImp>( mMesh_, hostEntity_->vertex( i ) );
    }

    template <int cc>
    std::enable_if_t< cc == 1, typename GridImp::template Codim<cc>::Entity >
    subEntity (unsigned int i) const {
      assert( i < subEntities( cc ) );
      // remark: the i-th edge in CGAL corresponds to the (dim-i)-th edge in DUNE
      return MMeshEntity<cc, dim, GridImp> ( mMesh_, typename GridImp::template HostGridEntity<1> ( hostEntity_, dim-i ) );
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
    const HostGridEntity& hostEntity () const
    {
      return hostEntity_;
    }

  private:
    //! the host entity of this entity
    HostGridEntity hostEntity_;

    //! the grid implementation
    const GridImp* mMesh_;

  }; // end of MMeshEntity codim = 0

} // namespace Dune

#endif
