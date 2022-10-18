// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_GRID_HH
#define DUNE_MMESH_INTERFACE_GRID_HH

/** \file
 * \brief The MMeshInterfaceGrid class
 */

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include <dune/common/deprecated.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/version.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/mmesh/grid/multiid.hh>

// The components of the MMesh interface
#include <dune/mmesh/grid/mmesh.hh>
#include "entity.hh"
#include "entityseed.hh"
#include "connectedcomponent.hh"
#include "dgfparser.hh"
#include "geometry.hh"
#include "gridfactory.hh"
#include "hierarchiciterator.hh"
#include "indexsets.hh"
#include "intersectioniterator.hh"
#include "intersections.hh"
#include "leafiterator.hh"
#include "traits.hh"

namespace Dune
{

#if HAVE_MPI
  using Comm = MPI_Comm;
#else
  using Comm = No_Comm;
#endif

  // MMeshInterfaceGrid family
  template<class MMesh>
  struct MMeshInterfaceGridFamily
  {

  public:

    typedef GridTraits<
        MMesh::dimension-1, // grid dimension
        MMesh::dimension, // world dimension
        MMeshInterfaceGrid<MMesh>,
        MMeshInterfaceGridGeometry,
        MMeshInterfaceGridEntity,
        MMeshInterfaceGridLeafIterator, // LevelIterator
        MMeshInterfaceGridLeafIntersection,
        MMeshInterfaceGridLeafIntersection, // LevelIntersection
        MMeshInterfaceGridLeafIntersectionIterator,
        MMeshInterfaceGridLeafIntersectionIterator, // LevelIntersectionIterator
        MMeshInterfaceGridHierarchicIterator,
        MMeshInterfaceGridLeafIterator,
        MMeshInterfaceGridLeafIndexSet< const MMeshInterfaceGrid<MMesh> >, // LevelIndexSet
        MMeshInterfaceGridLeafIndexSet< const MMeshInterfaceGrid<MMesh> >,
        MMeshInterfaceGridGlobalIdSet< const MMeshInterfaceGrid<MMesh> >,
        MMeshImpl::MultiId, // GlobalIdSet::IdType,
        MMeshInterfaceGridGlobalIdSet< const MMeshInterfaceGrid<MMesh> >, // LocalIdSet
        MMeshImpl::MultiId, // LocalIdSet::IdType,
        Communication< Comm >,
        DefaultLevelGridViewTraits,
        DefaultLeafGridViewTraits,
        MMeshInterfaceGridEntitySeed
        > Traits;
  };

  //**********************************************************************
  //
  // --MMeshInterfaceGrid class
  //
  //************************************************************************
  /*!
   * \brief Provides a DUNE grid interface class for the interface of a MMesh interface grid
   * \ingroup GridImplementations
   * \ingroup MMeshInterfaceGrid
   */
  template <class MMesh>
  class MMeshInterfaceGrid
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   : public GridDefaultImplementation< MMesh::dimension-1, MMesh::dimension,
                                       typename MMesh::FieldType,
                                       MMeshInterfaceGridFamily<MMesh> >
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
  public:
    static constexpr int dimension = MMesh::dimension-1;
    static constexpr int dimensionworld = MMesh::dimension;

    using FieldType = typename MMesh::FieldType;
    using IdType = typename MMesh::IdType;

    using BoundarySegments = std::unordered_map< IdType, std::size_t >;

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! the Traits
    typedef MMeshInterfaceGridFamily<MMesh> GridFamily;
    typedef typename GridFamily::Traits Traits;

    //! the grid implementation
    using GridImp = typename GridFamily::Traits::Grid;

    //! the unique pointer to the grid
    typedef std::unique_ptr< GridImp > GridPtrType;

    //! the underlying hostgrid
    using HostGridType = typename MMesh::HostGridType;

    //! the underlying mmesh
    using MMeshType = MMesh;

    //! the leaf iterator
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;

    //! The type used to store coordinates, inherited from the MMesh
    typedef FieldType ctype;

    //! The type used for coordinates
    typedef Dune::FieldVector<ctype, dimensionworld> GlobalCoordinate;

    //! The type of the underlying entities
    template<int cd>
    using MMeshInterfaceEntity = typename HostGridEntityChooser_<HostGridType, dimensionworld, cd+1>::type;

    //! The type of the underlying vertex handle
    using VertexHandle = MMeshInterfaceEntity<dimension>;

    //! The type of the underlying edge handle
    using EdgeHandle = MMeshInterfaceEntity<dimension-1>;

    //! The type of the element output
    using ElementOutput = std::list<MMeshInterfaceEntity<0>>;

    //! The type of a codim 0 entity
    using Entity = typename Traits::template Codim<0>::Entity;

    //! The type of a vertex
    using Vertex = typename Traits::template Codim<dimension>::Entity;

    //! The type of a connected component
    using ConnectedComponent = MMeshInterfaceConnectedComponent<const GridImp>;

    //! The type of the default remeshing indicator
    using RemeshingIndicator = RatioIndicator<GridImp>;

    //! Constructor
    explicit MMeshInterfaceGrid(MMesh* mMesh, BoundarySegments boundarySegments = {})
     : mMesh_(mMesh),
       boundarySegments_(boundarySegments),
       #ifdef HAVE_MPI
       comm_( MPIHelper::getCommunicator() )
      #endif
    {
      leafIndexSet_ = std::make_unique<MMeshInterfaceGridLeafIndexSet<const GridImp>>( this );
      globalIdSet_ = std::make_unique<MMeshInterfaceGridGlobalIdSet<const GridImp>>( this );
    }

    /** \brief Return maximum level defined in this grid.
     *
     * Levels are numbered 0 ... maxlevel with 0 the coarsest level.
     */
    int maxLevel() const {
      return 0;
    }

    /** \brief Number of grid entities per level and codim
     */
    int size (int level, int codim) const {
      // we only have one level
      assert(level == 0);
      return size(codim);
    }

    /** \brief returns the number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const {
      return localBoundarySegments_.size();
    }

    const BoundarySegments& boundarySegments() const
    {
      return localBoundarySegments_;
    }

    void addBoundarySegment ( const std::vector< std::size_t >& ids, std::size_t bndSegIdx )
    {
      boundarySegments_[ ids ] = bndSegIdx;
    }

    void setBoundarySegments( const BoundarySegments boundarySegments )
    {
      boundarySegments_ = boundarySegments;
    }

    //! number of leaf entities per codim in this process
    int size (int codim) const {
      return leafIndexSet().size(codim);
    }


    //! number of entities per level, codim and geometry type in this process
    int size (int level, GeometryType type) const {
      // we only have one level
      assert(level == 0);
      return size(type);
    }


    //! number of leaf entities per codim and geometry type in this process
    int size (GeometryType type) const
    {
      return leafIndexSet().size(type);
    }

    /** \brief Access to the GlobalIdSet */
    const typename Traits::GlobalIdSet& globalIdSet() const {
      return *globalIdSet_;
    }


    /** \brief Access to the LocalIdSet */
    const typename Traits::LocalIdSet& localIdSet() const {
      return *globalIdSet_;
    }


    /** \brief Access to the LevelIndexSets */
    const MMeshInterfaceGridLeafIndexSet<const GridImp>& levelIndexSet(int level) const
    {
      if (level != 0)
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
      return *leafIndexSet_;
    }

    /** \brief Access to the LeafIndexSet */
    const MMeshInterfaceGridLeafIndexSet<const GridImp>& leafIndexSet() const
    {
      return *leafIndexSet_;
    }

    /** \brief Create Entity from EntitySeed */
    template < class EntitySeed >
    typename Traits::template Codim<EntitySeed::codimension>::Entity
    entity(const EntitySeed& seed) const
    {
      typedef MMeshInterfaceGridEntity<
        EntitySeed::codimension,
        dimension,
        const typename Traits::Grid
        > EntityImp;

      auto hostEntity = seed.impl().hostEntity();
      assert( hostEntity != decltype(hostEntity)() );
      return EntityImp(this, hostEntity);
    }

    /** \brief Create Entity from a host entity */
    typename Traits::template Codim<dimension>::Entity entity(const MMeshInterfaceEntity<dimension>& hostEntity) const
    {
      return entity( typename Traits::template Codim<dimension>::EntitySeed( hostEntity ) );
    }

    /** \brief Create codim 0 entity from a host entity */
    typename Traits::template Codim<0>::Entity entity(const MMeshInterfaceEntity<0>& hostEntity) const
    {
      return entity( typename Traits::template Codim<0>::EntitySeed( hostEntity ) );
    }

    /** \brief Create codim 0 entity from a host entity */
    template<int dimw = dimension+1>
    std::enable_if_t<dimw == 3, typename Traits::template Codim<1>::Entity> entity(const MMeshInterfaceEntity<1>& hostEntity) const
    {
      return entity( typename Traits::template Codim<1>::EntitySeed( hostEntity ) );
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return MMeshInterfaceGridLeafIterator<codim,All_Partition, const GridImp>(this);
    }


    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return MMeshInterfaceGridLeafIterator<codim,All_Partition, const GridImp>(this, true);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      return MMeshInterfaceGridLeafIterator<codim,PiType, const GridImp>(this);
    }


    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
      return MMeshInterfaceGridLeafIterator<codim,PiType, const GridImp>(this, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return MMeshInterfaceGridLeafIterator<codim,All_Partition, const GridImp>(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return MMeshInterfaceGridLeafIterator<codim,All_Partition, const GridImp>(this, true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      return MMeshInterfaceGridLeafIterator<codim,PiType, const GridImp>(this);
    }


    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      return MMeshInterfaceGridLeafIterator<codim,PiType, const GridImp>(this, true);
    }

    /** \brief Global refine
     *
     * Marks all elements for refinement and adapts the grid
     */
    void globalRefine(int steps = 1)
    {
      for(int i = 0; i < steps; ++i)
      {
        // mark all elements
        for (const auto& element : elements(this->leafGridView()))
          mark( 1, element );

        preAdapt();
        adapt();
        postAdapt();
      }
    }


    /** \brief Mark entity for refinement
     *
     * This only works for entities of codim 0.
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e) const
    {
      if ( refCount != 0 )
      {
        mark_[ globalIdSet_->id( e ) ] = refCount;
        return true;
      }
      return false;
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark
     */
    int getMark(const typename Traits::template Codim<0>::Entity & e) const
    {
      auto entry = mark_.find( globalIdSet_->id( e ) );
      if ( entry != mark_.end() )
        return entry->second;
      else
        return 0;
    }

    /** \brief returns false, if at least one entity is marked for adaption */
    bool preAdapt()
    {
      return mark_.size() > 0;
    }

    /** \brief Triggers the grid adaptation process
      * \return if triangulation has changed
      * \note Refers to the adapt() of the underlying MMesh
      *       because it will definitely change as well
      */
    bool adapt()
    {
      return mMesh_->adapt();
    }

    //! Callback for the grid adaptation process with restrict/prolong
    template< class GridImp, class DataHandle >
    bool adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &handle )
    {
      preAdapt();
      sequence_ += 1;

      for (const auto& ielement : elements( this->leafGridView() ))
        if (ielement.isNew())
        {
          bool initialize = true;
          const auto &connComp = ielement.impl().connectedComponent();

          for ( const auto& old : connComp.children() )
          {
            const Entity& father = old;

            if (father.geometry().volume() > ielement.geometry().volume())
            {
              ielement.impl().bindFather( father );
              handle.prolongLocal( father, ielement, true );
            }
            else
            {
              father.impl().bindFather( ielement );
              handle.restrictLocal( ielement, father, initialize );
            }

            initialize = false;
          }
        }

      postAdapt();
      return true;
    }

    //! Clean up refinement markers
    void postAdapt()
    {
      mark_.clear();
      childrenConnectedComponentMap_.clear();
    }

  public:
    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize(int codim) const {
      return 0;
    }


    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const
    {
      std::size_t size = 0;

      if (codim == 0)
        for (const auto& e : elements(this->leafGridView()))
          if (e.partitionType() == GhostEntity)
            size++;

      if (codim == dimension)
        for (const auto& v : vertices(this->leafGridView()))
          if (v.partitionType() == GhostEntity)
            size++;

      return size;
    }


    /** \brief Size of the overlap on a given level */
    unsigned int overlapSize(int level, int codim) const {
      return overlapSize(codim);
    }


    /** \brief Size of the ghost cell layer on a given level */
    unsigned int ghostSize(int level, int codim) const {
      return ghostSize(codim);
    }


    /** \brief Distributes this grid over the available nodes in a distributed machine
     *
     * \param minlevel The coarsest grid level that gets distributed
     * \param maxlevel does currently get ignored
     */
    void loadBalance() {
      mMesh_->loadBalance();
    };

    template<class T>
    bool loadBalance( const T& t ) { return false; };

    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
      DUNE_THROW(NotImplemented, "MMeshInterfaceGrid::loadBalance()");
    }


    /** \brief Collective communication */
    const Communication< Comm >& comm () const
    {
      return comm_;
    }

    template< class Data, class InterfaceType, class CommunicationDirection >
    void communicate (
      Data &dataHandle,
      InterfaceType interface,
      CommunicationDirection direction,
      int level = 0 ) const
    {
      if (comm().size() <= 1)
        return;

#if HAVE_MPI
      if( (interface == InteriorBorder_All_Interface) || (interface == All_All_Interface) )
      {
        MMeshCommunication<GridImp, MMeshType> communication( getMMesh().partitionHelper() );
        const auto& gv = this->leafGridView();

        switch( direction )
        {
        case ForwardCommunication:
          communication( gv.template begin< 0, Interior_Partition >(), gv.template end< 0, Interior_Partition >(),
                         gv.template begin< 0, Ghost_Partition >(), gv.template end< 0, Ghost_Partition >(),
                         dataHandle, InteriorEntity, GhostEntity,
                         (interface == All_All_Interface) );
          break;

        case BackwardCommunication:
          communication( gv.template begin< 0, Ghost_Partition >(), gv.template end< 0, Ghost_Partition >(),
                         gv.template begin< 0, Interior_Partition >(), gv.template end< 0, Interior_Partition >(),
                         dataHandle, GhostEntity, InteriorEntity,
                         (interface == All_All_Interface) );
          break;
        }
      }
      else
        DUNE_THROW( NotImplemented, "Communication on interface type " << interface << " not implemented." );
#else
        DUNE_THROW( NotImplemented, "MPI not found!" );
#endif //HAVE_MPI
    }

    //! Return if interface segment is part of the interface
    bool isInterface( const MMeshInterfaceEntity<0>& segment ) const
    {
      int count = mMesh_->interfaceSegments().count( getVertexIds_( segment ) );
      assert( count <= 1 );
      return ( count > 0 );
    }

    //! Return if an edge is of the interface
    template< int d = dimension >
    std::enable_if_t< d == 2, bool > isInterface( const MMeshInterfaceEntity<1>& edge ) const
    {
      auto circulator = mMesh_->getHostGrid().incident_facets( edge );
      for ( std::size_t i = 0; i < CGAL::circulator_size( circulator ); ++i, ++circulator )
        if ( isInterface( *circulator ) )
          return true;
      return false;
    }

    //! Return if vertex is part of the interface
    bool isInterface( const VertexHandle& vertex ) const
    {
      return vertex->info().isInterface;
    }

    //! Return if intersection is part of the interface
    bool isInterface( const typename Traits::LeafIntersection& intersection ) const
    {
      return isInterface( intersection.impl().getHostIntersection() );
    }

    //! Return if dune entity is part of the interface
    template< class Entity >
    bool isInterface( const Entity& entity ) const
    {
      return isInterface( entity.impl().hostEntity() );
    }

    //! Return domain marker of entity
    template< class Entity >
    std::size_t domainMarker( const Entity& entity ) const
    {
      assert( isInterface( entity ) );
      return mMesh_->interfaceSegments()[ getVertexIds_( entity.impl().hostEntity() ) ];
    }

    //! Mark a set of children elements as refinement of a connected component
    void markAsRefined( const std::vector< std::vector< std::size_t > >& children, const ConnectedComponent connectedComponent )
    {
      for( const auto& child : children )
        childrenConnectedComponentMap_.insert( { child, connectedComponent } );
    }

    //! Return if an entity has a connected component
    bool hasConnectedComponent( const Entity& entity ) const
    {
      int count = childrenConnectedComponentMap_.count( getVertexIds_( entity.impl().hostEntity() ) );
      assert( count <= 1 );
      return (count == 1 );
    }

    //! Return the connected component for a given entity
    const ConnectedComponent& getConnectedComponent( const Entity& entity ) const
    {
      auto it = childrenConnectedComponentMap_.find( getVertexIds_( entity.impl().hostEntity() ) );
      assert( it != childrenConnectedComponentMap_.end() );
      return it->second;
    }

    //! Return a non-const reference to the connected component for a given entity
    ConnectedComponent& getConnectedComponent( const Entity& entity )
    {
      auto it = childrenConnectedComponentMap_.find( getVertexIds_( entity.impl().hostEntity() ) );
      assert( it != childrenConnectedComponentMap_.end() );
      return it->second;
    }

    // Return if MMeshInterfaceEntity can be mirrored
    bool canBeMirrored(const MMeshInterfaceEntity<0>& hostEntity) const
    {
      using HostGridEntity = typename GridImp::MMeshType::template HostGridEntity<0>;

      if (hostEntity.first->neighbor(hostEntity.second) == HostGridEntity())
        return false;

      if constexpr (dimensionworld == 2)
        return hostEntity.first->dimension() >= 1;

      return true;
    }

    //! Mirror a host entity (2d)
    template <int d = dimensionworld>
    std::enable_if_t< d == 2, const MMeshInterfaceEntity<0> >
    mirrorHostEntity(const MMeshInterfaceEntity<0>& other) const
    {
      return mMesh_->getHostGrid().mirror_edge( other );
    }

    //! Mirror a host entity (3d)
    template <int d = dimensionworld>
    std::enable_if_t< d == 3, const MMeshInterfaceEntity<0> >
    mirrorHostEntity(const MMeshInterfaceEntity<0>& other) const
    {
      return mMesh_->getHostGrid().mirror_facet( other );
    }

    // **********************************************************
    // End of Interface Methods
    // **********************************************************

  private:
    Communication<Comm> comm_;

    static inline auto getVertexIds_( const MMeshInterfaceEntity<0>& entity )
    {
      std::vector<std::size_t> ids( dimensionworld );
      for( int i = 0; i < dimensionworld; ++i )
        ids[ i ] = entity.first->vertex((entity.second+i+1)%(dimensionworld+1))->info().id;
      std::sort(ids.begin(), ids.end());
      return ids;
    }

  public:
    //! compute the grid ids
    void setIds()
    {
      globalIdSet_->update(this);
    }

    //! compute the grid indices
    void setIndices()
    {
      leafIndexSet_->update(this);

      localBoundarySegments_.clear();
      std::size_t count = 0;
      for (const auto& e : elements(this->leafGridView()))
        for (const auto& is : intersections(this->leafGridView(), e))
          if (is.boundary())
          {
            auto iid = globalIdSet_->id( entity( is.impl().getHostIntersection() ) );
            localBoundarySegments_.insert( std::make_pair(iid, count++) );
          }
    }

    //! Return reference to MMesh
    const MMesh& getMMesh() const
    {
      return *mMesh_;
    }

    //! Return reference to underlying CGAL triangulation
    const HostGridType& getHostGrid() const
    {
      return mMesh_->getHostGrid();
    }

    //! Return sequence number
    int sequence() const
    {
      return sequence_;
    }

    const auto& partitionHelper() const
    {
      return mMesh_->partitionHelper();
    }

  private:
    std::unique_ptr<MMeshInterfaceGridLeafIndexSet<const GridImp>> leafIndexSet_;

    std::unique_ptr<MMeshInterfaceGridGlobalIdSet<const GridImp>> globalIdSet_;

    std::unordered_map< std::vector< std::size_t >, ConnectedComponent, HashUIntVector > childrenConnectedComponentMap_;

  protected:
    //! The host grid which contains the actual grid hierarchy structure
    MMesh* mMesh_;

    BoundarySegments boundarySegments_, localBoundarySegments_;

    mutable std::unordered_map< MMeshImpl::MultiId, int > mark_;

    int sequence_ = 0;

  }; // end Class MMeshInterfaceGrid


  // DGFGridInfo for MMeshInterfaceGrid
  // ----------------------------------

  template<class MMesh>
  struct DGFGridInfo< MMeshInterfaceGrid<MMesh> >
  {
    static int refineStepsForHalf ()
    {
      return MMesh::dimension-1;
    }

    static double refineWeight ()
    {
      return -1;
    }
  };

  // Capabilites of MMeshInterfaceGrid
  // ---------------------------------

  /// @cond
  namespace Capabilities
  {
    /** \brief has only one geometry type for all entities
    \ingroup MMeshInterfaceGrid
    */
    template<class MMesh>
    struct hasSingleGeometryType<MMeshInterfaceGrid<MMesh>>
    {
      static const bool v = true;
      static const unsigned int topologyId = Dune::GeometryType::simplex;
    };

    /** \brief has entities for some codimensions
     * \ingroup MMeshInterfaceGrid
     */
    template<class MMesh, int codim>
    struct hasEntity<MMeshInterfaceGrid<MMesh>, codim>
    {
      static const bool v = (codim >= 0 || codim <= MMesh::dimension-1);
    };

    template<class MMesh, int codim>
    struct hasEntityIterator<MMeshInterfaceGrid<MMesh>, codim>
    {
      static const bool v = (codim >= 0 || codim <= MMesh::dimension-1);
    };

    /** \brief has conforming level grids
     * \ingroup MMeshInterfaceGrid
     */
    template<class MMesh>
    struct isLevelwiseConforming<MMeshInterfaceGrid<MMesh>>
    {
      static const bool v = true;
    };

    /** \brief has conforming leaf grids when host grid has
     * \ingroup MMeshInterfaceGrid
     */
    template<class MMesh>
    struct isLeafwiseConforming<MMeshInterfaceGrid<MMesh>>
    {
      static const bool v = true;
    };
  } // end namespace Capabilities
  /// @endcond

} // namespace Dune

#endif // DUNE_MMESH_INTERFACE_GRID_HH
