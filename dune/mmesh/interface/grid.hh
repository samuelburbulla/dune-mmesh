// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_GRID_HH
#define DUNE_MMESH_INTERFACE_GRID_HH

/** \file
 * \brief The MMesh class
 */

#include <config.h>

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <memory>

#include <dune/common/deprecated.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/version.hh>
#include <dune/common/to_unique_ptr.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

#include <CGAL/utility.h>

// The components of the MMesh interface
#include <dune/mmesh/mmesh.hh>
#include "entity.hh"
#include "entityseed.hh"
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

  // MMeshInterfaceGrid family
  template<int dim, class MMeshImp>
  struct MMeshInterfaceGridFamily
  {

  public:

    typedef GridTraits<
        dim-1, // grid dimension
        dim, // world dimension
        MMeshInterfaceGrid<MMeshImp, dim>,
        MMeshInterfaceGridGeometry,
        MMeshInterfaceGridEntity,
        MMeshInterfaceGridLeafIterator, // LevelIterator
        MMeshInterfaceGridLeafIntersection,
        MMeshInterfaceGridLeafIntersection, // LevelIntersection
        MMeshInterfaceGridLeafIntersectionIterator,
        MMeshInterfaceGridLeafIntersectionIterator, // LevelIntersectionIterator
        MMeshInterfaceGridHierarchicIterator,
        MMeshInterfaceGridLeafIterator,
        MMeshInterfaceGridLeafIndexSet< const MMeshInterfaceGrid<MMeshImp, dim> >, // LevelIndexSet
        MMeshInterfaceGridLeafIndexSet< const MMeshInterfaceGrid<MMeshImp, dim> >,
        MMeshInterfaceGridGlobalIdSet< const MMeshInterfaceGrid<MMeshImp, dim> >,
        std::size_t, // GlobalIdSet::IdType,
        MMeshInterfaceGridGlobalIdSet< const MMeshInterfaceGrid<MMeshImp, dim> >, // LocalIdSet
        std::size_t, // LocalIdSet::IdType,
        CollectiveCommunication< MMeshInterfaceGrid<MMeshImp, dim> >,
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
   * \brief Provides a DUNE grid interface class for the interface of a MIMesh
   * \ingroup GridImplementations
   * \ingroup MMeshInterfaceGrid
   *
   * \tparam MMeshImp The MMesh grid type for which this MMeshInterfaceGrid wraps the interface
   */
  template <class MMeshImp, int dim, class GridFamily>
  class MMeshInterfaceGrid
   : public GridDefaultImplementation< dim-1, dim,
                                       typename MMeshImp::FieldType,
                                       GridFamily >
  {
  public:
    static constexpr int dimension = dim-1;
    static constexpr int dimensionworld = dim;

    using FieldType = typename MMeshImp::FieldType;

    using BoundarySegments = std::map< std::vector< std::size_t >, std::size_t >;

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! the Traits
    typedef typename GridFamily::Traits Traits;

    //! the grid implementation
    using GridImp = typename GridFamily::Traits::Grid;

    #if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
      typedef ToUniquePtr< GridImp > GridPtrType;
    #else
      typedef GridImp* GridPtrType;
    #endif

    //! the underlying hostgrid
    using HostGridType = typename MMeshImp::HostGridType;

    //! the underlying mmesh
    using MMesh = MMeshImp;

    //! the leaf iterator
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;

    //! The type used to store coordinates, inherited from the MMeshImp
    typedef FieldType ctype;

    //! The type used for coordinates
    typedef Dune::FieldVector<ctype, dimensionworld> GlobalCoordinate;

    //! The type of the underlying entities
    template<int cd>
    using MMeshInterfaceEntity = typename HostGridEntityChooser_<HostGridType, dimensionworld, cd+1>::type;

    //! The type of the underlying vertex handle
    using VertexHandle = MMeshInterfaceEntity<dimension>;

    //! The type of the element output
    using ElementOutput = std::list<MMeshInterfaceEntity<0>>;

    //! The type of a codim 0 entity
    using Entity = typename Traits::template Codim<0>::Entity;

    //! The type of a vertex
    using Vertex = typename Traits::template Codim<dimension>::Entity;

    /** \brief Constructor
     *
     * \param hostgrid The host grid wrapped by the MMesh
     */
    explicit MMeshInterfaceGrid(const MMeshImp* mMesh, BoundarySegments boundarySegments = {})
     : mMesh_(mMesh),
       boundarySegments_(boundarySegments)
    {
      leafIndexSet_ = std::make_unique<MMeshInterfaceGridLeafIndexSet<const GridImp>>( this );
      globalIdSet_ = std::make_unique<MMeshInterfaceGridGlobalIdSet<const GridImp>>( this );
      localIdSet_ = std::make_unique<MMeshInterfaceGridGlobalIdSet<const GridImp>>( this );

      setIndices();
    }

    virtual ~MMeshInterfaceGrid() {};

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
      return boundarySegments().size();
    }

    const BoundarySegments& boundarySegments() const
    {
      return boundarySegments_;
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
      return *localIdSet_;
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
      DUNE_THROW( NotImplemented, "globalRefine() for MMeshInterfaceGrid" );
    }


    /** \brief Mark entity for refinement
     *
     * This only works for entities of codim 0.
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e) const
    {
      DUNE_THROW( NotImplemented, "mark() for MMeshInterfaceGrid" );
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark
     */
    int getMark(const typename Traits::template Codim<0>::Entity & e) const
    {
      DUNE_THROW( NotImplemented, "getMark() for MMeshInterfaceGrid" );
    }

    /** \brief returns false, if at least one entity is marked for adaption */
    bool preAdapt()
    {
      DUNE_THROW( NotImplemented, "preAdapt() for MMeshInterfaceGrid" );
    }

    //! Triggers the grid refinement process
    bool adapt()
    {
      DUNE_THROW( NotImplemented, "adapt() for MMeshInterfaceGrid" );
    }

  public:
    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize(int codim) const {
      return 0;
    }


    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const {
      return 0;
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
    void loadBalance()
    {};

    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement)
    {}

    template<class T>
    bool loadBalance( const T& t ) { return false; };


    /** \brief dummy collective communication */
    const CollectiveCommunication< GridImp >& comm () const
    {
      return ccobj;
    }

    template< class Data, class InterfaceType, class CommunicationDirection >
    void communicate (
      Data &data,
      InterfaceType iftype,
      CommunicationDirection dir,
      int level = 0 ) const
    {}

    //! Return if interface segment is part of the interface
    bool isInterface( const MMeshInterfaceEntity<0>& segment ) const
    {
      int count = mMesh_->interfaceSegments().count( getVertexIds_( segment ) );
      assert( count <= 1 );

      bool isOnBoundary = mMesh_->getHostGrid().is_infinite( (segment.first)->neighbor(segment.second) );
      return ( count == 1 && !isOnBoundary );
    }

    //! Return if vertex is part of the interface
    bool isInterface (const Vertex& e) const
    {
      return e.impl().hostEntity()->info().isInterface;
    }

    //! Return if entity is part of the interface
    bool isInterface (const Entity& e) const
    {
      return isInterface( e.impl().hostEntity() );
    }

    // **********************************************************
    // End of Interface Methods
    // **********************************************************

  private:
    CollectiveCommunication< GridImp > ccobj;

    static inline auto getVertexIds_( const MMeshInterfaceEntity<0>& entity )
    {
      std::vector<std::size_t> ids;
      for( int i = 0; i < dimensionworld; ++i )
        ids.push_back( entity.first->vertex((entity.second+i+1)%(dimensionworld+1))->info().id );
      std::sort(ids.begin(), ids.end());
      return ids;
    }

  public:
    //! compute the grid indices and ids
    void setIndices()
    {
      leafIndexSet_->update(this);
    }

    const MMeshImp& getMMesh() const
    {
      return *mMesh_;
    }

    const HostGridType& getHostGrid() const
    {
      return mMesh_->getHostGrid();
    }

  private:
    std::unique_ptr<MMeshInterfaceGridLeafIndexSet<const GridImp>> leafIndexSet_;

    std::unique_ptr<MMeshInterfaceGridGlobalIdSet<const GridImp>> globalIdSet_;

    std::unique_ptr<MMeshInterfaceGridGlobalIdSet<const GridImp>> localIdSet_;

  protected:
    //! The host grid which contains the actual grid hierarchy structure
    const MMeshImp* mMesh_;

    BoundarySegments boundarySegments_;

  }; // end Class MMeshInterfaceGrid


  /// @cond
  namespace Capabilities
  {
    /** \brief has entities for some codimensions
     * \ingroup MMeshInterfaceGrid
     */
    template<class MMeshImp, int dim, int codim>
    struct hasEntity<MMeshInterfaceGrid<MMeshImp, dim>, codim>
    {
      static const bool v = (codim >= 0 || codim <= dim-1);
    };

    template<class MMeshImp, int dim, int codim>
    struct hasEntityIterator<MMeshInterfaceGrid<MMeshImp, dim>, codim>
    {
      static const bool v = (codim >= 0 || codim <= dim-1);
    };

    /** \brief has conforming level grids
     * \ingroup MMeshInterfaceGrid
     */
    template<class MMeshImp, int dim>
    struct isLevelwiseConforming<MMeshInterfaceGrid<MMeshImp, dim>>
    {
      static const bool v = true;
    };

    /** \brief has conforming leaf grids when host grid has
     * \ingroup MMeshInterfaceGrid
     */
    template<class MMeshImp, int dim>
    struct isLeafwiseConforming<MMeshInterfaceGrid<MMeshImp, dim>>
    {
      static const bool v = true;
    };
  } // end namespace Capabilities
  /// @endcond

} // namespace Dune

#endif // DUNE_MMESH_INTERFACE_GRID_HH
