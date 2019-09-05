// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_MMESH_HH
#define DUNE_MMESH_MMESH_HH

/** \file
 * \brief The MMesh class
 */

#include <config.h>

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <memory>

// dune-common includes
#include <dune/common/deprecated.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/version.hh>
#include <dune/common/to_unique_ptr.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

// CGAL includes
#include <CGAL/utility.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Sphere_3.h>

// The components of the MMesh interface
#include "grid/common.hh"
#include "grid/incidentiterator.hh"
#include "grid/geometry.hh"
#include "grid/entity.hh"
#include "grid/entityseed.hh"
#include "grid/intersectioniterator.hh"
#include "grid/interfaceiterator.hh"
#include "grid/leafiterator.hh"
#include "grid/indexsets.hh"
#include "grid/hierarchiciterator.hh"
#include "grid/mmeshdefaults.hh"
#include "grid/pointfieldvector.hh"
#include "grid/rangegenerators.hh"
#include "interface/traits.hh"
// Further includes below!


namespace Dune
{
  // Forward declarations
  template<int dim, class HostGrid>
  struct MMeshFamily;

  template<class HostGrid, int dim, class GridFamily = MMeshFamily<dim, HostGrid>>
  class MMesh;

  // Type shortcut with default triangulation
  template<int dim>
  using MovingMesh = MMesh< typename MMeshDefaults::Delaunay<dim>::Triangulation, dim >;

  // MMesh family
  template<int dim, class HostGrid>
  struct MMeshFamily
  {

  public:

    using Traits = GridTraits<
        dim,
        dim,
        MMesh<HostGrid, dim>,
        MMeshGeometry,
        MMeshEntity,
        MMeshLeafIterator, // LevelIterator
        MMeshLeafIntersection,
        MMeshLeafIntersection, // LevelIntersection
        MMeshLeafIntersectionIterator,
        MMeshLeafIntersectionIterator, // LevelIntersectionIterator
        MMeshHierarchicIterator,
        MMeshLeafIterator,
        MMeshLeafIndexSet< const MMesh<HostGrid, dim> >, // LevelIndexSet
        MMeshLeafIndexSet< const MMesh<HostGrid, dim> >,
        MMeshGlobalIdSet< const MMesh<HostGrid, dim> >,
        std::size_t, // GlobalIdSet::IdType,
        MMeshLocalIdSet< const MMesh<HostGrid, dim> >,
        std::size_t, // LocalIdSet::IdType,
        CollectiveCommunication< MMesh<HostGrid, dim> >,
        DefaultLevelGridViewTraits,
        DefaultLeafGridViewTraits,
        MMeshEntitySeed
        >;
  };

  /// @cond
  template<class HostGrid, int dim, int codim> class HostGridEntityChooser_ { struct type {}; };

  template<class HG> class HostGridEntityChooser_<HG,2,0> { public: using type = typename HG::Face_handle; };
  template<class HG> class HostGridEntityChooser_<HG,2,1> { public: using type = std::pair<typename HG::Face_handle, int>; };
  template<class HG> class HostGridEntityChooser_<HG,2,2> { public: using type = typename HG::Vertex_handle; };

  template<class HG> class HostGridEntityChooser_<HG,3,0> { public: using type = typename HG::Cell_handle; };
  template<class HG> class HostGridEntityChooser_<HG,3,1> { public: using type = std::pair<typename HG::Cell_handle, int>; };
  template<class HG> class HostGridEntityChooser_<HG,3,2> { public: using type = CGAL::Triple<typename HG::Cell_handle, int, int>; };
  template<class HG> class HostGridEntityChooser_<HG,3,3> { public: using type = typename HG::Vertex_handle; };
  /// @endcond

  //**********************************************************************
  //
  // --MMeshBase class
  //
  //************************************************************************
  /*!
   * \brief Provides the base class of a grid wrapper of a CGAL triangulation
   * \ingroup GridImplementations
   * \ingroup MMesh
   *
   * \tparam HostGrid The CGAL host grid type wrapped by the MMesh
   */
  template <class HostGrid, int dim, class GridFamily>
  class MMesh
   : public GridDefaultImplementation< dim, dim,
                                       /*FieldType=*/typename HostGrid::Point::R::RT,
                                       GridFamily >
  {
  public:
    //! The world dimension
    static constexpr int dimension = dim;

    //! The this type
    using ThisType = MMesh<HostGrid, dim, GridFamily>;

    //! The hostgrid type
    using HostGridType = HostGrid;

    //! The point type
    using Point = typename HostGrid::Point;

    //! The field type
    using FieldType = typename Point::R::RT;

    //! The boundary segment map
    using BoundarySegments = std::unordered_map< std::vector< std::size_t >, std::size_t, HashUIntVector >;

    //! The interface segment set
    using InterfaceSegments = std::unordered_set< std::vector< std::size_t >, HashUIntVector >;

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! The Traits
    using Traits = typename GridFamily::Traits;

    //! The grid implementation
    using GridImp = typename GridFamily::Traits::Grid;

    #if DUNE_VERSION_NEWER(DUNE_GRID, 2, 7)
      typedef ToUniquePtr< GridImp > GridPtrType;
    #else
      typedef GridImp* GridPtrType;
    #endif

    //! The leaf iterator
    using LeafIterator = typename Traits::template Codim<0>::LeafIterator;

    //! The type used for coordinates
    using GlobalCoordinate = Dune::FieldVector<FieldType, dimension>;

    //! The type of the underlying entities
    template<int cd>
    using HostGridEntity = typename HostGridEntityChooser_<HostGridType, dimension, cd>::type;

    //! The type of the underlying element handle
    using ElementHandle = HostGridEntity<0>;

    //! The type of the underlying vertex handle
    using VertexHandle = HostGridEntity<dimension>;

    //! The type of the element output
    using ElementOutput = std::list<HostGridEntity<0>>;

    //! The type of an codim 0 entity
    using Entity = typename Traits::template Codim<0>::Entity;

    //! The type of a vertex
    using Vertex = typename Traits::template Codim<dimension>::Entity;

    //! The type of an interface element
    using InterfaceElement = typename Traits::template Codim<1>::Entity;

    //! The type of the interface grid
    using InterfaceGrid = MMeshInterfaceGrid<GridImp, dimension>;

    //! The type of an interface grid entity
    using InterfaceEntity = Dune::Entity<0, dimension-1, const InterfaceGrid, MMeshInterfaceGridEntity>;

    /** \brief Constructor
     *
     * \param hostgrid          The host grid wrapped by the MMesh
     * \param boundarySegments  The boundary segment index mapper
     * \param interfaceSegments The set of interface segments
     */
    explicit MMesh(HostGrid hostgrid,
                   BoundarySegments boundarySegments,
                   InterfaceSegments interfaceSegments_)
     : hostgrid_(hostgrid),
       boundarySegments_(boundarySegments),
       interfaceSegments_(interfaceSegments_)
    {
      leafIndexSet_ = std::make_unique<MMeshLeafIndexSet<const GridImp>>( This() );
      globalIdSet_ = std::make_unique<MMeshGlobalIdSet<const GridImp>>( This() );
      localIdSet_ = std::make_unique<MMeshLocalIdSet<const GridImp>>( This() );

      setIndices();

      interfaceGrid_ = std::make_shared<InterfaceGrid>( This() );
    }

    explicit MMesh(HostGrid hostgrid)
     : hostgrid_(hostgrid)
    {
      leafIndexSet_ = std::make_unique<MMeshLeafIndexSet<const GridImp>>( This() );
      globalIdSet_ = std::make_unique<MMeshGlobalIdSet<const GridImp>>( This() );
      localIdSet_ = std::make_unique<MMeshLocalIdSet<const GridImp>>( This() );

      setIndices();

      interfaceGrid_ = std::make_shared<InterfaceGrid>( This() );
    }

    virtual ~MMesh() {};

    //! This pointer to derived class
    const GridImp* This() const { return static_cast<const GridImp*>(this); }

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

    //! Returns the number of boundary segments within the macro grid
    size_t numBoundarySegments () const {
      return boundarySegments().size();
    }

    //! Returns the boundary segment to index map
    const BoundarySegments& boundarySegments() const
    {
      return boundarySegments_;
    }

    //! Returns the interface segment set
    const InterfaceSegments& interfaceSegments() const
    {
      return interfaceSegments_;
    }

    //! Returns the interface segment set
    InterfaceSegments interfaceSegments()
    {
      return interfaceSegments_;
    }

    void addInterfaceSegment( const std::vector< std::size_t >& v )
    {
      std::vector<std::size_t> sorted_vertices( v );
      std::sort(sorted_vertices.begin(), sorted_vertices.end());
      interfaceSegments_.insert( sorted_vertices );
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
    const MMeshLeafIndexSet<const GridImp>& levelIndexSet(int level) const
    {
      if (level != 0)
        DUNE_THROW(GridError, "levelIndexSet of nonexisting level " << level << " requested!");
      return *leafIndexSet_;
    }

    /** \brief Access to the LeafIndexSet */
    const MMeshLeafIndexSet<const GridImp>& leafIndexSet() const
    {
      return *leafIndexSet_;
    }

    /** \brief Create Entity from EntitySeed */
    template < class EntitySeed >
    typename Traits::template Codim<EntitySeed::codimension>::Entity
    entity(const EntitySeed& seed) const
    {
      using EntityImp = MMeshEntity<
        EntitySeed::codimension,
        dimension,
        const typename Traits::Grid
        >;

      auto hostEntity = seed.impl().hostEntity();
      assert( hostEntity != decltype(hostEntity)() );
      return EntityImp(This(), hostEntity);
    }

    //! Return the entity corresponding to a vertex handle
    Vertex entity(const HostGridEntity<dimension>& vertexHandle) const
    {
      return entity( typename Traits::template Codim<dimension>::EntitySeed( vertexHandle ) );
    }

    //! Return the entity corresponding to a element handle
    Entity entity(const HostGridEntity<0>& elementHandle) const
    {
      return entity( typename Traits::template Codim<0>::EntitySeed( elementHandle ) );
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This());
    }


    //! One past the end on This() level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This(), true);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This());
    }


    //! One past the end on This() level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This(), true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This());
    }


    //! One past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This(), true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This());
    }


    //! One past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This(), true);
    }

    //! Iterator to first interface entity
    template<int codim>
    MMeshInterfaceIterator<codim, const ThisType> interfaceBegin( bool includeBoundary = false ) const
    {
      using Impl = typename MMeshInterfaceIterator<codim, const ThisType>::Implementation;
      return MMeshInterfaceIterator<codim, const ThisType>( Impl(this, includeBoundary) );
    }

    //! One past the end of the sequence of interface entities
    template<int codim>
    MMeshInterfaceIterator<codim, const ThisType> interfaceEnd( bool includeBoundary = false ) const
    {
      using Impl = typename MMeshInterfaceIterator<codim, const ThisType>::Implementation;
      return MMeshInterfaceIterator<codim, const ThisType>( Impl(this, true, includeBoundary) );
    }

    //! Iterator to first interface entity
    MMeshInterfaceVertexIterator<const ThisType> interfaceVerticesBegin( bool includeBoundary = false ) const
    {
      using Impl = typename MMeshInterfaceVertexIterator<const ThisType>::Implementation;
      return MMeshInterfaceVertexIterator<const ThisType>( Impl(this, includeBoundary) );
    }

    //! One past the end of the sequence of interface entities
    MMeshInterfaceVertexIterator<const ThisType> interfaceVerticesEnd( bool includeBoundary = false ) const
    {
      using Impl = typename MMeshInterfaceVertexIterator<const ThisType>::Implementation;
      return MMeshInterfaceVertexIterator<const ThisType>( Impl(this, true, includeBoundary) );
    }

    /** \brief Global refine
     *
     * Marks all elements for refinement and adapts the grid
     */
    void globalRefine(int steps = 1)
    {
      DUNE_THROW( NotImplemented, "globalRefine for MMesh" );
    }

    /** \brief Mark entity for refinement
     *
     * This only works for entities of codim 0.
     */
    bool mark(int refCount, const typename Traits::template Codim<0>::Entity & e) const
    {
      DUNE_THROW( NotImplemented, "mark for MMesh" );
      return false;
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark
     */
    int getMark(const typename Traits::template Codim<0>::Entity & e) const
    {
      DUNE_THROW( NotImplemented, "getMark for MMesh" );
      return 0;
    }

    //! Returns false, if at least one entity is marked for adaption
    bool preAdapt()
    {
      DUNE_THROW( NotImplemented, "preAdapt for MMesh" );
      return false;
    }

    //! Triggers the grid refinement process
    bool adapt()
    {
      DUNE_THROW( NotImplemented, "adapt for MMesh" );
      return false;
    }

    /** \brief Clean up refinement markers */
    void postAdapt()
    {
      DUNE_THROW( NotImplemented, "postAdapt for MMesh" );
    }

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


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    //! get the host grid
    const HostGrid& getHostGrid() const
    {
      return hostgrid_;
    }

    //! get the host grid
    HostGrid& getHostGrid()
    {
      return hostgrid_;
    }

    //! get the interface grid
    const InterfaceGrid& interfaceGrid() const
    {
      return *interfaceGrid_;
    }

    //! get a pointer to the interface grid
    const std::shared_ptr<InterfaceGrid>& interfaceGridPtr()
    {
      return interfaceGrid_;
    }

    //! Return if element is part of the interface
    bool isInterface( const InterfaceElement& segment ) const
    {
      std::vector<std::size_t> indices;
      for( std::size_t i = 0; i < segment.subEntities(dimension); ++i )
        indices.push_back( segment.impl().template subEntity<dimension>(i).impl().hostEntity()->info().id );

      std::sort(indices.begin(), indices.end());

      int count = interfaceSegments_.count( indices );
      assert( count <= 1 );

      bool isOnBoundary = getHostGrid().is_infinite( (segment.impl().hostEntity().first)->neighbor(segment.impl().hostEntity().second) );
      return ( count == 1 && !isOnBoundary );
    }

    //! Return if entity shares a facet with the interface
    bool isOnInterface( const Entity& entity ) const
    {
      for( std::size_t i = 0; i < entity.subEntities(1); ++i )
        if( isInterface( entity.template subEntity<1>(i) ) )
          return true;

      return false;
    }

    //! Return a codim 1 entity as a interface grid codim 0 entity
    InterfaceEntity asInterfaceEntity( const InterfaceElement& segment ) const
    {
      return InterfaceEntity {{ interfaceGrid_.get(), segment.impl().hostEntity() }};
    }

    //! compute the grid indices and ids
    void setIndices()
    {
      setIds();
      leafIndexSet_->update(This());
    }

  private:
    //! compute the grid indices and ids
    void setIds()
    {
      localIdSet_->update(This());
      globalIdSet_->update(This());
    }

    CollectiveCommunication< GridImp > ccobj;

    std::unique_ptr<MMeshLeafIndexSet<const GridImp>> leafIndexSet_;

    std::unique_ptr<MMeshGlobalIdSet<const GridImp>> globalIdSet_;

    std::unique_ptr<MMeshLocalIdSet<const GridImp>> localIdSet_;

    std::shared_ptr<InterfaceGrid> interfaceGrid_;

    //! The host grid which contains the actual grid hierarchy structure
    HostGrid hostgrid_;

    BoundarySegments boundarySegments_;

    InterfaceSegments interfaceSegments_;

  }; // end Class MMesh

  /// @cond
  namespace Capabilities
  {
    template<class HostGrid, int dim, int codim>
    struct hasEntity<MMesh<HostGrid, dim>, codim>
    {
      static const bool v = (codim >= 0 && codim <= dim);
    };

    template<class HostGrid, int dim, int codim>
    struct hasEntityIterator<MMesh<HostGrid, dim>, codim>
    {
      static const bool v = (codim >= 0 && codim <= dim);
    };

    template<class HostGrid, int dim>
    struct isLevelwiseConforming<MMesh<HostGrid, dim>>
    {
      static const bool v = true;
    };

    template<class HostGrid, int dim>
    struct isLeafwiseConforming<MMesh<HostGrid, dim>>
    {
      static const bool v = true;
    };
  } // end namespace Capabilities
  /// @endcond

} // namespace Dune


template<int codim, int dim, class HostGrid>
std::ostream& operator<< ( std::ostream& stream, const Dune::Entity<codim, dim, const Dune::MMesh<HostGrid, dim>, Dune::MMeshEntity >& entity )
{
  return stream << "MMeshEntity<codim=" << codim << ", dim=" << dim << ", center=(" << entity.geometry().center() << ")>";
}

#include "grid/gridfactory.hh"
#include "grid/structuredgridfactory.hh"
#include "grid/dgfparser.hh"
#include "grid/gmshparser.hh"
#include "interface/grid.hh"

#endif // DUNE_MMESH_MMESH_HH
