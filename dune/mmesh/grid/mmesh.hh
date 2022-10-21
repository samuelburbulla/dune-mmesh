// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_MMESH_HH
#define DUNE_MMESH_GRID_MMESH_HH

/** \file
 * \brief The MMesh class
 */

#include <string>
#include <unordered_set>
#include <unordered_map>
#include <memory>

// dune-common includes
#include <dune/common/deprecated.hh>
#include <dune/common/parallel/communication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/version.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/adaptcallback.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

// The components of the MMesh interface
#include "common.hh"
#include "connectedcomponent.hh"
#include "cutsettriangulation.hh"
#include "incidentiterator.hh"
#include "geometry.hh"
#include "entity.hh"
#include "entityseed.hh"
#include "intersectioniterator.hh"
#include "interfaceiterator.hh"
#include "leafiterator.hh"
#include "indexsets.hh"
#include "hierarchiciterator.hh"
#include "pointfieldvector.hh"
#include "rangegenerators.hh"
#include "../remeshing/distance.hh"
#include "../remeshing/longestedgerefinement.hh"
#include "../remeshing/ratioindicator.hh"
#include "../interface/traits.hh"
#include "../misc/boundaryidprovider.hh"
#include "../misc/twistutility.hh"
// Further includes below!

#if HAVE_MPI
  #include <dune/common/parallel/mpicommunication.hh>
  #include "../misc/communication.hh"
  #include "../misc/partitionhelper.hh"
  using Comm = MPI_Comm;
#else
  using Comm = Dune::No_Comm;
#endif

namespace Dune
{

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
        MMeshImpl::MultiId,
        MMeshGlobalIdSet< const MMesh<HostGrid, dim> >, // LocalIdSet
        MMeshImpl::MultiId,
        Communication<Comm>,
        DefaultLevelGridViewTraits,
        DefaultLeafGridViewTraits,
        MMeshEntitySeed
        >;
  };

  /// @cond
  /*!
   * \brief Determine the correct entity type in the CGAL host grid
   * \ingroup MMesh
   */
  template<class HostGrid, int dim, int codim> class HostGridEntityChooser_ { struct type {}; };

  template<class HG> class HostGridEntityChooser_<HG,2,0> { public: using type = typename HG::Face_handle; };
  template<class HG> class HostGridEntityChooser_<HG,2,1> { public: using type = std::pair<typename HG::Face_handle, int>; };
  template<class HG> class HostGridEntityChooser_<HG,2,2> { public: using type = typename HG::Vertex_handle; };

  template<class HG> class HostGridEntityChooser_<HG,3,0> { public: using type = typename HG::Cell_handle; };
  template<class HG> class HostGridEntityChooser_<HG,3,1> { public: using type = std::pair<typename HG::Cell_handle, int>; };
  template<class HG> class HostGridEntityChooser_<HG,3,2> { public: using type = CGAL::Triple<typename HG::Cell_handle, int, int>; };
  template<class HG> class HostGridEntityChooser_<HG,3,3> { public: using type = typename HG::Vertex_handle; };

  //! The refinement insertion point struct
  template<class Point, class Edge, class IdType, class VertexHandle, class InterfaceGridConnectedComponent>
  struct RefinementInsertionPointStruct
  {
    Point point;
    std::size_t insertionLevel = 0;
    Edge edge;
    IdType edgeId;
    VertexHandle v0, v1;
    bool isInterface = false;
    InterfaceGridConnectedComponent connectedcomponent;
  };
  /// @endcond

  //**********************************************************************
  //
  // --MMesh class
  //
  //************************************************************************
  /*!
   * \brief The MMesh class templatized by the CGAL host grid type and the dimension.
   * \ingroup GridImplementations
   * \ingroup MMesh
   */
  template < class HostGrid, int dim >
  class MMesh
#ifndef DOXYGEN_SHOULD_SKIP_THIS
   : public GridDefaultImplementation< dim, dim,
                                       double, ///*FieldType=*/typename HostGrid::Point::R::RT,
                                       MMeshFamily<dim, HostGrid> >
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
  {
  public:
    //! The world dimension
    static constexpr int dimension = dim;

    //! The hostgrid type
    using HostGridType = HostGrid;

    //! The grid family type
    using GridFamily = MMeshFamily<dim, HostGrid>;

    //! The point type
    using Point = typename HostGrid::Point;

    //! The field type
    using FieldType = typename Point::R::RT;

    //! The type of an id
    using IdType = MMeshImpl::MultiId;

    //! The boundary segment map
    using BoundarySegments = std::unordered_map< IdType, std::size_t >;

    //! The boundary id map
    using BoundaryIds = std::unordered_map< std::size_t, std::size_t >;

    //! The interface segment set
    using InterfaceSegments = std::unordered_map< IdType, std::size_t >;

    //**********************************************************
    // The Interface Methods
    //**********************************************************

    //! The Traits
    using Traits = typename GridFamily::Traits;

    //! The grid implementation
    using GridImp = typename GridFamily::Traits::Grid;

    //! The unique pointer to the grid
    using GridPtrType = std::unique_ptr< GridImp >;

    //! The leaf iterator
    using LeafIterator = typename Traits::template Codim<0>::LeafIterator;

    //! The type used for coordinates
    using GlobalCoordinate = Dune::FieldVector<FieldType, dimension>;

    //! The type of the underlying entities
    template<int cd>
    using HostGridEntity = typename HostGridEntityChooser_<HostGridType, dimension, cd>::type;

    //! The type of the underlying element handle
    using ElementHandle = HostGridEntity<0>;

    //! The type of the underlying edge handle
    using FacetHandle = HostGridEntity<1>;

    //! The type of the underlying edge handle
    using EdgeHandle = HostGridEntity<dimension-1>;

    //! The type of the underlying vertex handle
    using VertexHandle = HostGridEntity<dimension>;

    //! The type of the element output
    using ElementOutput = std::list<HostGridEntity<0>>;

    //! The type of the boundary edges output
    using BoundaryEdgesOutput = std::list<FacetHandle>;

    //! The type of a codim 0 entity
    using Entity = typename Traits::template Codim<0>::Entity;

    //! The type of a codim 1 entity ('facet')
    using Facet = typename Traits::template Codim<1>::Entity;

    //! The type of a codim dim-1 entity ('edge')
    using Edge = typename Traits::template Codim<dimension-1>::Entity;

    //! The type of a codim dim entity ('vertex')
    using Vertex = typename Traits::template Codim<dimension>::Entity;

    //! The type of a caching entity
    using CachingEntity = MMeshCachingEntity< 0, dimension, const GridImp >;

    //! The type used to store connected components of entities
    using ConnectedComponent = MMeshConnectedComponent< const GridImp >;

    //! The type of an interface element
    using InterfaceElement = typename Traits::template Codim<1>::Entity;

    //! The type of an intersection
    using Intersection = typename Traits::LeafIntersection;

    //! The type of the interface grid
    using InterfaceGrid = MMeshInterfaceGrid<GridImp>;

    //! The type of an interface grid entity
    using InterfaceEntity = typename InterfaceGrid::Traits::template Codim<0>::Entity;

    //! The type of a connected component of interface grid entities
    using InterfaceGridConnectedComponent = MMeshInterfaceConnectedComponent<const InterfaceGrid>;

    //! The type of a refinement insertion point
    using RefinementInsertionPoint = RefinementInsertionPointStruct<Point, Edge, IdType, VertexHandle, InterfaceGridConnectedComponent>;

    //! The type of the employed remeshing indicator
    using RemeshingIndicator = RatioIndicator<GridImp>;

    //! The type of the employed refinement strategy
    using RefinementStrategy = LongestEdgeRefinement<GridImp>;

    //! The type of the employed refinement strategy for the interface grid
    using InterfaceRefinementStrategy = LongestEdgeRefinement<InterfaceGrid>;

    //! Constructor that takes a CGAL triangulation
    explicit MMesh(HostGrid hostgrid)
     : MMesh(hostgrid, {}, {}, {}, {}) {}

    //! Constructor that takes additional about interface and boundary.
    explicit MMesh(HostGrid hostgrid,
                   BoundarySegments boundarySegments,
                   BoundarySegments interfaceBoundarySegments,
                   BoundaryIds boundaryIds,
                   InterfaceSegments interfaceSegments)
     : hostgrid_(hostgrid),
       boundarySegments_(boundarySegments),
       boundaryIds_(boundaryIds),
       interfaceSegments_(interfaceSegments),
#ifdef HAVE_MPI
       comm_( MPIHelper::getCommunicator() ),
#endif
       partitionHelper_(*this)
    {
      leafIndexSet_ = std::make_unique<MMeshLeafIndexSet<const GridImp>>( This() );
      globalIdSet_ = std::make_unique<MMeshGlobalIdSet<const GridImp>>( This() );
      globalIdSet_->update(This());

      interfaceGrid_ = std::make_shared<InterfaceGrid>( This(), interfaceBoundarySegments );
      loadBalance();
      indicator_.init(*this);
    }

    //! This pointer to derived class
    const GridImp* This() const { return static_cast<const GridImp*>(this); }
    GridImp* This() { return static_cast<GridImp*>(this); }

  public:
    //! update the grid indices and ids
    void update()
    {
      setIds();
      setIndices();
      interfaceGrid_->setIds();
      interfaceGrid_->setIndices();
    }

  private:
    //! compute the grid ids
    void setIds()
    {
      globalIdSet_->update(This());
    }

    //! compute the grid indices
    void setIndices()
    {
      leafIndexSet_->update(This());

      if (comm().size() == 1)
        localBoundarySegments_ = boundarySegments_;
      else
      {
        localBoundarySegments_.clear();
        std::size_t count = 0;
        for (const auto& e : elements(this->leafGridView()))
          for (const auto& is : intersections(this->leafGridView(), e))
            if (is.boundary())
            {
              IdType iid = globalIdSet_->template id<1>( entity( is.impl().getHostIntersection() ) );
              localBoundarySegments_.insert( std::make_pair(iid, count++) );
            }
      }
    }

  public:
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

    //! returns the number of boundary segments within the macro grid
    size_t numBoundarySegments () const {
      return localBoundarySegments_.size();
    }

    //! returns the boundary segment to index map
    const BoundarySegments& boundarySegments() const
    {
      return localBoundarySegments_;
    }

    //! returns the boundary segment index to boundary id map
    const BoundaryIds& boundaryIds() const
    {
      return boundaryIds_;
    }

    //! returns the interface segment set
    const InterfaceSegments& interfaceSegments() const
    {
      return interfaceSegments_;
    }

    //! returns the interface segment set
    InterfaceSegments& interfaceSegments()
    {
      return interfaceSegments_;
    }

    //! Add an intersection to the interface
    void addInterface( const Intersection& intersection, const std::size_t marker = 1 )
    {
      if (isInterface( intersection ))
        return;

      const auto& facet = entity( intersection.impl().getHostIntersection() );
      std::vector<std::size_t> ids;
      for( std::size_t i = 0; i < facet.subEntities(dim); ++i )
      {
        const auto& vertex = facet.impl().template subEntity<dim>(i);
        vertex.impl().hostEntity()->info().isInterface = true;
        ids.push_back( globalIdSet().id( vertex ).vt()[0] );
      }
      std::sort(ids.begin(), ids.end());
      interfaceSegments_.insert( std::make_pair(ids, marker) );

      // Add interface element to connected component in order to mark element as new
      const auto& ientity = asInterfaceEntity( intersection );
      InterfaceGridConnectedComponent connectedComponent ( ientity );
      interfaceGrid_->markAsRefined( {ids}, connectedComponent );

      interfaceGrid_->setIds();
      interfaceGrid_->setIndices();

      if (comm().size() > 0)
        partitionHelper_.computeInterfacePartitions();
    }

    //! Add wrapped intersection to the interface
    template <class I>
    void addInterface( const I& intersection, const std::size_t marker = 1 )
    {
      addInterface(intersection.impl().hostIntersection(), marker);
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

    InterfaceElement entity(const FacetHandle& facetHandle) const
    {
      return entity( typename Traits::template Codim<1>::EntitySeed( facetHandle ) );
    }

    template< int d = dim >
    std::enable_if_t< d == 3, Edge >
    entity(const HostGridEntity<2>& edgeHandle) const
    {
      return entity( typename Traits::template Codim<2>::EntitySeed( edgeHandle ) );
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This());
    }


    //! one past the end on this level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This(), true);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This());
    }


    //! one past the end on this level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lend (int level) const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This(), true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafbegin() const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This());
    }


    //! one past the end of the sequence of leaf entities
    template<int codim>
    typename Traits::template Codim<codim>::LeafIterator leafend() const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This(), true);
    }


    //! Iterator to first leaf entity of given codim
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafbegin() const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This());
    }


    //! one past the end of the sequence of leaf entities
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LeafIterator leafend() const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This(), true);
    }

    typename Traits::LevelIntersectionIterator ilevelbegin( const typename Traits::template Codim< 0 >::Entity& entity ) const
    {
      return entity.impl().ilevelbegin();
    }

    typename Traits::LevelIntersectionIterator ilevelend( const typename Traits::template Codim< 0 >::Entity& entity ) const
    {
      return entity.impl().ilevelend();
    }

    typename Traits::LeafIntersectionIterator ileafbegin( const typename Traits::template Codim< 0 >::Entity& entity ) const
    {
      return entity.impl().ileafbegin();
    }

    typename Traits::LeafIntersectionIterator ileafend( const typename Traits::template Codim< 0 >::Entity& entity ) const
    {
      return entity.impl().ileafend();
    }

    //! iterator to first interface entity
    template<int codim>
    MMeshInterfaceIterator<codim, const GridImp> interfaceBegin( bool includeBoundary = false ) const
    {
      using Impl = typename MMeshInterfaceIterator<codim, const GridImp>::Implementation;
      return MMeshInterfaceIterator<codim, const GridImp>( Impl(this, false) );
    }

    //! one past the end of the sequence of interface entities
    template<int codim>
    MMeshInterfaceIterator<codim, const GridImp> interfaceEnd( bool includeBoundary = false ) const
    {
      using Impl = typename MMeshInterfaceIterator<codim, const GridImp>::Implementation;
      return MMeshInterfaceIterator<codim, const GridImp>( Impl(this, true, includeBoundary) );
    }

    //! iterator to first interface entity
    MMeshInterfaceVertexIterator<const GridImp> interfaceVerticesBegin( bool includeBoundary = false ) const
    {
      using Impl = typename MMeshInterfaceVertexIterator< const GridImp>::Implementation;
      return MMeshInterfaceVertexIterator<const GridImp>( Impl(this, includeBoundary) );
  }

    //! one past the end of the sequence of interface entities
    MMeshInterfaceVertexIterator<const GridImp> interfaceVerticesEnd( bool includeBoundary = false ) const
    {
      using Impl = typename MMeshInterfaceVertexIterator<const GridImp>::Implementation;
      return MMeshInterfaceVertexIterator<const GridImp>( Impl(this, true, includeBoundary) );
    }

    //! Return if vertex is part of the interface
    bool isInterface( const Vertex& vertex ) const
    {
      return vertex.impl().hostEntity()->info().isInterface;
    }

    //! Return if intersection is part of the interface
    bool isInterface( const Intersection& intersection ) const
    {
      return isInterface( entity( intersection.impl().getHostIntersection() ) );
    }

    //! Return if intersection is part of the interface
    template< class OtherIntersection >
    bool isInterface( const OtherIntersection& intersection ) const
    {
      return isInterface( intersection.impl().hostIntersection() );
    }

    //! Return if element is part of the interface
    bool isInterface( const InterfaceElement& segment ) const
    {
      static std::vector<std::size_t> ids( segment.subEntities(dimension) );
      for( std::size_t i = 0; i < segment.subEntities(dimension); ++i )
      {
        const auto& vertex = segment.impl().template subEntity<dimension>(i);
        if ( !vertex.impl().isInterface() )
          return false;

        ids[i] = this->globalIdSet().id( vertex ).vt()[0];
      }

      std::sort(ids.begin(), ids.end());

      int count = interfaceSegments_.count( ids );
      assert( count <= 1 );
      return ( count > 0 );
    }

    //! Return if edge in 3d is part of an interface segment
    template< int d = dim >
    std::enable_if_t< d == 3, bool >
    isInterface( const Edge& edge ) const
    {
      std::vector<std::size_t> ids;
      for( std::size_t i = 0; i < edge.subEntities(dimension); ++i )
      {
        const auto& vertex = edge.impl().template subEntity<dimension>(i);

        if( !isInterface( vertex ) )
          return false;

        ids.push_back( this->globalIdSet().id( vertex ).vt()[0] );
      }

      std::sort(ids.begin(), ids.end());

      for ( const auto& iseg : interfaceSegments_ )
      {
        const auto& seg = iseg.first.vt();
        if ( (seg[0] == ids[0] && seg[1] == ids[1])
          || (seg[1] == ids[0] && seg[2] == ids[1])
          || (seg[0] == ids[0] && seg[2] == ids[1]) )
          return true;
      }

      return false;
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

    //! Return an intersection as a interface grid codim 0 entity
    InterfaceEntity asInterfaceEntity( const Intersection& intersection ) const
    {
      return InterfaceEntity {{ interfaceGrid_.get(), intersection.impl().getHostIntersection() }};
    }

    //! Return an intersection as a interface grid codim 0 entity
    template< class OtherIntersection >
    InterfaceEntity asInterfaceEntity( const OtherIntersection& intersection ) const
    {
      return asInterfaceEntity( intersection.impl().hostIntersection() );
    }

    //! Return an interface entity as intersection of a MMesh entity
    Intersection asIntersection( const InterfaceEntity& interfaceEntity ) const
    {
      const auto& host = interfaceEntity.impl().hostEntity();
      return asIntersection(host);
    }

    //! Return a facet as intersection
    Intersection asIntersection( const Facet& facet ) const
    {
      const auto& host = facet.impl().hostEntity();
      return asIntersection(host);
    }

    //! Return a host facet as intersection
    Intersection asIntersection( const FacetHandle& host ) const
    {
      FacetHandle hostFacet = host;

      // make sure to get the intersection seen from inside at domain boundary
      if ( hostgrid_.is_infinite(hostFacet.first) )
        hostFacet = interfaceGrid().mirrorHostEntity( hostFacet );

      return MMeshLeafIntersection<const GridImp> ( This(), hostFacet.first, hostFacet.second );
    }

    //! Locate an entity by coordinate using CGAL's locate
    const Entity locate ( const GlobalCoordinate& p, const Entity& element = {} ) const
    {
      return entity( hostgrid_.inexact_locate( makePoint( p ), element.impl().hostEntity() ) );
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
      e.impl().mark( refCount );
      if(refCount > 0) ++refineMarked_;
      if(refCount < 0) ++coarsenMarked_;
      return true;
    }

    /** \brief Mark elements for adaption using the default remeshing indicator
     * \return if elements have been marked.
     */
    bool markElements()
    {
      bool change = false;
      indicator_.update();

      for (const auto& element : elements( this->leafGridView() ))
      {
        mark( indicator_(element), element );
        change |= indicator_(element) != 0;
      }

      for (const auto& ielement : elements( interfaceGrid_->leafGridView() ))
      {
        interfaceGrid_->mark( indicator_(ielement), ielement );
        change |= indicator_(ielement) != 0;
      }

      return change;
    }

    /** \brief Return refinement mark for entity
     *
     * \return refinement mark
     */
    int getMark(const typename Traits::template Codim<0>::Entity & e) const
    {
      return e.impl().getMark();
    }

    /** \brief returns false, if at least one entity is marked for adaption */
    bool preAdapt()
    {
      return (interfaceGrid_->preAdapt()) || (coarsenMarked_ > 0) || (remove_.size() > 0);
    }

    /** \brief Refine edge manually */
    void refineEdge(const Entity& entity, const std::size_t edgeIndex, const double where = 0.5)
    {
      RefinementInsertionPoint ip;
      ip.edge = entity.template subEntity<dim-1>(edgeIndex);
      ip.edgeId = globalIdSet().id( ip.edge );
      ip.point = makePoint( where * ip.edge.geometry().corner(0) + (1 - where) * ip.edge.geometry().corner(1) );
      ip.v0 = ip.edge.impl().template subEntity<dim>(0).impl().hostEntity();
      ip.v1 = ip.edge.impl().template subEntity<dim>(1).impl().hostEntity();
      ip.insertionLevel = ip.edge.impl().insertionLevel() + 1;

      if ( isInterface( ip.edge ) )
      {
        ip.isInterface = true;
        if constexpr (dim != 3)
        {
          InterfaceEntity component {{ interfaceGrid_.get(), ip.edge.impl().hostEntity() }};
          ip.connectedcomponent = InterfaceGridConnectedComponent( component );
        }
      }

      if ( inserted_.insert( ip.edgeId ).second )
      {
        insert_.push_back( ip );
        if (verbose_)
          std::cout << "Insert vertex manually: " << ip.point << std::endl;
      }
    }

    /** \brief Insert vertex in cell manually */
    void insertVertexInCell(const GlobalCoordinate& position)
    {
      RefinementInsertionPoint ip;
      ip.point = makePoint( position );

      insert_.push_back( ip );
      if (verbose_)
        std::cout << "Insert vertex in cell manually: " << ip.point << std::endl;
    }

    /** \brief Remove interface vertex manually */
    void removeVertex(const typename InterfaceGrid::Traits::template Codim<dim-1>::Entity& interfaceVertex)
    {
      const Vertex& vertex = entity( interfaceVertex.impl().hostEntity() );
      removeVertex( vertex );
    }

    /** \brief Remove vertex manually */
    void removeVertex(const Vertex& vertex)
    {
      if ( removed_.insert( globalIdSet().id( vertex ) ).second )
      {
        remove_.push_back( vertex.impl().hostEntity() );
        if (verbose_)
          std::cout << "Remove vertex manually: " << vertex.geometry().center() << std::endl;
      }
    }

    /** \brief Triggers the grid adaptation process
      * \return if triangulation has changed
      */
    bool adapt()
    {
      // Obtain the adaption points
      for ( const auto& element : elements( this->leafGridView() ) )
      {
        int mark = element.impl().getMark();

        // refine
        if (mark == 1)
        {
          std::pair<Edge, GlobalCoordinate> pair = RefinementStrategy::refinement(element);

          RefinementInsertionPoint ip;
          ip.edge = pair.first;
          ip.edgeId = globalIdSet().id( ip.edge );
          ip.point = makePoint( pair.second );
          ip.v0 = ip.edge.impl().template subEntity<dim>(0).impl().hostEntity();
          ip.v1 = ip.edge.impl().template subEntity<dim>(1).impl().hostEntity();
          ip.insertionLevel = ip.edge.impl().insertionLevel() + 1;

          if ( !isInterface( ip.edge ) )
            if ( inserted_.insert( ip.edgeId ).second )
            {
              insert_.push_back( ip );
              if (verbose_)
                std::cout << "Insert vertex because of marked cell: " << ip.point << std::endl;
            }
        }

        // coarsen
        else if (mark == -1)
        {
          auto vertex = RefinementStrategy::coarsening(element);

          if ( vertex != decltype(vertex)() )
            if ( removed_.insert( globalIdSet().id( vertex ) ).second )
            {
              remove_.push_back( vertex.impl().hostEntity() );
              if (verbose_)
                std::cout << "Remove vertex because of marked cell: " << vertex.geometry().center() << std::endl;
            }
        }
      }

      // Obtain the adaption points for the interface
      for ( const auto& element : elements( this->interfaceGrid().leafGridView() ) )
      {
        int mark = this->interfaceGrid().getMark( element );

        // refine
        if (mark == 1)
        {
          auto pair = InterfaceRefinementStrategy::refinement(element);

          RefinementInsertionPoint ip;
          ip.edge = entity( pair.first.impl().hostEntity() );
          ip.edgeId = globalIdSet().id( ip.edge );
          ip.point = makePoint( pair.second );
          ip.v0 = ip.edge.impl().template subEntity<dim>(0).impl().hostEntity();
          ip.v1 = ip.edge.impl().template subEntity<dim>(1).impl().hostEntity();
          ip.insertionLevel = ip.edge.impl().insertionLevel() + 1;
          ip.isInterface = true;
          ip.connectedcomponent = InterfaceGridConnectedComponent( element );

          if ( inserted_.insert( ip.edgeId ).second )
          {
            insert_.push_back( ip );
            if (verbose_)
              std::cout << "Insert interface vertex because of marked cell: " << ip.point << std::endl;
          }
        }

        // coarsen
        else if (mark == -1)
        {
          auto ivertex = InterfaceRefinementStrategy::coarsening(element);
          if ( ivertex != decltype(ivertex)() )
          {
            const Vertex& vertex = entity( ivertex.impl().hostEntity() );

            if ( removed_.insert( globalIdSet().id( vertex ) ).second )
            {
              remove_.push_back( vertex.impl().hostEntity() );
              if (verbose_)
                std::cout << "Remove interface vertex because of marked cell: " << vertex.geometry().center() << std::endl;
            }
          }
        }
      }

      if (verbose_)
        std::cout << "- insert " << insert_.size() << "\t remove " << remove_.size() << std::endl;

      return adapt_(true);
    }

    /** \brief Mark elements such that after movement of interface vertices no cell degenerates
     * \param shifts Vector that maps interface vertex index to global coordinate
     * \return If elements have been marked.
     */
    bool ensureInterfaceMovement( std::vector<GlobalCoordinate> shifts )
    {
      assert( shifts.size() == this->interfaceGrid().leafIndexSet().size(dimension-1) );
      bool change = false;

      // check if grid is still valid
      for ( const auto& element : elements( this->leafGridView() ) )
        if ( signedVolume_( element ) <= 0.0 )
          DUNE_THROW( GridError, "A cell has a negative volume! Maybe the interface has been moved too far?" );

      // temporarily move vertices
      moveInterface( shifts );

      for ( const auto& element : elements( this->leafGridView() ) )
        if ( signedVolume_( element ) <= 0.0 )
        {
          // disable removing of vertices in 3d
          if constexpr ( dim == 3 )
          {
            DUNE_THROW( GridError, "Interface could not be moved, because the removal of vertices is not supported in 3D!" );
          }
          else
          {
            bool allInterface = true;
            for( std::size_t j = 0; j < element.subEntities(dimension); ++j )
            {
              const auto& v = element.template subEntity<dimension>(j);

              // if vertex is part of interface, we have to do more
              if ( !v.impl().isInterface() )
              {
                if ( RefinementStrategy::boundaryFlag( v ) == 1 )
                  continue;

                bool inserted = removed_.insert( globalIdSet().id( v ) ).second;
                if ( inserted )
                {
                  remove_.push_back( v.impl().hostEntity() );
                  change = true;
                  if (verbose_)
                    std::cout << "Remove vertex because of negative volume: " << v.geometry().center() << std::endl;
                }
                allInterface = false;
              }
            }

            // if all vertices are part of the interface, we can try to compute an intersection
            if ( allInterface )
            {
              Edge interfaceEdge;
              for( std::size_t j = 0; j < element.subEntities(1); ++j )
              {
                const auto& edge = element.template subEntity<1>(j);
                if( isInterface( edge ) )
                {
                  if (interfaceEdge != Edge()) // with more than one interface edge we don't know what to do
                    DUNE_THROW( GridError, "Interface could not be moved because a constrained cell with more than one interface edge would have negative volume!" );

                  interfaceEdge = edge;
                }
              }

               // with no interface edge we can try to coarsen the interface
              if (interfaceEdge == Edge())
              {
                bool allNonRemovable = true;
                for( std::size_t j = 0; j < element.subEntities(dimension); ++j )
                {
                  const auto v = element.template subEntity<dimension>(j);
                  const auto iv = interfaceGrid().entity( v.impl().hostEntity() );
                  if ( RefinementStrategy::isRemoveable( iv ) )
                  {
                    allNonRemovable = false;
                    bool inserted = removed_.insert( globalIdSet().id( v ) ).second;
                    if ( inserted )
                    {
                      remove_.push_back( v.impl().hostEntity() );
                      change = true;
                      if (verbose_)
                        std::cout << "Remove interface vertex because of negative volume: " << v.geometry().center() << std::endl;
                    }
                  }
                }
                if ( allNonRemovable )
                  DUNE_THROW( GridError, "Interface could not be moved because a constrained cell without any interface edge would have negative volume!" );
              }
              else
              {
                // compute intersection
                Vertex thirdVertex;
                for( std::size_t j = 0; j < element.subEntities(dimension); ++j )
                {
                  const auto& v = element.template subEntity<dimension>(j);
                  if ( v != interfaceEdge.impl().template subEntity<dimension>(0)
                    && v != interfaceEdge.impl().template subEntity<dimension>(1) )
                    thirdVertex = v;
                }

                Vertex adjVertex;
                InterfaceEntity crossingEdge;

                const auto& iThirdVertex = interfaceGrid().entity( thirdVertex.impl().hostEntity() );
                for ( const auto& e : incidentInterfaceElements( iThirdVertex ) )
                {
                  if ( crossingEdge != InterfaceEntity() )
                    DUNE_THROW( GridError, "Interface could not be moved because two interfaces cross!" );

                  crossingEdge = e;

                  if ( e.impl().template subEntity<1>(0) == iThirdVertex )
                    adjVertex = entity( e.impl().template subEntity<1>(1).impl().hostEntity() );
                  else
                    adjVertex = entity( e.impl().template subEntity<1>(0).impl().hostEntity() );
                }

                // compute intersection
                const auto iegeo = interfaceEdge.geometry();
                const auto& c0 = iegeo.corner(0);
                const auto& c1 = iegeo.corner(1);
                const auto cegeo = crossingEdge.geometry();
                GlobalCoordinate x = Dune::PolygonCutting<double, GlobalCoordinate>::lineIntersectionPoint( c0, c1, cegeo.corner(0), cegeo.corner(1) );

                // check if intersection point is in interfaceEdge
                if ( (c0 - x)*(c1 - x) > 0. )
                  continue;

                RefinementInsertionPoint ip;
                ip.edge = interfaceEdge;
                ip.edgeId = globalIdSet().id( interfaceEdge );
                ip.point = makePoint( x );
                ip.v0 = thirdVertex.impl().hostEntity();
                ip.insertionLevel = thirdVertex.impl().insertionLevel() + 1;
                ip.isInterface = true;
                ip.connectedcomponent = InterfaceGridConnectedComponent(
                  (interfaceEdge.geometry().volume() > crossingEdge.geometry().volume()) // take the edge with the larger volume
                  ? interfaceGrid().entity( interfaceEdge.impl().hostEntity() ) : crossingEdge
                );

                if (verbose_)
                  std::cout << "Insert interface intersection point: " << ip.point << std::endl;

                insert_.push_back( ip );
                change = true;

                // move third vertex back to old position
                const auto& idx = interfaceGrid().leafIndexSet().index(iThirdVertex);
                thirdVertex.impl().hostEntity()->point() = makePoint( thirdVertex.geometry().center() - shifts[idx] );
                shifts[idx] = GlobalCoordinate( 0.0 );
              }
            }
          }
        }

      // move vertices back
      for ( GlobalCoordinate& s : shifts )
        s *= -1.0;
      moveInterface( shifts );

      return change;
    }

    /** \brief Mark elements such that after movement of vertices no cell degenerates
     * \param shifts Vector that maps vertex index to GlobalCoordinate
     * \return If elements have been marked.
     */
    bool ensureVertexMovement( std::vector<GlobalCoordinate> shifts )
    {
      assert( shifts.size() == this->leafIndexSet().size(dimension) );
      bool change = false;

      // temporarily move vertices
      moveVertices( shifts );

      for ( const auto& element : elements( this->leafGridView() ) )
        if ( signedVolume_( element ) <= 0.0 )
        {
          // disable removing of vertices in 3d
          if constexpr ( dim == 3 )
          {
            DUNE_THROW( GridError, "Vertices could not be moved, because the removal of vertices is not supported in 3D!" );
          }
          else
          {
            for( std::size_t j = 0; j < element.subEntities(dimension); ++j )
            {
              const auto& v = element.template subEntity<dimension>(j);

              // if vertex is part of interface, we have to do more
              if ( !v.impl().isInterface() )
              {
                if ( RefinementStrategy::boundaryFlag( v ) == 1 )
                  continue;

                bool inserted = removed_.insert( globalIdSet().id( v ) ).second;
                if ( inserted )
                {
                  remove_.push_back( v.impl().hostEntity() );
                  change = true;
                  if (verbose_)
                    std::cout << "Remove vertex because of negative volume: " << v.geometry().center() << std::endl;
                }
              }
            }
          }
        }

      // move vertices back
      for ( GlobalCoordinate& s : shifts )
        s *= -1.0;
      moveVertices( shifts );

      return change;
    }

    //! Callback for the grid adaptation process with restrict/prolong
    template< class GridImp, class DataHandle >
    bool adapt ( AdaptDataHandleInterface< GridImp, DataHandle > &handle )
    {
      preAdapt();
      adapt();
      sequence_ += 1;

      for (const auto& element : elements( this->leafGridView() ))
        if (element.isNew())
        {
          bool initialize = true;

          for ( const auto& old : element.impl().connectedComponent().children() )
          {
            const Entity& father = old;

            MMeshImpl::CutSetTriangulation<Entity> cutSetTriangulation( old, element );
            for ( auto& middle : cutSetTriangulation.triangles() )
            {
              middle.impl().bindFather( father );
              handle.prolongLocal( father, middle, true );
              middle.impl().bindFather( element );
              handle.restrictLocal( element, middle, initialize );

              initialize = false;
            }
          }
        }

      postAdapt();
      return true;
    }

  private:
    template<int d = dim>
    std::enable_if_t< d == 2, FieldType >
    signedVolume_( const Entity& element ) const
    {
      return hostgrid_.triangle( element.impl().hostEntity() ).area();
    }

    template<int d = dim>
    std::enable_if_t< d == 3, FieldType >
    signedVolume_( const Entity& element ) const
    {
      return hostgrid_.tetrahedron( element.impl().hostEntity() ).volume();
    }

    bool adapt_( bool buildComponents = true )
    {
      if (insert_.size() == 0 && remove_.size() == 0)
        return false;

      std::vector<std::size_t> insertComponentIds;
      std::vector<std::size_t> removeComponentIds;
      static constexpr bool writeComponents = verbose_; // for debugging

      if ( buildComponents )
      {
        // build the connected components
        componentCount_ = connectedComponents_.size()+1; // 0 is the default component

        // we have to mark the cells repeatedly for the case that conflict zones overlap
        bool markAgain = true;
        std::size_t componentNumber;
        while( markAgain )
        {
          markAgain = false;
          insertComponentIds.clear();
          removeComponentIds.clear();
          insertComponentIds.reserve(insert_.size());
          removeComponentIds.reserve(remove_.size());

          // mark elements as mightVanish
          for ( const auto& ip : insert_ )
          {
            if ( ip.edgeId != IdType() )
              markAgain |= markElementsForInsertion_( ip.edge, componentNumber );
            else
              markAgain |= markElementForInsertion_( ip.point, componentNumber );

            insertComponentIds.push_back( componentNumber-1 );
          }

          // mark elements as mightVanish
          for ( const auto& vh : remove_ )
          {
            markAgain |= markElementsForRemoval_( vh, componentNumber );
            removeComponentIds.push_back( componentNumber-1 );
          }
        }

        buildConnectedComponents_();

        // plot connected components
        if( writeComponents )
          writeComponents_();
      }

      // actually insert the points
      std::vector<VertexHandle> newVertices;
      for ( const auto& ip : insert_ )
      {
        VertexHandle vh;
        bool connect = false;
        if ( ip.edgeId != IdType() )
        {
          auto eh = ip.edge.impl().hostEntity();

          // if edge id has changed, we have to update the edge handle
          if ( this->globalIdSet().id( ip.edge ) != ip.edgeId )
          {
            // add empty entry to newVertices to match indices with insertComponentIds
            newVertices.emplace_back();

            continue;
            // getEdge_( ip, eh ); // TODO this does not work correct yet
          }

          if ( ip.v0 != ip.edge.impl().template subEntity<dim>(0).impl().hostEntity()
            && ip.v0 != ip.edge.impl().template subEntity<dim>(1).impl().hostEntity() )
            connect = true;

          if ( ip.isInterface == true )
            vh = insertInInterface_( ip );
          else
            vh = insertInEdge_( ip.point, eh );

          vh->info().insertionLevel = ip.insertionLevel;
        }
        else
        {
          vh = insertInCell_( ip.point );

          // check if edge is really part of the triangulation
          if ( ip.v0 != VertexHandle() && !getHostGrid().tds().is_edge( ip.v0, vh ) ) // TODO 3D
          {
            // try again with half distance
            if constexpr (dimension == 2)
              hostgrid_.remove( vh );

            GlobalCoordinate x;
            x = makeFieldVector( ip.point );
            x -= makeFieldVector( ip.v0->point() );
            x *= 0.5;
            x += makeFieldVector( ip.v0->point() );
            vh = insertInCell_( makePoint( x ) );

            if ( !getHostGrid().tds().is_edge( ip.v0, vh ) ) // TODO 3D
            {
              // DUNE_THROW( GridError, "Edge added to interface is not part of the new triangulation: " << ip.v0->point() << " to " << vh->point() );
              std::cerr << "Error: Edge added to interface is not part of the new triangulation: " << ip.v0->point() << " to " << vh->point() << std::endl;
              continue;
            }
          }

          globalIdSet_->setNextId( vh );
          vh->info().insertionLevel = ip.insertionLevel;

          if ( ip.isInterface )
            connect = true;
        }

        // connect vertex and ip.v0 with interface
        if ( connect )
        {
          std::size_t id = vh->info().id;
          vh->info().isInterface = true;
          std::vector<std::size_t> ids;
          ids.push_back( id );
          ids.push_back( globalIdSet().id( entity( ip.v0 ) ).vt()[0] );
          std::sort(ids.begin(), ids.end());
          interfaceSegments_.insert( std::make_pair(ids, 1) ); // TODO: compute interface marker corresponding to ip.v0

          // pass this refinement information to the interface grid
          interfaceGrid_->markAsRefined( /*children*/ {ids}, ip.connectedcomponent );

          // set boundary segment to the same as ip.v0 if possible, default to 0
          auto v0id = interfaceGrid_->globalIdSet().id( interfaceGrid_->entity( ip.v0 ) ).vt()[0];
          auto it = interfaceGrid_->boundarySegments().find( { v0id } );
          if ( it != interfaceGrid_->boundarySegments().end() )
            interfaceGrid_->addBoundarySegment( { id }, it->second );
          else
            interfaceGrid_->addBoundarySegment( { id }, 0 );
        }

        // store vertex handles for later use
        newVertices.push_back( vh );
      }

      // actually remove the points
      int ci = 0;
      for ( const auto& vh : remove_ )
      {
        ElementOutput elements;
        if ( vh->info().isInterface )
          elements = removeFromInterface_( vh );
        else
          hostgrid_.removeAndGiveNewElements( vh, elements );

        // flag all elements inside conflict area as new and map connected component
        if ( buildComponents )
          markElementsAfterRemoval_( elements, removeComponentIds[ci] );
        ci++;
      }

      // first update ids
      setIds();
      interfaceGrid_->setIds();

      // then, update partitions
      if (comm().size() > 1)
        partitionHelper_.updatePartitions();

      // afterwards, update index sets
      setIndices();

      if ( buildComponents )
      {
        // flag incident elements as new and map connected component
        for ( std::size_t i = 0; i < newVertices.size(); ++i )
          if ( newVertices[i] != VertexHandle() )
            markElementsAfterInsertion_( newVertices[i], insertComponentIds[i] );
      }

      // update interface grid
      interfaceGrid_->setIndices();

      if ( buildComponents )
      {
        if( writeComponents )
          writeComponents_();
      }

      loadBalance();
      return newVertices.size() > 0;
    }

    template<int d = dim>
    std::enable_if_t< d == 2, void >
    getEdge_ ( const RefinementInsertionPoint &ip, EdgeHandle& eh ) const
    {
      assert( hostgrid_.is_edge( ip.v0, ip.v1 ) );
      hostgrid_.is_edge( ip.v0, ip.v1, eh.first, eh.second );
    }

    template<int d = dim>
    std::enable_if_t< d == 3, void >
    getEdge_ ( const RefinementInsertionPoint &ip, EdgeHandle& eh ) const
    {
      assert( hostgrid_.is_edge( ip.v0, ip.v1, eh.first, eh.second, eh.third ) );
      hostgrid_.is_edge( ip.v0, ip.v1, eh.first, eh.second, eh.third );
    }

    template<int d = dim>
    std::enable_if_t< d == 2, VertexHandle >
    insertInEdge_ ( const Point &point, const EdgeHandle& eh )
    {
      return hostgrid_.insert_in_edge( point, eh.first, eh.second );
    }

    template<int d = dim>
    std::enable_if_t< d == 3, VertexHandle >
    insertInEdge_ ( const Point &point, const EdgeHandle& eh )
    {
      return hostgrid_.insert_in_edge( point, eh.first, eh.second, eh.third );
    }

    template<int d = dim>
    std::enable_if_t< d == 2, VertexHandle >
    insertInCell_ ( const Point &point )
    {
      auto face = hostgrid_.locate( point );
      return hostgrid_.insert_in_face( point, face );
    }

    template<int d = dim>
    std::enable_if_t< d == 3, VertexHandle >
    insertInCell_ ( const Point &point )
    {
      auto cell = hostgrid_.locate( point );
      return hostgrid_.insert_in_cell( point, cell );
    }

  public:
    /** \brief Move interface vertices
     * \param shifts Vector that maps interface vertex indices to GlobalCoordinate
     */
    void moveInterface( const std::vector<GlobalCoordinate>& shifts )
    {
      const auto& iindexSet = this->interfaceGrid().leafIndexSet();
      assert( shifts.size() == iindexSet.size(dimension-1) );

      for( const auto& vertex : vertices( this->interfaceGrid().leafGridView() ) )
      {
        const VertexHandle& vh = vertex.impl().hostEntity();
        vh->point() = makePoint( vertex.geometry().center() + shifts[iindexSet.index(vertex)] );
      }
    }

    /** \brief Move vertices
     * \param shifts Vector that maps interface vertex indices to GlobalCoordinate
     */
    void moveVertices( const std::vector<GlobalCoordinate>& shifts )
    {
      const auto& indexSet = this->leafIndexSet();
      assert( shifts.size() == indexSet.size(dimension) );

      for( const auto& vertex : vertices( this->leafGridView() ) )
      {
        const VertexHandle& vh = vertex.impl().hostEntity();
        vh->point() = makePoint( vertex.geometry().center() + shifts[indexSet.index(vertex)] );
      }
    }

    //! Insert p into the triangulation and add a new interface segment between p and vertex
    template< typename Vertex >
    void addToInterface( const Vertex& vertex, const GlobalCoordinate& p )
    {
      RefinementInsertionPoint ip;
      ip.edgeId = IdType();
      ip.point = makePoint( p );
      ip.v0 = vertex.impl().hostEntity();
      ip.insertionLevel = vertex.impl().insertionLevel() + 1;
      ip.isInterface = true;

      // take some incident element as father
      for ( const auto& elem : incidentInterfaceElements( vertex ) )
        if (!elem.isNew())
        {
          ip.connectedcomponent = InterfaceGridConnectedComponent( elem );
          break;
        }

      insert_.push_back( ip );
      if (verbose_)
        std::cout << "Add vertex to interface: " << ip.point << std::endl;
    }

    //! Clean up refinement markers
    void postAdapt()
    {
      for( const Entity& entity : elements( this->leafGridView() ) )
      {
        mark( 0, entity );
        entity.impl().setIsNew( false );
        entity.impl().setWillVanish( false );
        entity.impl().hostEntity()->info().componentNumber = 0;
      }

      refineMarked_ = 0;
      coarsenMarked_ = 0;
      createdEntityConnectedComponentMap_.clear();
      vanishingEntityConnectedComponentMap_.clear();
      connectedComponents_.clear();

      insert_.clear();
      remove_.clear();
      inserted_.clear();
      removed_.clear();
    }

  private:
    void writeComponents_() const
    {
      static int performCount = 0;
      const auto& gridView = this->leafGridView();
      const auto& indexSet = this->leafIndexSet();

      std::vector<std::size_t> component( indexSet.size(0), 0 );
      std::vector<std::size_t> findcomponent( indexSet.size(0), 0 );
      std::vector<bool> isnew( indexSet.size(0), 0 );
      std::vector<bool> mightvanish( indexSet.size(0), 0 );

      for( const auto& e : elements( gridView ) )
      {
        component[indexSet.index(e)] = e.impl().hostEntity()->info().componentNumber;
        isnew[indexSet.index(e)] = e.isNew();
        mightvanish[indexSet.index(e)] = e.mightVanish();
      }

      VTKWriter<typename Traits::LeafGridView> vtkWriter( gridView );
      vtkWriter.addCellData( component, "component" );
      vtkWriter.addCellData( isnew, "isnew" );
      vtkWriter.addCellData( mightvanish, "mightvanish" );
      vtkWriter.write("mmesh-components-" + std::to_string( performCount ));

      VTKWriter<typename InterfaceGrid::LeafGridView> ivtkWriter( interfaceGrid_->leafGridView() );
      ivtkWriter.write("mmesh-components-interface-" + std::to_string( performCount++ ));
    }

    void buildConnectedComponents_()
    {
      for( const Entity& entity : elements( this->leafGridView() ) )
      {
        if( entity.mightVanish() )
        {
          const auto& hostEntity = entity.impl().hostEntity();
          const IdType id = this->globalIdSet().id( entity );

          bool found = false;
          const std::size_t componentId = hostEntity->info().componentNumber - 1;

          // check if component for entity exists already and add entity if necessary
          if( componentId < connectedComponents_.size() )
          {
            // sth. changed, we have to update the component to include this entity
            if( !entity.isNew() )
            {
              ConnectedComponent& component = connectedComponents_[ componentId ];
              if ( !component.hasEntity( entity ) )
                component.update( entity );
            }

            vanishingEntityConnectedComponentMap_.insert( std::make_pair( id, componentId ) );
          }
          // else create new component
          else
          {
            connectedComponents_[ componentId ] = ConnectedComponent( This(), entity );
            vanishingEntityConnectedComponentMap_.insert( std::make_pair( id, componentId ) );
          }
        }
      }
    }

    //! Insert a point to the triangulation and connect it to the vertices of the given interface element
    VertexHandle insertInInterface_( const RefinementInsertionPoint& ip )
    {
      assert( isInterface(ip.edge) );

      // get the vertex ids of the interface segment
      std::vector<std::size_t> ids;
      for( std::size_t i = 0; i < ip.edge.subEntities(dim); ++i )
        ids.push_back( globalIdSet().id( ip.edge.impl().template subEntity<dim>(i) ).vt()[0] );
      std::sort(ids.begin(), ids.end());

      // erase old interface segment
      std::size_t marker = interfaceSegments_[ids];
      interfaceSegments_.erase( ids );

      // insert the point
      auto eh = ip.edge.impl().hostEntity();
      const auto& vh = insertInEdge_( ip.point, eh );
      vh->info().isInterface = true;

      // get the (next) id for vh
      std::size_t id = globalIdSet_->setNextId( vh );

      // insert the new interface segments
      std::vector< std::vector<std::size_t> > allNewIds;
      for( int i = 0; i < dimension; ++i )
      {
        std::vector<std::size_t> newIds;
        newIds.push_back( id );
        for( int j = 0; j < dimension-1; ++j )
          newIds.push_back( ids[(i+j)%dimension] );

        std::sort(newIds.begin(), newIds.end());
        interfaceSegments_.insert( std::make_pair(newIds, marker) );
        allNewIds.push_back( newIds );
      }

      // pass this refinement information to the interface grid
      interfaceGrid_->markAsRefined( /*children*/ allNewIds, ip.connectedcomponent );

      return vh;
    }

    //! Remove a vertex from the triangulation and connect the corresponding interface elements
    template< int d = dimension >
    std::enable_if_t< d == 2, ElementOutput> removeFromInterface_( const VertexHandle& vh )
    {
      std::size_t id = vh->info().id;
      InterfaceGridConnectedComponent connectedComponent;

      // find and remove interface segments
      std::size_t marker = 1;
      std::vector<VertexHandle> otherVhs;
      for ( const auto& e : incidentInterfaceElements( interfaceGrid_->entity( vh ) ) )
      {
        // TODO: handle the case of overlapping components
        assert( !interfaceGrid_->hasConnectedComponent( e ) );
        connectedComponent.add( e );

        const auto v0 = e.template subEntity<1>(0).impl().hostEntity();
        const auto v1 = e.template subEntity<1>(1).impl().hostEntity();

        VertexHandle other = (v0 != vh) ? v0 : v1;

        std::vector<std::size_t> ids( 2 );
        ids[ 0 ] = id;
        ids[ 1 ] = other->info().id;
        std::sort(ids.begin(), ids.end());

        auto it = interfaceSegments_.find( ids );
        if ( it != interfaceSegments_.end() )
        {
          marker = interfaceSegments_[ ids ];
          interfaceSegments_.erase( it );
          otherVhs.push_back( other );
        }
      }

      if( otherVhs.size() != 2 )
        return {}; // otherwise, we remove a tip or a junction

      std::list<EdgeHandle> hole;
      hostgrid_.make_hole(vh, hole);

      ElementOutput elements;
      hostgrid_.remeshHoleConstrained(vh, hole, elements, otherVhs);

      // add the new interface segment
      std::vector<std::size_t> ids {{ otherVhs[0]->info().id, otherVhs[1]->info().id }};
      std::sort(ids.begin(), ids.end());
      interfaceSegments_.insert( std::make_pair(ids, marker) );

      // pass this refinement information to the interface grid
      interfaceGrid_->markAsRefined( /*children*/ { ids }, connectedComponent );

      return elements;
    }

    template< int d = dimension >
    std::enable_if_t< d == 3, ElementOutput> removeFromInterface_( const VertexHandle& vh )
    {
      DUNE_THROW( NotImplemented, "removeFromInterface() in 3d" );
    }

  public:
    /*@}*/

    /** \brief Size of the overlap on the leaf level */
    unsigned int overlapSize(int codim) const {
      return 0;
    }


    /** \brief Size of the ghost cell layer on the leaf level */
    unsigned int ghostSize(int codim) const {
      std::size_t size = 0;

      if (codim == 0)
        for (const auto& e : elements(this->leafGridView(), Partitions::ghost))
          size++;

      if (codim == 1)
        for (const auto& f : facets(this->leafGridView()))
          if (f.partitionType() == GhostEntity)
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
     */
    void loadBalance()
    {
      partitionHelper_.distribute();
      update();
    };

    /** \brief Distributes this grid over the available nodes in a distributed machine
     *
     * \param minlevel The coarsest grid level that gets distributed
     * \param maxlevel does currently get ignored
     */
    template<class T>
    bool loadBalance( const T& t ) {
      DUNE_THROW(NotImplemented, "MMesh::loadBalance(t)");
    };

    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement){
      DUNE_THROW(NotImplemented, "MMesh::loadBalance()");
    }


    /** \brief Communication */
    const Communication<Comm>& comm () const
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
        MMeshCommunication<GridImp, GridImp> communication( partitionHelper_ );
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


    // **********************************************************
    // End of Interface Methods
    // **********************************************************

    //! Get reference to the underlying CGAL triangulation.
    const HostGrid& getHostGrid() const
    {
      return hostgrid_;
    }

    //! Get non-const reference to the underlying CGAL triangulation.
    HostGrid& getHostGrid()
    {
      return hostgrid_;
    }

    //! Get reference to the interface grid.
    const InterfaceGrid& interfaceGrid() const
    {
      return *interfaceGrid_;
    }

    //! Get a non-const reference to the interface grid.
    InterfaceGrid& interfaceGrid()
    {
      return *interfaceGrid_;
    }

    //! Get a pointer to the interface grid
    const auto& interfaceGridPtr()
    {
      return interfaceGrid_;
    }

    const ConnectedComponent& getConnectedComponent( const Entity& entity ) const
    {
      const auto& it = createdEntityConnectedComponentMap_.find( this->globalIdSet().id( entity ) );
      assert( it->second < connectedComponents_.size() );
      assert( it != createdEntityConnectedComponentMap_.end() );
      assert( connectedComponents_.at( it->second ).size() > 0 );
      return connectedComponents_.at( it->second );
    }

    const RemeshingIndicator& indicator() const
    {
      return indicator_;
    }

    RemeshingIndicator& indicator()
    {
      return indicator_;
    }

    const auto& distance() const
    {
      return indicator_.distance();
    }

    int sequence() const
    {
      return sequence_;
    }

    const PartitionHelper<GridImp>& partitionHelper() const
    {
      return partitionHelper_;
    }

  private:
    // count how much elements where marked
    mutable int coarsenMarked_;
    mutable int refineMarked_;
    mutable std::size_t componentCount_;

    //! The storage of the connected components of entities
    std::unordered_map< IdType, ConnectedComponent > connectedComponents_;

    //! Maps the entities to its connected components
    std::unordered_map< IdType, std::size_t > vanishingEntityConnectedComponentMap_;
    std::unordered_map< IdType, std::size_t > createdEntityConnectedComponentMap_;

    Communication<Comm> comm_;
    PartitionHelper<GridImp> partitionHelper_;

    std::unique_ptr<MMeshLeafIndexSet<const GridImp>> leafIndexSet_;
    std::unique_ptr<MMeshGlobalIdSet<const GridImp>> globalIdSet_;
    std::shared_ptr<InterfaceGrid> interfaceGrid_;

    std::vector<RefinementInsertionPoint> insert_;
    std::unordered_set< IdType > inserted_;
    std::vector<VertexHandle> remove_;
    std::unordered_set< IdType > removed_;

    //! The host grid which contains the actual grid hierarchy structure
    HostGrid hostgrid_;
    BoundarySegments boundarySegments_, localBoundarySegments_;
    BoundaryIds boundaryIds_;
    InterfaceSegments interfaceSegments_;
    RemeshingIndicator indicator_;

    static const bool verbose_ = false;
    int sequence_ = 0;

  private:
    //! Flag all elements in conflict as mightVanish
    bool markElementsForInsertion_ ( const Edge& edge, std::size_t& componentNumber )
    {
      ElementOutput elements;
      getIncidentToEdge_( edge, elements );

      bool markAgain = getComponentNumber_( elements, componentNumber );

      // set componentNumber and mightVanish
      for(const auto& element : elements)
      {
        element->info().componentNumber = componentNumber;
        element->info().mightVanish = true;
      }

      return markAgain;
    }

    /// @cond
    template< int d = dimension >
    std::enable_if_t< d == 2, void > getIncidentToEdge_( const Edge& edge, ElementOutput& elements ) const
    {
      const auto& eh = edge.impl().hostEntity();
      elements.push_back( eh.first );
      elements.push_back( eh.first->neighbor( eh.second ) );
    }

    template< int d = dimension >
    std::enable_if_t< d == 3, void > getIncidentToEdge_( const Edge& edge, ElementOutput& elements ) const
    {
      const auto& eh = edge.impl().hostEntity();
      auto cit = this->getHostGrid().incident_cells( eh );
      for ( std::size_t i = 0; i < CGAL::circulator_size(cit); ++i, ++cit )
        elements.push_back( cit );
    }
    /// @endcond

    //! Flag element in conflict with point
    bool markElementForInsertion_ ( const Point& point, std::size_t& componentNumber )
    {
      ElementOutput elements;
      elements.push_back( this->getHostGrid().locate( point ) );

      bool markAgain = getComponentNumber_( elements, componentNumber );

      // set componentNumber and mightVanish
      for(const auto& element : elements)
      {
        element->info().componentNumber = componentNumber;
        element->info().mightVanish = true;
      }

      return markAgain;
    }

    //! Flag all incident elements as new
    void markElementsAfterInsertion_ ( const HostGridEntity<dimension>& vh, const std::size_t componentId )
    {
      for( const auto& element : incidentElements( entity( vh ) ) )
      {
        element.impl().hostEntity()->info().isNew = true;
        element.impl().hostEntity()->info().componentNumber = componentId + 1;
        const IdType id = this->globalIdSet().id( element );
        this->createdEntityConnectedComponentMap_.insert( std::make_pair( id, componentId ) );
      }
    }

    //! Flag all incident elements as mightVanish
    bool markElementsForRemoval_ ( const HostGridEntity<dimension>& vh, std::size_t& componentNumber )
    {
      ElementOutput elements;
      for( const auto& element : incidentElements( entity( vh ) ) )
        elements.push_back( element.impl().hostEntity() );

      bool markAgain = getComponentNumber_( elements, componentNumber );

      // set componentNumber and mightVanish
      for( const auto& element : elements )
      {
        element->info().componentNumber = componentNumber;
        element->info().mightVanish = true;
      }

      return markAgain;
    }

    //! Flag all elements in conflict as new
    void markElementsAfterRemoval_ ( ElementOutput& elements, const std::size_t componentId )
    {
      // flag all elements in conflict as mightVanish
      for(const auto& element : elements)
      {
        element->info().isNew = true;
        element->info().componentNumber = componentId + 1;

        const auto id = this->globalIdSet().id( this->entity( element ) );
        this->createdEntityConnectedComponentMap_.insert( std::make_pair( id, componentId ) );
      }
    }

    //! Get the component number for a given set of elements
    bool getComponentNumber_( const ElementOutput& elements, std::size_t& componentNumber ) const
    {
      bool markAgain = false;
      componentNumber = 0;

      // search for a componentNumber
      for( const auto& element : elements )
        if( element->info().componentNumber > 0 )
        {
          // check if we find different component numbers
          if( componentNumber != 0 && element->info().componentNumber != componentNumber )
          {
            markAgain = true;
            componentNumber = std::min( element->info().componentNumber, componentNumber );
          }
          else
            componentNumber = element->info().componentNumber;
        }

      // if no componentNumber was found, create a new one
      if ( componentNumber == 0 )
        componentNumber = componentCount_++;

      return markAgain;
    }

  }; // end Class MMesh


  // Capabilites of MMesh
  // --------------------

  /// @cond
  namespace Capabilities
  {
    /** \brief MMesh has only one geometry type for all entities
    \ingroup MMesh
    */
    template<class HostGrid, int dim>
    struct hasSingleGeometryType<MMesh<HostGrid, dim>>
    {
      static const bool v = true;
      static const unsigned int topologyId = Dune::GeometryType::simplex;
    };

    /** \brief MMesh can communicate on codim 0
    \ingroup MMesh
    */
    template<class HostGrid, int dim, int codim>
    struct canCommunicate<MMesh<HostGrid, dim>, codim>
    {
      static const bool v = false;
    };

    template<class HostGrid, int dim>
    struct canCommunicate<MMesh<HostGrid, dim>, 0>
    {
      static const bool v = true;
    };

    /** \brief MMesh has entities for some codimensions
     * \ingroup MMesh
     */
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

    /** \brief MMesh has conforming level grids
     * \ingroup MMesh
     */
    template<class HostGrid, int dim>
    struct isLevelwiseConforming<MMesh<HostGrid, dim>>
    {
      static const bool v = true;
    };

    /** \brief MMesh has conforming leaf grids when host grid has
     * \ingroup MMesh
     */
    template<class HostGrid, int dim>
    struct isLeafwiseConforming<MMesh<HostGrid, dim>>
    {
      static const bool v = true;
    };

    /** \brief MMesh has backup and restore facilities
     * \ingroup MMesh
     */
    template<class HostGrid, int dim>
    struct hasBackupRestoreFacilities<MMesh<HostGrid, dim>>
    {
      static const bool v = false;
    };

  } // end namespace Capabilities
  /// @endcond

} // namespace Dune


template<int codim, int dim, class HostGrid>
std::ostream& operator<< ( std::ostream& stream, const Dune::Entity<codim, dim, const Dune::MMesh<HostGrid, dim>, Dune::MMeshEntity >& entity )
{
  return stream << "MMeshEntity<codim=" << codim << ", dim=" << dim << ", center=(" << entity.geometry().center() << ")>";
}

#include "gridfactory.hh"
#include "structuredgridfactory.hh"
#include "dgfparser.hh"
#include "gmshgridfactory.hh"
#include "gmshreader.hh"
#include "../interface/grid.hh"
#include "../misc/capabilities.hh"
#include "../misc/persistentcontainer.hh"

#endif // DUNE_MMESH_GRID_MMESH_HH
