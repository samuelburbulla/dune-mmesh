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

#include <dune/common/deprecated.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/version.hh>
#include <dune/common/to_unique_ptr.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/grid.hh>

#include <CGAL/utility.h>

// The components of the MMesh interface
#include "grid/incidentiterator.hh"
#include "grid/geometry.hh"
#include "grid/entity.hh"
#include "grid/entityseed.hh"
#include "grid/intersectioniterator.hh"
#include "grid/leafiterator.hh"
#include "grid/indexsets.hh"
#include "grid/hierarchiciterator.hh"
#include "grid/mmeshdefaults.hh"
#include "grid/pointfieldvector.hh"
#include "grid/rangegenerators.hh"
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

    typedef GridTraits<
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
        > Traits;
  };

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
    static constexpr int dimension = dim;
    using HostGridType = HostGrid;

    using Point = typename HostGrid::Point;
    using FieldType = typename Point::R::RT;

    using BoundarySegments = std::map< std::vector< unsigned int >, std::size_t >;

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

    //! the leaf iterator
    typedef typename Traits::template Codim<0>::LeafIterator LeafIterator;

    //! The type used to store coordinates, inherited from the HostGrid
    typedef FieldType ctype;

    //! The type used for coordinates
    typedef Dune::FieldVector<ctype, dimension> GlobalCoordinate;

    //! The type of the underlying entities
    template<int cd>
    using HostGridEntity = typename HostGridEntityChooser_<HostGridType, dimension, cd>::type;
    using VertexHandle = HostGridEntity<dimension>;

    //! The type of the element output
    using ElementOutput = std::list<HostGridEntity<0>>;

    //! The type used to store connected components of entities
    using Entity = typename Traits::template Codim<0>::Entity;
    using Edge = typename Traits::template Codim<dimension-1>::Entity;
    using Vertex = typename Traits::template Codim<dimension>::Entity;

    /** \brief Constructor
     *
     * \param hostgrid The host grid wrapped by the MMesh
     */
    explicit MMesh(HostGrid hostgrid, BoundarySegments boundarySegments)
     : hostgrid_(hostgrid),
       boundarySegments_(boundarySegments),
       numBoundarySegments_(boundarySegments.size())
    {
      leafIndexSet_ = std::make_unique<MMeshLeafIndexSet<const GridImp>>( This() );
      globalIdSet_ = std::make_unique<MMeshGlobalIdSet<const GridImp>>( This() );
      localIdSet_ = std::make_unique<MMeshLocalIdSet<const GridImp>>( This() );

      setIndices();
    }

    explicit MMesh(HostGrid hostgrid)
     : hostgrid_(hostgrid)
    {
      for ( const auto& element : elements( This()->leafGridView() ) )
        for ( const auto& intersection : intersections( This()->leafGridView(), element ) )
          if ( intersection.boundary() )
            numBoundarySegments_++;

      leafIndexSet_ = std::make_unique<MMeshLeafIndexSet<const GridImp>>( This() );
      globalIdSet_ = std::make_unique<MMeshGlobalIdSet<const GridImp>>( This() );
      localIdSet_ = std::make_unique<MMeshLocalIdSet<const GridImp>>( This() );

      setIndices();
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

    /** \brief returns the number of boundary segments within the macro grid
     */
    size_t numBoundarySegments () const {
      return numBoundarySegments_;
    }

    const BoundarySegments& boundarySegments() const {
      return boundarySegments_;
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
    int size (GeometryType type) const {
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
    const MMeshLeafIndexSet<const GridImp>& leafIndexSet() const {
      return *leafIndexSet_;
    }

    /** \brief Create Entity from EntitySeed */
    template < class EntitySeed >
    typename Traits::template Codim<EntitySeed::codimension>::Entity
    entity(const EntitySeed& seed) const
    {
      typedef MMeshEntity<
        EntitySeed::codimension,
        dimension,
        const typename Traits::Grid
        > EntityImp;

      auto hostEntity = seed.impl().hostEntity();
      assert( hostEntity != decltype(hostEntity)() );
      return EntityImp(This(), hostEntity);
    }

    //! Return the entity corresponding to a vertex handle
    Vertex entity(const HostGridEntity<dimension>& vertexHandle) const {
      return entity( typename Traits::template Codim<dimension>::EntitySeed( vertexHandle ) );
    }

    //! Iterator to first entity of given codim on level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lbegin (int level) const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This());
    }


    //! one past the end on This() level
    template<int codim>
    typename Traits::template Codim<codim>::LevelIterator lend (int level) const {
      return MMeshLeafIterator<codim,All_Partition, const GridImp>(This(), true);
    }


    //! Iterator to first entity of given codim on level
    template<int codim, PartitionIteratorType PiType>
    typename Traits::template Codim<codim>::template Partition<PiType>::LevelIterator lbegin (int level) const {
      return MMeshLeafIterator<codim,PiType, const GridImp>(This());
    }


    //! one past the end on This() level
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

    /** \brief returns false, if at least one entity is marked for adaption */
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
    void loadBalance() {};

    template<class T>
    bool loadBalance( const T& t ) { return false; };

    void loadBalance(int strategy, int minlevel, int depth, int maxlevel, int minelement) {}


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

    const HostGrid& getHostGrid() const
    {
      return hostgrid_;
    }

    HostGrid& getHostGrid()
    {
      return hostgrid_;
    }

  private:
    CollectiveCommunication< GridImp > ccobj;

  public:
    //! compute the grid indices and ids
    void setIndices()
    {
      setIds();
      leafIndexSet_->update(This());
    }

  protected:
    //! compute the grid indices and ids
    void setIds()
    {
      localIdSet_->update(This());
      globalIdSet_->update(This());
    }

  private:
    std::unique_ptr<MMeshLeafIndexSet<const GridImp>> leafIndexSet_;

    std::unique_ptr<MMeshGlobalIdSet<const GridImp>> globalIdSet_;

    std::unique_ptr<MMeshLocalIdSet<const GridImp>> localIdSet_;

  protected:
    //! The host grid which contains the actual grid hierarchy structure
    HostGrid hostgrid_;

    BoundarySegments boundarySegments_;

    std::size_t numBoundarySegments_;

  }; // end Class MMesh


  namespace Capabilities
  {
    /** \brief has entities for some codimensions
     * \ingroup MMesh
     */
    template<class HostGrid, int dim, int codim>
    struct hasEntity<MMesh<HostGrid, dim>, codim>
    {
      static const bool v = (codim == 0 || codim == dim);
    };

    template<class HostGrid, int dim, int codim>
    struct hasEntityIterator<MMesh<HostGrid, dim>, codim>
    {
      static const bool v = (codim == 0 || codim == dim);
    };

    /** \brief has conforming level grids
     * \ingroup MMesh
     */
    template<class HostGrid, int dim>
    struct isLevelwiseConforming<MMesh<HostGrid, dim>>
    {
      static const bool v = true;
    };

    /** \brief has conforming leaf grids when host grid has
     * \ingroup MMesh
     */
    template<class HostGrid, int dim>
    struct isLeafwiseConforming<MMesh<HostGrid, dim>>
    {
      static const bool v = true;
    };
  } // end namespace Capabilities

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

#endif // DUNE_MMESH_MMESH_HH
