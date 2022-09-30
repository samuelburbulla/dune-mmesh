// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MMESH_GRID_IMPLICITGRIDFACTORY_HH
#define DUNE_MMESH_GRID_IMPLICITGRIDFACTORY_HH

// System includes
#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>
#include <unordered_set>

// Dune includes
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/common/boundarysegment.hh>

// MMesh includes
#include <dune/mmesh/grid/mmesh.hh>
#include <dune/mmesh/grid/common.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>

namespace Dune
{

  /** \brief specialization of the implicit GridFactory for MMesh
   *
   *  \ingroup GridFactory
   */

  //! The implicit grid factory for MMesh
  template< class Grid >
  class MMeshImplicitGridFactory
    : public GridFactoryInterface< Grid >
  {
    typedef MMeshImplicitGridFactory< Grid > This;

  public:
    //! type of (scalar) coordinates
    typedef typename Grid::ctype ctype;
    //! type of the hostgrid
    typedef typename Grid::HostGridType HostGrid;

    //! dimension of the grid
    static const int dimension = Grid::dimension;
    //! dimension of the world
    static const int dimensionworld = Grid::dimensionworld;

    //! type of vector for world coordinates
    typedef FieldVector< ctype, dimensionworld > WorldVector;
    //! type of matrix from world coordinates to world coordinates
    typedef FieldMatrix< ctype, dimensionworld, dimensionworld > WorldMatrix;

    //! type of a Dune boundary segment
    typedef Dune::BoundarySegment< dimension, dimensionworld > BoundarySegment;
    //! type of an id
    typedef typename Grid::IdType IdType;

    //! type of the boundary segment id map
    typedef std::unordered_map< IdType, std::size_t > BoundarySegments;
    typedef std::unordered_map< std::size_t, std::size_t > BoundaryIds;

    //! type of the interface segment set
    typedef std::unordered_map< IdType, std::size_t > InterfaceSegments;

    template< int codim >
    struct Codim
    {
      typedef typename Grid::template Codim< codim >::Entity Entity;
    };

  private:
    typedef typename HostGrid::Point Point;
    typedef typename HostGrid::Vertex_handle Vertex_handle;
    typedef typename Grid::template HostGridEntity<0> Element_handle;
    typedef typename Grid::template HostGridEntity<1> Face_handle;

  public:
    //! are boundary ids supported by this factory?
    static const bool supportsBoundaryIds = true;
    //! the factory is not able to create periodic meshes
    static const bool supportPeriodicity = false;

    /** default constructor */
    MMeshImplicitGridFactory ()
    {}

    /** \brief insert an element into the macro grid
     *
     *  \param[in]  type      GeometryType of the new element
     *  \param[in]  v         indices of the element vertices (starting with 0)
     *
     *  \note The implicit grid factory ignores the insertion of specific elements.
     *  The CGAL backend constructs the elements by inserting vertices implicitly.
     */
    void insertElement ( const GeometryType &type,
                         const std::vector< unsigned int > &v )
    {
      Element_handle fh;
      if( isElement( v, fh ) )
        fh->info().insertionIndex = countElements;

      // Increase element count in each case
      countElements++;
    };

    void insertElement ( const GeometryType &type,
                         const std::vector< unsigned int > &v,
                         const size_t domainMarker )
    {
      insertElement( type, v );
    }

    /** \brief Returns if there is a face with the given vertices in the triangulation
     *  \ingroup 2D
     *
     *  \param[in]  v        indices of the element vertices
     */
    template< int d = dimension >
    std::enable_if_t< d == 2, bool >
    isElement( const std::vector< unsigned int >& v, Element_handle& fh ) const
    {
        assert( v.size() == dimension+1 );
        return tr_.is_face( vhs_[v[0]], vhs_[v[1]], vhs_[v[2]], fh );
    }

    /** \brief Returns if there is a cell with the given vertices in the triangulation
     *  \ingroup 3D
     *
     *  \param[in]  v        indices of the element vertices
     */
    template< int d = dimension >
    std::enable_if_t< d == 3, bool >
    isElement( const std::vector< unsigned int >& v, Element_handle& fh ) const
    {
        assert( v.size() == dimension+1 );
        return tr_.is_cell( vhs_[v[0]], vhs_[v[1]], vhs_[v[2]], vhs_[v[3]], fh );
    }

    /** \brief insert boundary segment
     *
     *  \param[in]  vertices  Vertices
     */
    void insertBoundarySegment(const std::vector<unsigned int>& vertices)
    {
      assert( vertices.size() == dimension );

      std::vector< std::size_t > sorted_vertices;
      for ( const auto& v : vertices )
        sorted_vertices.push_back( vhs_[v]->info().id );
      std::sort(sorted_vertices.begin(), sorted_vertices.end());

      if( boundarySegments_.find( sorted_vertices ) != boundarySegments_.end() )
        DUNE_THROW( GridError, "A boundary segment was inserted twice." );

      boundarySegments_.insert( std::make_pair( sorted_vertices, countBoundarySegments++ ) );
    }

    void insertBoundarySegment(const std::vector<unsigned int>& vertices,
                               const std::shared_ptr<Dune::BoundarySegment<dimension,dimension>>& boundarySegment)
    {
      DUNE_THROW( NotImplemented, "insertBoundarySegments with Dune::BoundarySegment" );
    }

    void insertInterfaceBoundarySegment ( const std::vector< unsigned int >& vertices )
    {
      assert( vertices.size() == dimension-1 );

      std::vector< std::size_t > sorted_vertices;
      for ( const auto& v : vertices )
        sorted_vertices.push_back( vhs_[v]->info().id );
      std::sort(sorted_vertices.begin(), sorted_vertices.end());

      if( boundarySegments_.find( sorted_vertices ) != boundarySegments_.end() )
          DUNE_THROW( GridError, "A boundary segment was inserted twice." );

      interfaceBoundarySegments_.insert( std::make_pair( sorted_vertices, countInterfaceBoundarySegments++ ) );
    }

    /** \brief Insert a vertex into the macro grid
     *
     *  \param[in]  pos  position of the vertex (in world coordinates)
     *  \note This method assumes that the vertices are inserted consecutively
     *        with respect to their index.
     */
    void insertVertex ( const WorldVector &pos )
    {
      // Insert vertex
      Vertex_handle vh = tr_.insert( makePoint( pos ) );

      // Store insertion index
      vh->info().id = countVertices;
      vh->info().idWasSet = true;

      // Increase vertex counter
      countVertices++;

      // Store vertex handle for later use
      vhs_.push_back( vh );
    }

    /** \brief insert an interface into the macro grid
     *
     *  \param[in]  vertices  indices of the interface vertices (starting with 0)
     *  \param[in]  marker    marker value of the interface segment (default 1)
     */
    void insertInterface ( const std::vector< unsigned int > &vertices, const std::size_t marker = 1 )
    {
      assert( vertices.size() == dimension );

      std::vector< std::size_t > sorted_vertices;
      for( const auto& v : vertices )
        sorted_vertices.push_back( vhs_[v]->info().id );
      std::sort(sorted_vertices.begin(), sorted_vertices.end());
      interfaceSegments_.insert( std::make_pair( sorted_vertices, marker ) );
    }

    /** \brief return index of inserted vertex within the macro grid
     *
     *  \param[in]  pos  position of the vertex (in world coordinates)
     */
    unsigned int insertionIndex ( const typename Codim<0>::Entity &entity ) const
    {
      return entity.impl().hostEntity()->info().insertionIndex;
    }

    /** \brief return insertion index of entity
     *
     *  \param[in]  entity  Entity
     */
    unsigned int insertionIndex ( const typename Codim< dimension >::Entity &entity ) const
    {
      return entity.impl().hostEntity()->info().id;
    }

    /** \brief return insertion index of boundary intersection
     *
     *  \param[in]  intersection  Leaf intersection
     */
    unsigned int insertionIndex ( const typename Grid::LeafIntersection &intersection ) const
    {
      return intersection.impl().boundaryId();
    }

    //! returns the boundary segment to index map
    const BoundarySegments& boundarySegments() const
    {
      return boundarySegments_;
    }

    //! returns the boundary segment index to boundary id map
    const BoundaryIds& boundaryIds() const
    {
      return boundaryIds_;
    }

    //! add a boundary id
    void addBoundaryId( std::size_t boundarySegmentIndex, std::size_t boundaryId )
    {
      boundaryIds_.insert( std::make_pair( boundarySegmentIndex, boundaryId ) );
    }

    /** \brief finalize grid creation and hand over the grid
     *
     *  This version of createGrid is original to the MMesh grid factroy,
     *  allowing to specity a grid name.
     *
     *  \returns a pointer to the newly created grid
     */
    std::unique_ptr<Grid> createGrid ()
    {
      return std::make_unique<Grid> (
        std::move(tr_),
        std::move(boundarySegments_),
        std::move(interfaceBoundarySegments_),
        std::move(boundaryIds_),
        std::move(interfaceSegments_)
      );
    }

  private:
    /** \brief Convert FieldVector to CGAL Point
     *  \ingroup 2D
     */
    template< int dim = dimension >
    std::enable_if_t< dim == 2, Point >
    makePoint( const WorldVector& v ) const
    {
      return Point ( v[ 0 ], v[ 1 ] );
    }

    /** \brief Convert FieldVector to CGAL Point
     *  \ingroup 3D
     */
    template< int dim = dimension >
    std::enable_if_t< dim == 3, Point >
    makePoint( const WorldVector& v ) const
    {
      return Point ( v[ 0 ], v[ 1 ], v[ 2 ] );
    }

    //! Private members
    HostGrid tr_;
    std::vector< Vertex_handle > vhs_;
    BoundarySegments boundarySegments_, interfaceBoundarySegments_;
    BoundaryIds boundaryIds_;
    InterfaceSegments interfaceSegments_;
    std::size_t countVertices = 0, countElements = 0;
    std::size_t countBoundarySegments = 0, countInterfaceBoundarySegments = 0;
  };

} // end namespace Dune

#endif
