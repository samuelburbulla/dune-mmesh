// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_INTERSECTIONS_HH
#define DUNE_MMESH_GRID_INTERSECTIONS_HH

#include <unordered_set>

// Dune MMesh includes
#include <dune/mmesh/grid/entity.hh>
#include <dune/mmesh/grid/leafiterator.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>
#include <dune/mmesh/misc/twistutility.hh>

// local includes

/** \file
 * \brief The MMeshLeafIntersection and MMeshLevelIntersection classes
 */

namespace Dune
{

  // External forward declarations
  template< class Grid >
  struct HostGridAccess;

  /** \brief An intersection with a leaf neighbor element
   * \ingroup MMesh
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class MMeshLeafIntersection
  {
    friend class MMeshLeafIntersectionIterator<GridImp>;

    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::template HostGridEntity<0> HostGridEntity;
    typedef typename GridImp::template HostGridEntity<1> HostLeafIntersection;

  public:
    enum {dimension=GridImp::dimension};
    enum {dimensionworld=GridImp::dimensionworld};
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry::Implementation LocalGeometryImpl;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimworld> NormalVector;

    MMeshLeafIntersection()
    {}

    MMeshLeafIntersection(const GridImp* mMesh,
                          const HostGridEntity& hostEntity,
                          const int index)
      : mMesh_(mMesh)
      , hostIntersection_(hostEntity, index)
    {}

    MMeshLeafIntersection(const GridImp* mMesh,
                          HostLeafIntersection&& hostIntersection)
      : mMesh_(mMesh)
      , hostIntersection_(std::move(hostIntersection))
    {}

    //! returns true if the host entities are equal
    bool equals(const MMeshLeafIntersection& other) const
    {
      return hostIntersection_ == other.hostIntersection_;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const {
      return MMeshEntity<0,dim,GridImp>(mMesh_, hostIntersection_.first);
    }

    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const {
      assert( neighbor() );
      auto neighborHostEntity = (hostIntersection_.first)->neighbor(hostIntersection_.second);
      return MMeshEntity<0,dim,GridImp>(mMesh_, neighborHostEntity);
    }

    //! return true if intersection is with boundary.
    bool boundary () const {
      auto neighborHostEntity = (hostIntersection_.first)->neighbor(hostIntersection_.second);
      return mMesh_->getHostGrid().is_infinite(neighborHostEntity);
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const {
      return unitOuterNormal( FieldVector<ctype, dim-1> ( 0 ) );
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor() const {
      return !boundary();
    }

    //! return the boundary segment index
    std::size_t boundarySegmentIndex() const
    {
      HostGridEntity cell  = hostIntersection_.first;
      const auto& facetIdx = hostIntersection_.second;

      std::vector< std::size_t > vertices;
      for( std::size_t i = 0; i < dim; ++i )
        vertices.push_back( cell->vertex( (facetIdx+i+1)%(dim+1) )->info().id );
      std::sort(vertices.begin(), vertices.end());

      auto it = mMesh_->boundarySegments().find( vertices );
      if( it == mMesh_->boundarySegments().end() )
        return 0; // default to 0

      return it->second;
    }

    //! return the boundary id
    std::size_t boundaryId() const
    {
      auto it = mMesh_->boundaryIds().find( boundarySegmentIndex() );
      if( it == mMesh_->boundaryIds().end() )
        return 0; // default

      return it->second;
    }

    //! Return true if this is a conforming intersection
    bool conforming () const {
      // we are always conforming
      return true;
    }

    //! Geometry type of an intersection
    GeometryType type () const {
      return GeometryTypes::simplex(dim-1);
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      LocalGeometryImpl impl( indexInInside() );
      return LocalGeometry( impl );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      assert( neighbor() );
      LocalGeometryImpl impl( indexInOutside() );
      return LocalGeometry( impl );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    Geometry geometry () const
    {
      return Geometry( hostIntersection_ );
    }

    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const {
      return MMeshImpl::cgalFacetToDuneFacet<dim, HostLeafIntersection>( hostIntersection_ );
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const {
      const auto& neighbor = hostIntersection_.first->neighbor( hostIntersection_.second );
      const auto& second = mMesh_->getHostGrid().mirror_index( hostIntersection_.first, hostIntersection_.second );
      HostLeafIntersection facetFromOutside ( { neighbor, second } );
      return MMeshImpl::cgalFacetToDuneFacet<dim, HostLeafIntersection>( facetFromOutside );
    }

    //! return outer normal
    FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return integrationOuterNormal( local );
    }

    //! return outer normal multiplied by the integration element
    template< int d = dim >
    typename std::enable_if_t< d == 2, NormalVector >
    integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const
    {
      HostGridEntity face = hostIntersection_.first;
      const auto& edgeIdx = hostIntersection_.second;

      const auto& p1 = face->vertex( face->cw ( edgeIdx ) )->point();
      const auto& p2 = face->vertex( face->ccw( edgeIdx ) )->point();

      return NormalVector ( { p1.y() - p2.y(), - p1.x() + p2.x() } );
    }

    template< int d = dim >
    typename std::enable_if_t< d == 3, NormalVector >
    integrationOuterNormal (const FieldVector<ctype, dim-1>& local) const
    {
      HostGridEntity cell  = hostIntersection_.first;
      const auto& facetIdx = hostIntersection_.second;

      const auto& p1 = cell->vertex( (facetIdx+1)%4 )->point();
      const auto& p2 = cell->vertex( (facetIdx+2)%4 )->point();
      const auto& p3 = cell->vertex( (facetIdx+3)%4 )->point();

      return makeFieldVector(
        ( facetIdx % 2 == 0) ?
          CGAL::normal(p1, p2, p3) :
          CGAL::normal(p1, p3, p2)
      );
    }

    //! return unit outer normal
    FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      NormalVector n = integrationOuterNormal( local );
      n /= n.two_norm();
      return n;
    }

    const HostLeafIntersection& getHostIntersection() const
    {
      return hostIntersection_;
    }

  private:
    //! the host intersection
    const GridImp* mMesh_;
    HostLeafIntersection hostIntersection_;
  };

}  // namespace Dune

#endif
