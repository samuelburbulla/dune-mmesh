// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_INTERSECTIONS_HH
#define DUNE_MMESH_INTERFACE_INTERSECTIONS_HH

#include <unordered_set>

// Dune MMesh includes
#include <dune/mmesh/grid/entity.hh>
#include <dune/mmesh/grid/leafiterator.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>

// local includes

/** \file
 * \brief The MMeshInterfaceGridLeafIntersection and MMeshLevelIntersection classes
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
  class MMeshInterfaceGridLeafIntersection
  {
    friend class MMeshInterfaceGridLeafIntersectionIterator<GridImp>;

    friend struct HostGridAccess< typename std::remove_const< GridImp >::type >;

    enum {dimension=GridImp::dimension};

    enum {dimensionworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::template MMeshInterfaceEntity<0> MMeshInterfaceEntity;
    typedef typename GridImp::template MMeshInterfaceEntity<1> HostLeafIntersection;

    using LocalIndexMap = std::unordered_map< std::size_t, std::size_t >;
    typedef typename GridImp::MMeshType::template Codim<dimensionworld>::Entity MMeshVertex;

  public:
    typedef typename GridImp::template Codim<1>::Geometry Geometry;
    typedef typename GridImp::template Codim<1>::LocalGeometry LocalGeometry;
    typedef typename GridImp::template Codim<0>::Entity Entity;
    typedef FieldVector<ctype, dimensionworld> NormalVector;

    MMeshInterfaceGridLeafIntersection()
    {}

    MMeshInterfaceGridLeafIntersection(const GridImp* grid,
                          const MMeshInterfaceEntity& hostEntity,
                          const std::size_t index,
                          const std::size_t nbIdx)
      : grid_(grid),
        interfaceEntity_(hostEntity),
        index_(index),
        nbIdx_(nbIdx)
    {
      const auto& indexSet = grid_->leafIndexSet();

      std::array< std::size_t, dimension > ids;
      try {
        for( int i = 0; i < dimension; ++i )
          ids[i] = indexSet.vertexIndexMap().at(
            interfaceEntity_.first->vertex((interfaceEntity_.second+i+( (i == 1 && index_ == 2) ? index_+2 : index_+1 ))%(dimensionworld+1))->info().index
          );
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }
      std::sort(ids.begin(), ids.end());
      try {
        localIndexMap_ = indexSet.indexMap().at( ids );
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }
      assert( 0 <= nbIdx_ && ( nbIdx_ < numOutside() || boundary() ) );
    }

    //! returns true if the host entities are equal
    bool equals(const MMeshInterfaceGridLeafIntersection& other) const
    {
      return interfaceEntity_ == other.interfaceEntity_
          && index_ == other.index_
          && nbIdx_ == other.nbIdx_;
    }

    //! return Entity on the inside of this intersection
    //! (that is the Entity where we started this Iterator)
    Entity inside() const {
      return MMeshInterfaceGridEntity<0,dimension,GridImp>(grid_, interfaceEntity_);
    }

    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    std::size_t numOutside() const {
      return localIndexMap_.size() - 1;
    }

    //! return Entity on the outside of this intersection
    //! (that is the neighboring Entity)
    Entity outside() const
    {
      assert( neighbor() );

      const auto& vertexIndexMap = grid_->leafIndexSet().vertexIndexMap();

      std::size_t myLastVertexIndex;
      try {
        myLastVertexIndex = vertexIndexMap.at(
          interfaceEntity_.first->vertex((interfaceEntity_.second+((index_ == 0) ? index_+dimensionworld : index_))%(dimensionworld+1))->info().index
        );
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }

      // get the i-th index map entry
      auto it = localIndexMap_.begin();
      // exclude myself
      while ( it->first == myLastVertexIndex )
        ++it;

      for( std::size_t count = 0; count < nbIdx_; count++ )
      {
        ++it;

        // exclude myself
        if ( it->first == myLastVertexIndex )
          count--;
      }

      // obtain neighbor by searching in incident elements of intersection vertex
      const auto& vertex0 = interfaceEntity_.first->vertex((interfaceEntity_.second+index_+1)%(dimensionworld+1));

      std::array< std::size_t, dimensionworld > vIdx;
      try {
        vIdx[0] = vertexIndexMap.at( vertex0->info().index );
      } catch (std::exception &e) {
        DUNE_THROW(InvalidStateException, e.what());
      }
      vIdx[1] = it->first;

      if ( dimensionworld == 3 )
        try {
          vIdx[2] = vertexIndexMap.at( interfaceEntity_.first->vertex((interfaceEntity_.second+((index_ < 2) ? index_+2 : 1))%(dimensionworld+1))->info().index );
        } catch (std::exception &e) {
          DUNE_THROW(InvalidStateException, e.what());
        }
      std::sort(vIdx.begin(), vIdx.end());

      // search in incident facets for the right entity
      for( const auto& facet : incidentFacets( MMeshVertex {{ &grid_->getMMesh(), vertex0 }} ) )
      {
        std::array< std::size_t, dimensionworld > tmpIdx;
        bool notInterface = false;
        for( std::size_t i = 0; i < facet.subEntities(dimensionworld); ++i )
        {
          std::size_t idx = facet.impl().template subEntity<dimensionworld>(i).impl().hostEntity()->info().index;
          auto it = vertexIndexMap.find( idx );
          if( it != vertexIndexMap.end() )
            tmpIdx[i] = it->second;
          else
            notInterface = true;
        }
        if( notInterface )
          continue;
        std::sort(tmpIdx.begin(), tmpIdx.end());

        if ( vIdx == tmpIdx )
          return MMeshInterfaceGridEntity<0,dimension,GridImp>(grid_, facet.impl().hostEntity());
      }
      DUNE_THROW( InvalidStateException, "Neighbor could not be determined!" );
    }

    //! return true if intersection is with boundary.
    bool boundary () const {
      return localIndexMap_.size() == 1;
    }

    /** \brief Return unit outer normal (length == 1)
     *
     *   The returned vector is the normal at the center() of the
     *     intersection's geometry.
     *       It is scaled to have unit length. */
    NormalVector centerUnitOuterNormal () const {
      return unitOuterNormal( FieldVector<ctype, dimension-1> ( 0 ) );
    }

    //! return true if across the edge an neighbor on this level exists
    bool neighbor() const {
      return !boundary();
    }

    //! return the boundary segment index
    size_t boundarySegmentIndex() const
    {
      auto face = interfaceEntity_.first;
      const auto& edgeIdx = interfaceEntity_.second;

      std::vector< std::size_t > vertices;
      vertices.push_back( face->vertex( (edgeIdx+index_+1)%(dimensionworld+1) )->info().id );

      if( dimension == 2 )
      {
        vertices.push_back( face->vertex( (edgeIdx+((index_ < 2) ? index_+2 : 1))%(dimensionworld+1) )->info().id );
        std::sort(vertices.begin(), vertices.end());
      }

      auto it = grid_->boundarySegments().find( vertices );
      if( it == grid_->boundarySegments().end() )
      {
        std::cerr << "BoundarySegmentIndex was not found, default to 0." << std::endl;
        return 0;
      }

      return it->second;
    }

    //! Return true if this is a conforming intersection
    bool conforming () const {
      // we are always conforming
      return true;
    }

    //! Geometry type of an intersection
    GeometryType type () const {
      return GeometryTypes::simplex(dimension-1);
    }

    //! intersection of codimension 1 of this neighbor with element where
    //! iteration started.
    //! Here returned element is in LOCAL coordinates of the element
    //! where iteration started.
    LocalGeometry geometryInInside () const
    {
      return LocalGeometry( index_ );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in LOCAL coordinates of neighbor
    LocalGeometry geometryInOutside () const
    {
      assert( neighbor() );
      return LocalGeometry( indexInOutside() );
    }

    //! intersection of codimension 1 of this neighbor with element where iteration started.
    //! Here returned element is in GLOBAL coordinates of the element where iteration started.
    template< int d = dimension >
    typename std::enable_if_t< d == 1, Geometry >
    geometry () const
    {
      return Geometry( interfaceEntity_.first->vertex( (interfaceEntity_.second+index_+1)%(dimensionworld+1) ) );
    }

    template< int d = dimension >
    typename std::enable_if_t< d == 2, Geometry >
    geometry () const
    {
      int vIdx0 = (interfaceEntity_.second+index_+1)%(dimensionworld+1);
      int vIdx1 = (interfaceEntity_.second+((index_ < 2) ? index_+2 : 1))%(dimensionworld+1);
      return Geometry( {{ interfaceEntity_.first, vIdx0, vIdx1 }} );
    }

    //! local number of codim 1 entity in self where intersection is contained in
    int indexInInside () const
    {
      return index_;
    }

    //! local number of codim 1 entity in neighbor where intersection is contained
    int indexInOutside () const
    {
      assert( dimensionworld == 2 );
      const auto neighbor = outside();
      const auto& vertex0 = interfaceEntity_.first->vertex((interfaceEntity_.second+index_+1)%(dimensionworld+1));

      if ( neighbor.template subEntity<dimension>( 0 ).impl().hostEntity() == vertex0 )
        return 0;
      if ( neighbor.template subEntity<dimension>( 1 ).impl().hostEntity() == vertex0 )
        return 1;

      DUNE_THROW( InvalidStateException, "indexInOutside() could not be determined!" );
    }

    //! return outer normal
    FieldVector<ctype, GridImp::dimensionworld> outerNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      return integrationOuterNormal( local );
    }

    //! return outer normal multiplied by the integration element
    template< int d = dimension >
    typename std::enable_if_t< d == 1, NormalVector >
    integrationOuterNormal (const FieldVector<ctype, dimension-1>& local) const
    {
      auto face = interfaceEntity_.first;
      const auto& edgeIdx = interfaceEntity_.second;

      const auto& p1 = face->vertex( (edgeIdx+index_+1)%(dimensionworld+1) )->point();
      const auto& p2 = face->vertex( (edgeIdx-index_+2)%(dimensionworld+1) )->point();

      NormalVector n ( { p1.x() - p2.x(), p1.y() - p2.y() } );
      n /= n.two_norm();
      return n;
    }

    template< int d = dimension >
    typename std::enable_if_t< d == 2, NormalVector >
    integrationOuterNormal (const FieldVector<ctype, dimension-1>& local) const
    {
      const auto& cell = interfaceEntity_.first;
      const auto& j = interfaceEntity_.second;
      auto a = makeFieldVector( cell->vertex( (j+index_+1)%(dimensionworld+1) )->point() );
      auto b = makeFieldVector( cell->vertex( (j+((index_ == 2) ? index_+3 : index_+2))%(dimensionworld+1) )->point() );
      auto v = makeFieldVector( cell->vertex( (j+((index_ == 0) ? index_+3 : index_))%(dimensionworld+1) )->point() );

      // return vector that is orthogonal to edge a-b and triangle a-b-v
      NormalVector ab = b - a;
      NormalVector vb = b - v;

      NormalVector elementNormal;
      elementNormal[0] = ab[1]*vb[2] - ab[2]*vb[1];
      elementNormal[1] = ab[2]*vb[0] - ab[0]*vb[2];
      elementNormal[2] = ab[0]*vb[1] - ab[1]*vb[0];

      NormalVector outerNormal;
      outerNormal[0] = ab[1]*elementNormal[2] - ab[2]*elementNormal[1];
      outerNormal[1] = ab[2]*elementNormal[0] - ab[0]*elementNormal[2];
      outerNormal[2] = ab[0]*elementNormal[1] - ab[1]*elementNormal[0];

      if( outerNormal * vb < 0.0 )
        outerNormal *= -1.0;

      return outerNormal;
    }

    //! return unit outer normal
    FieldVector<ctype, GridImp::dimensionworld> unitOuterNormal (const FieldVector<ctype, GridImp::dimension-1>& local) const {
      NormalVector n = integrationOuterNormal( local );
      n /= n.two_norm();
      return n;
    }

    const auto getHostVertex() const
    {
      return interfaceEntity_.first->vertex((interfaceEntity_.second+index_+1)%(dimensionworld+1));
    }

    const MMeshInterfaceEntity& getHostIntersection() const
    {
      return interfaceEntity_;
    }

  private:
    //! the host intersection
    const GridImp* grid_;
    MMeshInterfaceEntity interfaceEntity_;
    std::size_t index_;
    std::size_t nbIdx_;
    LocalIndexMap localIndexMap_;
  };

}  // namespace Dune

#endif
