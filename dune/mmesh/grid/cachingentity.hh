// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_CACHINGENTITY_HH
#define DUNE_MMESH_GRID_CACHINGENTITY_HH

/** \file
 * \brief The MMeshCachingEntity class
 */

#include <set>

// Dune includes
#include <dune/grid/common/grid.hh>

// MMesh includes
#include <dune/mmesh/grid/multiid.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>
#include <dune/mmesh/grid/polygoncutting.hh>

// CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/number_utils.h>
#include <CGAL/intersections.h>
#include <CGAL/Point_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/squared_distance_3.h>

namespace Dune
{

  //**********************************************************************
  //
  // --MMeshCachingEntity
  // --Entity
  //
  /** \brief The implementation of caching entities in a MMesh
   *   \ingroup MMesh
   *  The caching entity copys the CGAL face object instead of holding a Face_handle pointer.
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   *
   */
  template<int dim, class GridImp>
  class MMeshCachingEntity<0,dim,GridImp> :
    public MMeshEntity<0,dim,GridImp>
  {
    template <class GridImp_>
    friend class MMeshLeafIndexSet;

    template <class GridImp_>
    friend class MMeshLocalIdSet;

    template <class GridImp_>
    friend class MMeshGlobalIdSet;

  private:
    // this type
    typedef MMeshCachingEntity<0,dim,GridImp> ThisType;

    // base type
    typedef MMeshEntity<0,dim,GridImp> BaseType;

    // type of scalars
    typedef typename GridImp::ctype ctype;

    // type of the host grid
    typedef typename GridImp::HostGridType HostGrid;

    // equivalent entity in the host grid as pointer
    typedef typename GridImp::template HostGridEntity<0> HostGridEntity;

    // standard MMesh entity implementation
    typedef typename GridImp::template Codim<0>::Entity MMeshEntityType;

    // type of ids
    typedef MMeshImpl::MultiId IdType;

  public:
    // geometry type
    typedef AffineGeometry<ctype, dim, dim> Geometry;

    // local geometry type
    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    // type of global coordinate
    typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

    MMeshCachingEntity() = delete;

    explicit MMeshCachingEntity(const GridImp* mMesh, const HostGridEntity& hostEntity)
      : BaseType(mMesh, hostEntity, mMesh->globalIdSet().id( mMesh->entity( hostEntity ) ))
    {
      for( int i = 0; i < dim+1; ++i )
        this->vertex_[i] = makeFieldVector( hostEntity->vertex(i)->point() );
    }

    //! returns true if host entities are equal
    bool equals(const MMeshCachingEntity& other) const
    {
      return this->id_ == other.id_;
    }

    //! returns true if host entities are equal
    bool operator==(const MMeshCachingEntity& other) const
    {
      return this->equals(other);
    }

    //! returns true if caching entity has same id like mmesh entity
    bool operator==(const MMeshEntityType& entity) const
    {
      return this->id_ == this->mMesh_->globalIdSet().id( entity );
    }

    //! returns true if id of other is greater
    bool operator<(const MMeshCachingEntity& other) const
    {
      return this->id_ < other.id_;
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return true;
    }

    //! returns true if this entity is new after adaptation
    const bool isNew () const
    {
      return false;
    }

    //! returns true if this entity will vanish after adaptation
    const bool mightVanish () const
    {
      return true;
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
      return Geometry( GeometryTypes::simplex(dim), this->vertex_ );
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

    //! returns true if Entity has no children
    bool isLeaf() const {
      return false;
    }

    //calculates the intersection volume with another MMesh entity
    template<int d = dim>
    std::enable_if_t<d == 2, ctype>
    intersectionVolume ( const MMeshEntityType& entity ) const
    {
      std::array<GlobalCoordinate, 3> entityPoints;

      for ( int i = 0; i < 3; ++i )
      {
        entityPoints[i] =
          makeFieldVector( entity.impl().hostEntity()->vertex(i)->point() );
      }

      using PC = Dune::PolygonCutting<ctype, GlobalCoordinate>;
      return PC::polygonArea( PC::sutherlandHodgman(this->vertex_, entityPoints) );
    }

    template<int d = dim>
    std::enable_if_t<d == 3, ctype>
    intersectionVolume( const MMeshEntityType& entity ) const
    {
      // inexact types
      using Epick = CGAL::Exact_predicates_inexact_constructions_kernel;
      using Point = CGAL::Point_3<Epick>;
      using Triangle = CGAL::Triangle_3<Epick>;
      using Tetrahedron = CGAL::Tetrahedron_3<Epick>;

      // exact types
      using Epeck = CGAL::Exact_predicates_exact_constructions_kernel;
      using ExactPoint = CGAL::Point_3<Epeck>;
      using ExactSegment = CGAL::Segment_3<Epeck>;
      using ExactTriangle = CGAL::Triangle_3<Epeck>;

      // type converters
      typedef CGAL::Cartesian_converter<Epick, Epeck> Epick_to_Epeck;
      typedef CGAL::Cartesian_converter<Epeck, Epick> Epeck_to_Epick;
      Epick_to_Epeck to_exact;
      Epeck_to_Epick to_inexact;

      // store and convert points
      std::array<Point, 4> piA, piB;
      std::array<ExactPoint, 4> peA, peB;
      for ( int i = 0; i < 4; ++i )
      {
        const auto& p1 = makePoint( this->vertex_[i] );
        piA[i] = p1;
        peA[i] = to_exact( p1 );

        const auto& p2 = entity.impl().hostEntity()->vertex(i)->point();
        piB[i] = p2;
        peB[i] = to_exact( p2 );
      }

      // First, check if at least the bounding boxes intersect
      CGAL::Bbox_3 bbox1 = bbox_3(peA.begin(), peA.end());
      CGAL::Bbox_3 bbox2 = bbox_3(peB.begin(), peB.end());
      if( !CGAL::do_overlap( bbox1, bbox2 ) )
        return 0.0;

      // collect all points of the convex hull of the endpoints
      std::vector<Point> endpoints;

      const Tetrahedron tetraA( piA[0], piA[1], piA[2], piA[3] );

      // compute the intersection points
      for( int i = 0; i < 4; ++i )
      {
        const Triangle triangleB ( piB[i], piB[(i+1)%4], piB[(i+2)%4] );

        // test intersection inexactly
        if( CGAL::do_intersect( triangleB, tetraA ) )
        {
          const ExactTriangle triangleEB ( peB[i], peB[(i+1)%4], peB[(i+2)%4] );

          // compute the endpoints of all intersections of all triangles
          for( int j = 0; j < 4; ++j )
          {
            const ExactTriangle triangleEA ( peA[j], peA[(j+1)%4], peA[(j+2)%4] );

            // store the intersection points as inexact points
            const auto result = CGAL::intersection( triangleEA, triangleEB );
            if (result)
            {
              // intersection is PointList
              using PointList = std::vector<ExactPoint>;
              if (const PointList* list = boost::get<PointList>(&*result))
                for ( const auto& p : *list )
                endpoints.push_back( to_inexact( p ) );

              // intersection is triangle
              else if (const ExactTriangle* t = boost::get<ExactTriangle>(&*result))
                for( int k = 0; k < 3; ++k )
                  endpoints.push_back( to_inexact( t->vertex(k) ) );

              // intersection is segment
              else if (const ExactSegment* s = boost::get<ExactSegment>(&*result))
                for( int k = 0; k < 2; ++k )
                  endpoints.push_back( to_inexact( s->vertex(k) ) );

              // intersection is point
              else if (const ExactPoint* p = boost::get<ExactPoint>(&*result))
                endpoints.push_back( to_inexact( *p ) );
            }
          }
        }
      }

      // also add the interior points of tetraA in tetraB
      const Tetrahedron tetraB( piB[0], piB[1], piB[2], piB[3] );
      for( int i = 0; i < 4; ++i )
      {
        if( tetraB.has_on_bounded_side( piA[i] ) )
          endpoints.push_back( piA[i] );

        if( tetraA.has_on_bounded_side( piB[i] ) )
          endpoints.push_back( piB[i] );
      }

      std::sort(endpoints.begin(), endpoints.end());
      auto last = std::unique(endpoints.begin(), endpoints.end(), [](Point& a, Point& b){ return CGAL::squared_distance( a, b ) < 1e-14; });
      endpoints.erase(last, endpoints.end());

      // check if there is an intersection volume
      if( endpoints.size() <= 3 )
        return 0;

      // build convex hull and compute volume
      using Polyhedron = CGAL::Polyhedron_3<Epick>;
      Polyhedron poly;
      CGAL::convex_hull_3(endpoints.begin(), endpoints.end(), poly);
      assert( CGAL::is_closed( poly ) );
      return CGAL::to_double( CGAL::Polygon_mesh_processing::volume( poly ) );
    }

  }; // end of MMeshCachingEntity codim = 0

} // namespace Dune

#endif
