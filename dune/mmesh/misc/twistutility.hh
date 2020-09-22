#ifndef DUNE_MMESH_MISC_TWISTUTILITY_HH
#define DUNE_MMESH_MISC_TWISTUTILITY_HH

#include <dune/grid/test/checktwists.hh>

namespace Dune
{

  template< class Grid >
  class MMeshInterfaceGrid;

  //! Helper functions to compute twist brute-force
  namespace _MMeshTwistImpl
  {
    template< class Intersection >
    int computeTwist(const Intersection& intersection, const bool outside = true)
    {
      static constexpr int dimension = Intersection::dimensionworld;

      typedef typename Intersection::Entity Entity;
      typedef typename Entity::Geometry Geometry;
      typedef typename Geometry::ctype ctype;

      const ctype tolerance = std::numeric_limits< ctype >::epsilon();

      const Entity entityIn = intersection.inside();
      const Entity entityOut = intersection.outside();

      auto ref = referenceElement( entityIn.geometry() );

      const int n = outside ? intersection.indexInOutside() : intersection.indexInInside();
      const auto geo = outside ? entityOut.geometry() : entityIn.geometry();
      const auto lgeo = outside ? intersection.geometryInOutside() : intersection.geometryInInside();

      const auto x0 = geo.corner( ref.subEntity( n, 1, 0, dimension ) );
      const auto x1 = geo.corner( ref.subEntity( n, 1, 1, dimension ) );

      int numCorners = dimension;
      for( int i = 0; i < numCorners; ++i )
      {
        if( (geo.global( lgeo.corner(i) ) - x0).two_norm() <= tolerance )
        {
          if( (geo.global( lgeo.corner((i+1)%numCorners) ) - x1).two_norm() <= tolerance )
            return i;
          else
            return -i-1;
        }
      }

      DUNE_THROW(InvalidStateException, "Twist not found.");
      return 0;
    }
  }

  namespace MMeshTwist
  {
    template< class LeafIntersection >
    static inline int twistInSelf(const LeafIntersection& intersection)
    {
      static constexpr int dim = LeafIntersection::dimensionworld;
      if constexpr (dim == 2)
        return intersection.indexInInside() % dim;
      else
        return 0;
    }

    template< class LeafIntersection >
    static inline int twistInNeighbor(const LeafIntersection& intersection )
    {
      static constexpr int dim = LeafIntersection::dimensionworld;
      if constexpr (dim == 2)
        return 1 - intersection.indexInOutside() % 2;
      else
        return _MMeshTwistImpl::computeTwist(intersection, true);
    }
  }

  namespace MMeshInterfaceTwist
  {
    template< class LeafIntersection >
    static inline int twistInSelf(const LeafIntersection& intersection)
    {
      static constexpr int dim = LeafIntersection::dimensionworld;
      if constexpr (dim == 2)
        return 0;
      else
        return _MMeshTwistImpl::computeTwist(intersection, false);
    }

    template< class LeafIntersection >
    static inline int twistInNeighbor(const LeafIntersection& intersection )
    {
      static constexpr int dim = LeafIntersection::dimensionworld;
      if constexpr (dim == 2)
        return 0;
      else
        return _MMeshTwistImpl::computeTwist(intersection, true);
    }
  }

  namespace Fem
  {

    template< class Grid >
    struct TwistUtility;

    /** \brief Specialization of TwistUtility for MMesh.
    */
    template< class HostGrid, int dim >
    struct TwistUtility< MMesh< HostGrid, dim > >
    {
      typedef MMesh< HostGrid, dim > GridType;
      typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
      typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
      typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    public:
      //! \brief return twist for inner face
      static inline int twistInSelf(const GridType& grid, const LeafIntersection& intersection)
      {
        return MMeshTwist::twistInSelf(intersection);
      }

      //! \brief return twist for outer face
      static inline int twistInNeighbor(const GridType& grid, const LeafIntersection& intersection )
      {
        return MMeshTwist::twistInNeighbor(intersection);
      }

      //! \brief return element geometry type of inside or outside entity
      template <class Intersection>
      static inline GeometryType
      elementGeometry(const Intersection& intersection, const bool inside)
      {
        return GeometryType( Dune::GeometryType::simplex, dim );
      }

    private:
      TwistUtility(const TwistUtility&);
      TwistUtility& operator=(const TwistUtility&);
    };


    /** \brief Specialization of TwistUtility for MMesh InterfaceGrid.
    */
    template< class MMesh >
    struct TwistUtility< MMeshInterfaceGrid<MMesh> >
    {
      typedef MMeshInterfaceGrid<MMesh> GridType;
      static constexpr int dim = GridType::dimension;
      typedef typename GridType::Traits::LeafIntersectionIterator LeafIntersectionIterator;
      typedef typename LeafIntersectionIterator::Intersection LeafIntersection;
      typedef typename GridType::Traits::LevelIntersectionIterator LevelIntersectionIterator;
      typedef typename LevelIntersectionIterator::Intersection LevelIntersection;

    public:
      //! \brief return twist for inner face
      static inline int twistInSelf(const GridType& grid, const LeafIntersection& intersection)
      {
        return MMeshInterfaceTwist::twistInSelf(intersection);
      }

      //! \brief return twist for outer face
      static inline int twistInNeighbor(const GridType& grid, const LeafIntersection& intersection )
      {
        return MMeshInterfaceTwist::twistInNeighbor(intersection);
      }

      //! \brief return element geometry type of inside or outside entity
      template <class Intersection>
      static inline GeometryType
      elementGeometry(const Intersection& intersection, const bool inside)
      {
        return GeometryType( Dune::GeometryType::simplex, dim );
      }

    private:
      TwistUtility(const TwistUtility&);
      TwistUtility& operator=(const TwistUtility&);
    };

  }  // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_MMESH_MISC_TWISTUTILITY_HH
