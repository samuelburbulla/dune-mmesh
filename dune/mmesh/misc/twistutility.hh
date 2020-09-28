#ifndef DUNE_MMESH_MISC_TWISTUTILITY_HH
#define DUNE_MMESH_MISC_TWISTUTILITY_HH

namespace Dune
{

  template< class Grid >
  class MMeshInterfaceGrid;

  namespace MMeshTwist
  {
    template< class LeafIntersection >
    static inline int twistInSelf(const LeafIntersection& intersection)
    {
      return 0;
    }

    template< class LeafIntersection >
    static inline int twistInNeighbor(const LeafIntersection& intersection )
    {
      return 0;
    }
  }

  namespace MMeshInterfaceTwist
  {
    template< class LeafIntersection >
    static inline int twistInSelf(const LeafIntersection& intersection)
    {
      return 0;
    }

    template< class LeafIntersection >
    static inline int twistInNeighbor(const LeafIntersection& intersection )
    {
      return 0;
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
