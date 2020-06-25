#ifndef DUNE_MMESH_MISC_TWISTUTILITY_HH
#define DUNE_MMESH_MISC_TWISTUTILITY_HH

#if HAVE_DUNE_FEM

namespace Dune
{

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
        return intersection.indexInInside() % dim;
      }

      //! \brief return twist for outer face
      static inline int twistInNeighbor(const GridType& grid, const LeafIntersection& intersection )
      {
        // TODO: 3d
        return 1 - intersection.indexInOutside() % 2;
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
        return intersection.indexInInside() % dim;
      }

      //! \brief return twist for outer face
      static inline int twistInNeighbor(const GridType& grid, const LeafIntersection& intersection )
      {
        // TODO: 3d
        return 1 - intersection.indexInOutside() % 2;
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

#endif // #if HAVE_DUNE_FEM

#endif // #ifndef DUNE_MMESH_MISC_TWISTUTILITY_HH
