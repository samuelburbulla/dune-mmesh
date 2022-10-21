// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_INTERSECTIONITERATOR_HH
#define DUNE_MMESH_GRID_INTERSECTIONITERATOR_HH

/** \file
 * \brief The intersection iterator for the MMesh class
 */

// Dune includes
#include <dune/grid/common/intersection.hh>

// Dune MMesh includes
#include <dune/mmesh/grid/intersections.hh>
#include <dune/mmesh/grid/entity.hh>

namespace Dune
{

  /** \brief Iterator over all element neighbors
   * \ingroup MMesh
   * Mesh entities of codimension 0 ("elements") allow to visit all neighbors, where
   * a neighbor is an entity of codimension 0 which has a common entity of codimension 1
   * These neighbors are accessed via a IntersectionIterator. This allows the implement
   * non-matching meshes. The number of neighbors may be different from the number
   * of an element!
   */
  template<class GridImp>
  class MMeshLeafIntersectionIterator
  {

    enum {dim=GridImp::dimension};

    enum {dimworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::template HostGridEntity<0> HostGridEntity;

  public:
    typedef Dune::Intersection<const GridImp, Dune::MMeshLeafIntersection<GridImp> > Intersection;

    //! default constructor
    MMeshLeafIntersectionIterator() : i_( 0 ) {}

    //! constructor for (begin) iterator
    MMeshLeafIntersectionIterator(const GridImp* mMesh,
                                  const HostGridEntity& hostEntity)
      : mMesh_(mMesh)
      , hostEntity_(hostEntity)
      , i_( 0 )
    {
      while( proceed() )
        increment();
    }

    //! constructor for end iterator
    MMeshLeafIntersectionIterator(const GridImp* mMesh,
                                  const HostGridEntity& hostEntity,
                                  bool endDummy)
      : mMesh_(mMesh)
      , hostEntity_(hostEntity)
      , i_( dim + 1 )
    {}

    //! returns if iterators reference same intersection
    bool equals(const MMeshLeafIntersectionIterator& other) const {
      return i_ == other.i_ && hostEntity_ == other.hostEntity_;
    }

    //! prefix increment
    void increment() {
      do {
        ++i_;
      } while( proceed() );
    }

    //! \brief dereferencing
    Intersection dereference() const {
      // remark: the i-th intersection in CGAL corresponds to the (dim-i)-th edge in DUNE
      return MMeshLeafIntersection<GridImp>(mMesh_, hostEntity_, dim-i_);
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (i_ == dim + 1)
        return false;

      // on ghosts return intersection only if adjacent entity is interior
      if (hostEntity_->info().partition == 2)
      {
        const auto& neighborHostEntity = hostEntity_->neighbor(dim-i_);

        if (mMesh_->getHostGrid().is_infinite(neighborHostEntity))
          return true;

        // skip ghost neighbors
        if (neighborHostEntity->info().partition != 0)
          return true;
      }
      return false;
    }

    const GridImp* mMesh_;
    HostGridEntity hostEntity_;
    int i_;
  };

}  // namespace Dune

#endif
