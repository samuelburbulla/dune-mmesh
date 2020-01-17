// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_INTERSECTIONITERATOR_HH
#define DUNE_MMESH_INTERFACE_INTERSECTIONITERATOR_HH

/** \file
 * \brief The intersection iterator for the MMesh class
 */

// Dune includes
#include <dune/grid/common/intersection.hh>

// Dune MMesh includes
#include <dune/mmesh/interface/intersections.hh>
#include <dune/mmesh/interface/entity.hh>

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
  class MMeshInterfaceGridLeafIntersectionIterator
  {

    enum {dimension=GridImp::dimension};

    enum {dimensionworld=GridImp::dimensionworld};

    // The type used to store coordinates
    typedef typename GridImp::ctype ctype;

    typedef typename GridImp::template MMeshInterfaceEntity<0> MMeshInterfaceEntity;

  public:
    typedef Dune::Intersection<const GridImp, Dune::MMeshInterfaceGridLeafIntersection<GridImp> > Intersection;

    //! default constructor
    MMeshInterfaceGridLeafIntersectionIterator() {}

    //! constructor for (begin) iterator
    MMeshInterfaceGridLeafIntersectionIterator(const GridImp* grid,
                                  const MMeshInterfaceEntity& hostEntity)
      : grid_(grid)
      , hostEntity_(hostEntity)
      , i_( 0 )
      , nbIdx_( 0 )
    {
      const auto& indexSet = grid_->leafIndexSet();

      for( int d = 0; d < dimension+1; ++d )
      {
        std::array< std::size_t, dimension > ids;
        try {
          for( int i = 0; i < dimension; ++i )
            ids[i] = indexSet.vertexIndexMap().at(
              hostEntity_.first->vertex((hostEntity_.second+i+( (i == 1 && d == 2) ? d+2 : d+1 ))%(dimensionworld+1))->info().index
            );
        } catch (std::exception &e) {
          DUNE_THROW(InvalidStateException, e.what());
        }

        try {
          std::sort(ids.begin(), ids.end());
          maxNbIdx_[d] = std::max( 0, (int)indexSet.indexMap().at( ids ).size() - 2 );
        } catch (std::exception &e) {
          DUNE_THROW(InvalidStateException, e.what());
        }
      }
    }

    //! constructor for end iterator
    MMeshInterfaceGridLeafIntersectionIterator(const GridImp* grid,
                                  const MMeshInterfaceEntity& hostEntity,
                                  bool endDummy)
      : grid_(grid)
      , hostEntity_(hostEntity)
      , i_( dimension + 1 )
      , nbIdx_( 0 )
    {}

    //! returns if iterators reference same intersection
    bool equals(const MMeshInterfaceGridLeafIntersectionIterator& other) const {
      return hostEntity_ == other.hostEntity_
          && i_ == other.i_
          && nbIdx_ == other.nbIdx_;
    }

    //! prefix increment
    void increment()
    {
      if ( nbIdx_ < maxNbIdx_[i_] )
        nbIdx_++;
      else
      {
        ++i_;
        nbIdx_ = 0;
      }
    }

    //! dereferencing
    Intersection dereference() const {
      return MMeshInterfaceGridLeafIntersection<GridImp>(grid_, hostEntity_, i_, nbIdx_);
    }

  private:
    const GridImp* grid_;
    MMeshInterfaceEntity hostEntity_;
    std::size_t i_;
    std::size_t nbIdx_;
    std::array<std::size_t, dimension+1> maxNbIdx_;
  };

}  // namespace Dune

#endif
