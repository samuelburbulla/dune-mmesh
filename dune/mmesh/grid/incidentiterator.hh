// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_MMESHINCIDENTITERATOR_HH
#define DUNE_GRID_MMESHINCIDENTITERATOR_HH

// CGAL includes
#include <CGAL/circulator.h>

/** \file
 * \brief The MMeshIncidentIterator class
 */

// Dune includes
#include <dune/grid/common/gridenums.hh>

namespace Dune
{
  //! Forward declarations
  template<class GridImp, int dim>
  class MMeshIncidentIteratorImp;

  //!  The Entity Iterator alias
  template<class Grid>
  using MMeshIncidentIterator = EntityIterator<0, Grid, MMeshIncidentIteratorImp<Grid, Grid::dimension>>;

  /** \brief Iterator over all incident entities of a given codimension.
   *  \ingroup MMesh
   */

  //! 2D
  template<class GridImp>
  class MMeshIncidentIteratorImp<GridImp, 2>
  {
  private:
    //! The type of the requested vertex entity
    typedef typename GridImp::template HostGridEntity<GridImp::dimension> HostGridVertex;
    typedef typename GridImp::template HostGridEntity<0> HostGridEntity;

    //! The type of the element circulator
   using Circulator = typename GridImp::Element_circulator;
   using ElementContainer = std::vector<HostGridEntity>;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshIncidentIteratorImp(const GridImp* mMesh, const HostGridVertex& hostEntity)
    : mMesh_(mMesh),
      i_(0)
    {
      Circulator circulator = mMesh->getHostGrid().incident_faces(hostEntity);
      for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
        if (!mMesh->getHostGrid().is_infinite(circulator))
          elementContainer_.push_back( circulator );
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshIncidentIteratorImp(const GridImp* mMesh, const HostGridVertex& hostEntity, bool endDummy) :
      mMesh_(mMesh),
      i_(0)
    {
      Circulator circulator = mMesh->getHostGrid().incident_faces(hostEntity);
      for ( std::size_t i = 0; i < CGAL::circulator_size(circulator); ++i, ++circulator )
        if (!mMesh->getHostGrid().is_infinite(circulator))
          ++i_;
    }

    //! prefix increment
    void increment() {
      ++i_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, elementContainer_[i_] }};
    }

    //! equality
    bool equals(const MMeshIncidentIteratorImp& iter) const {
      return i_ == iter.i_;
    }

  private:
    const GridImp* mMesh_;
    ElementContainer elementContainer_;
    std::size_t i_;
  };


  //! 3D
  template<class GridImp>
  class MMeshIncidentIteratorImp<GridImp, 3>
  {
  private:
    //! The type of the requested vertex entity
    typedef typename GridImp::template HostGridEntity<GridImp::dimension> HostGridVertex;
    typedef typename GridImp::template HostGridEntity<0> HostGridEntity;

    //! The type of the element output
    using ElementOutput = typename GridImp::ElementOutput;

    //! The type of the element container
    using ElementContainer = std::vector<HostGridEntity>;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshIncidentIteratorImp(const GridImp* mMesh, const HostGridVertex& hostEntity)
    : mMesh_(mMesh),
      i_(0)
    {
      ElementOutput elements;
      mMesh_->getHostGrid().finite_incident_cells( hostEntity, std::back_inserter(elements) );

      typename ElementOutput::iterator fit;
      for(fit = elements.begin(); fit != elements.end(); fit++)
        elementContainer_.push_back( *fit );
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshIncidentIteratorImp(const GridImp* mMesh, const HostGridVertex& hostEntity, bool endDummy) :
      mMesh_(mMesh),
      i_(0)
    {
      ElementOutput elements;
      mMesh_->getHostGrid().finite_incident_cells( hostEntity, std::back_inserter(elements) );

      typename ElementOutput::iterator fit;
      for(fit = elements.begin(); fit != elements.end(); fit++)
          ++i_;
    }

    //! prefix increment
    void increment() {
      ++i_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, elementContainer_[i_] }};
    }

    //! equality
    bool equals(const MMeshIncidentIteratorImp& iter) const {
      return i_ == iter.i_;
    }

  private:
    const GridImp* mMesh_;
    ElementContainer elementContainer_;
    std::size_t i_;
  };

}  // namespace Dune

#endif
