// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_HIERARCHICITERATOR_HH
#define DUNE_MMESH_INTERFACE_HIERARCHICITERATOR_HH

/** \file
 * \brief The MMeshInterfaceGridHierarchicIterator class
 */

namespace Dune
{

  //**********************************************************************
  //
  /** \brief Iterator over the descendants of an entity.
   * \ingroup MMesh
     Mesh entities of codimension 0 ("elements") allow to visit all entities of
     codimension 0 obtained through nested, hierarchic refinement of the entity.
     Iteration over this set of entities is provided by the HierarchicIterator,
     starting from a given entity.
   */
  template<class GridImp>
  class MMeshInterfaceGridHierarchicIterator
  {
  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    //! the default constructor of the iterator
    explicit MMeshInterfaceGridHierarchicIterator(const GridImp* mMesh, const Entity& startEntity, int maxLevel) :
      mMesh_(mMesh),
      i_(0),
      startEntity_(startEntity)
    {}

    //! the constructor of the end iterator
    explicit MMeshInterfaceGridHierarchicIterator(const GridImp* mMesh, const Entity& startEntity, int maxLevel, bool endDummy) :
      mMesh_(mMesh),
      i_(1),
      startEntity_(startEntity)
    {}

    //! increment iterator
    void increment()
    {
      ++i_;
    }

    //! dereference iterator
    Entity dereference() const {
      return startEntity_;
    }

    //! compare iterators
    bool equals(const MMeshInterfaceGridHierarchicIterator& other) const {
      return startEntity_ == other.startEntity_ && i_ == other.i_;
    }

  private:
    const GridImp* mMesh_;
    int i_;
    Entity startEntity_;

  };

} // end namespace Dune

#endif
