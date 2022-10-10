// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_MMESHLEAFITERATOR_HH
#define DUNE_MMESH_INTERFACE_MMESHLEAFITERATOR_HH


/** \file
 * \brief The MMeshInterfaceGridLeafIterator class
 */

// Dune includes
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  /** \brief Iterator over all entities of a given codimension and level of a grid (2D).
   *  \ingroup MMesh
   */

  template<int codim, PartitionIteratorType pitype, class GridImp, typename Enable = void>
  class MMeshInterfaceGridLeafIteratorImp {};

  template<int codim, PartitionIteratorType pitype, class GridImp>
  using MMeshInterfaceGridLeafIterator = MMeshInterfaceGridLeafIteratorImp<codim, pitype, GridImp>;

  /** \brief MMeshInterfaceGridLeafIteratorImp for 2D
   */

  template<PartitionIteratorType pitype, class GridImp>
  class MMeshInterfaceGridLeafIteratorImp<0, pitype, GridImp, std::enable_if_t<GridImp::dimensionworld == 2>>
  {
  private:
    //! The type of the underlying entity iterator
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_faces_iterator;
    //! The type of the underlying interface host entity
   using HostGridFacet = typename GridImp::MMeshType::FacetHandle;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshInterfaceGridLeafIteratorImp() :
      mMesh_(nullptr)
    {}

    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(pitype == Interior_Partition
        ? mMesh->partitionHelper().leafInteriorBegin()
        : mMesh->getHostGrid().finite_faces_begin()),
      face_(0)
    {
      if( proceed() )
        increment();
    }

    /** \brief Constructor which creates the end iterator
     *  \param mMesh    Pointer to grid instance
     *  \param endDummy Here only to distinguish it from the other constructor
     */
    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(pitype == Interior_Partition
        ? mMesh->partitionHelper().leafInteriorEnd()
        : mMesh->getHostGrid().finite_faces_end()),
      face_(0)
    {}

    //! prefix increment
    void increment() {
      do {
        face_++;
        if (face_ == 3)
        {
          ++hostLeafIterator_;
          face_ = 0;
        }
      }
      while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, HostGridFacet( hostLeafIterator_, face_ ) }};
    }

    //! equality
    bool equals(const MMeshInterfaceGridLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_ && face_ == i.face_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      const auto endIterator = (pitype == Interior_Partition ? mMesh_->partitionHelper().leafInteriorEnd() : mMesh_->getHostGrid().finite_faces_end());
      if (hostLeafIterator_ == endIterator)
        return false;

      HostGridFacet facet ( hostLeafIterator_, face_ );
      if (!mMesh_->isInterface( facet ))
        return true;

      const auto mirrored = hostLeafIterator_->neighbor( face_ );
      if ( hostLeafIterator_->info().insertionIndex > mirrored->info().insertionIndex )
        return true;

      return !mMesh_->partitionHelper().contains(pitype, dereference());
    }

    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
    int face_;
  };


  template<PartitionIteratorType pitype, class GridImp>
  class MMeshInterfaceGridLeafIteratorImp<1, pitype, GridImp, std::enable_if_t<GridImp::dimensionworld == 2>>
  {
  private:
    //! The type of the underlying entities
    typedef typename GridImp::HostGridType::Vertex_iterator HostGridLeafIterator;

  public:
    enum {codimension = GridImp::dimension};

    typedef typename GridImp::template Codim<codimension>::Entity Entity;

    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_begin()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_vertices_end())
    {
      while( proceed() )
        ++hostLeafIterator_;
    }

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_end()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_vertices_end())
    {}

    //! prefix increment
    void increment() {
      ++hostLeafIterator_;

      while( proceed() )
        ++hostLeafIterator_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshInterfaceGridLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (hostLeafIterator_ == hostLeafIteratorEnd_)
        return false;
      if (!hostLeafIterator_->info().isInterface)
        return true;
      return !mMesh_->partitionHelper().contains(pitype, dereference());
    }

    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
    HostGridLeafIterator hostLeafIteratorEnd_;
  };


  /** \brief MMeshInterfaceGridLeafIteratorImp for 3D
   *  \ingroup MMesh
   */

  template<PartitionIteratorType pitype, class GridImp>
  class MMeshInterfaceGridLeafIteratorImp<0, pitype, GridImp, std::enable_if_t<GridImp::dimensionworld == 3>>
  {
  private:
    //! The type of the underlying entity iterator
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_cells_iterator;
    //! The type of the underlying interface host entity
   using HostGridFacet = typename GridImp::MMeshType::FacetHandle;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshInterfaceGridLeafIteratorImp() :
      mMesh_(nullptr)
    {}

    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(pitype == Interior_Partition
        ? mMesh->partitionHelper().leafInteriorBegin()
        : mMesh->getHostGrid().finite_cells_begin()),
      face_(0)
    {
      if( proceed() )
        increment();
    }

    /** \brief Constructor which creates the end iterator
     *  \param mMesh    Pointer to grid instance
     *  \param endDummy Here only to distinguish it from the other constructor
     */
    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(pitype == Interior_Partition
        ? mMesh->partitionHelper().leafInteriorEnd()
        : mMesh->getHostGrid().finite_cells_end()),
      face_(0)
    {}

    //! prefix increment
    void increment() {
      do {
        face_++;
        if (face_ == 4)
        {
          ++hostLeafIterator_;
          face_ = 0;
        }
      }
      while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, HostGridFacet( hostLeafIterator_, face_ ) }};
    }

    //! equality
    bool equals(const MMeshInterfaceGridLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_ && face_ == i.face_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      const auto endIterator = (pitype == Interior_Partition ? mMesh_->partitionHelper().leafInteriorEnd() : mMesh_->getHostGrid().finite_cells_end());
      if (hostLeafIterator_ == endIterator)
        return false;

      HostGridFacet facet ( hostLeafIterator_, face_ );
      if (!mMesh_->isInterface( facet ))
        return true;

      const auto mirrored = hostLeafIterator_->neighbor( face_ );
      if ( hostLeafIterator_->info().insertionIndex > mirrored->info().insertionIndex )
        return true;

      return !mMesh_->partitionHelper().contains(pitype, dereference());
    }

    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
    int face_;
  };

  template<PartitionIteratorType pitype, class GridImp>
  class MMeshInterfaceGridLeafIteratorImp<1, pitype, GridImp, std::enable_if_t<GridImp::dimensionworld == 3>>
  {
  private:
    //! The type of the underlying entities
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_edges_iterator;

  public:
    enum {codimension = 1};

    typedef typename GridImp::template Codim<1>::Entity Entity;

    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_edges_begin()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_edges_end())
    {
      while( proceed() )
        ++hostLeafIterator_;
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_edges_end()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_edges_end())
    {}

    //! prefix increment
    void increment() {
      ++hostLeafIterator_;

      while( proceed() )
        ++hostLeafIterator_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, *hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshInterfaceGridLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (hostLeafIterator_ == hostLeafIteratorEnd_)
        return false;
      if (!mMesh_->isInterface( *hostLeafIterator_ ))
        return true;
      return !mMesh_->partitionHelper().contains(pitype, dereference());
    }

    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
    HostGridLeafIterator hostLeafIteratorEnd_;
  };

  template<PartitionIteratorType pitype, class GridImp>
  class MMeshInterfaceGridLeafIteratorImp<2, pitype, GridImp, std::enable_if_t<GridImp::dimensionworld == 3>>
  {
  private:
    //! The type of the underlying entities
    typedef typename GridImp::HostGridType::Vertex_iterator HostGridLeafIterator;

  public:
    enum {codimension = GridImp::dimension};

    typedef typename GridImp::template Codim<codimension>::Entity Entity;

    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_begin()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_vertices_end())
    {
      while( proceed() )
        ++hostLeafIterator_;
    }

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_end()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_vertices_end())
    {}

    //! prefix increment
    void increment() {
      ++hostLeafIterator_;

      while( proceed() )
        ++hostLeafIterator_;
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshInterfaceGridLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (hostLeafIterator_ == hostLeafIteratorEnd_)
        return false;
      if (!hostLeafIterator_->info().isInterface)
        return true;
      return !mMesh_->partitionHelper().contains(pitype, dereference());
    }

    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
    HostGridLeafIterator hostLeafIteratorEnd_;
  };

}  // namespace Dune

#endif
