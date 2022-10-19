// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_MMESHLEAFITERATOR_HH
#define DUNE_MMESH_GRID_MMESHLEAFITERATOR_HH


/** \file
 * \brief The MMeshLeafIterator class
 */

// Dune includes
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  /** \brief Iterator over all entities of a given codimension and level of a grid (2D).
   *  \ingroup MMesh
   */

  template<int codim, PartitionIteratorType pitype, class GridImp, typename Enable = void>
  class MMeshLeafIteratorImp {};

  template<int codim, PartitionIteratorType pitype, class GridImp>
  using MMeshLeafIterator = MMeshLeafIteratorImp<codim, pitype, GridImp>;

  /** \brief MMeshLeafIterator for 2D
   */

  template<PartitionIteratorType pitype, class GridImp>
  class MMeshLeafIteratorImp<0, pitype, GridImp, std::enable_if_t<GridImp::dimension == 2>>
  {
  public:
    //! The type of the underlying entities
    using HostGridLeafIterator = typename GridImp::HostGridType::Finite_faces_iterator;
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshLeafIteratorImp() :
      mMesh_(nullptr)
    {}

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(pitype == Interior_Partition
        ? mMesh->partitionHelper().leafInteriorBegin()
        : mMesh->getHostGrid().finite_faces_begin())
    {
      while( proceed() )
        increment();
    }

    /** \brief Constructor which creates the end iterator
     *  \param mMesh    Pointer to grid instance
     *  \param endDummy Here only to distinguish it from the other constructor
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(pitype == Interior_Partition
        ? mMesh->partitionHelper().leafInteriorEnd()
        : mMesh->getHostGrid().finite_faces_end())
    {}

    //! prefix increment
    void increment() {
      do {
        ++hostLeafIterator_;
      } while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      const auto endIterator = (pitype == Interior_Partition ? mMesh_->partitionHelper().leafInteriorEnd() : mMesh_->getHostGrid().finite_faces_end());
      if (hostLeafIterator_ == endIterator)
        return false;
      return !mMesh_->partitionHelper().contains(pitype, hostLeafIterator_->info().partition);
    }

    const GridImp* mMesh_;
    HostGridLeafIterator hostLeafIterator_;
  };

  template<PartitionIteratorType pitype, class GridImp>
  class MMeshLeafIteratorImp<1, pitype, GridImp, std::enable_if_t<GridImp::dimension == 2>>
  {
  private:
    //! The type of the underlying entities
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_edges_iterator;

  public:
    enum {codimension = 1};

    typedef typename GridImp::template Codim<1>::Entity Entity;

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_edges_begin())
    {
      while( proceed() )
        increment();
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_edges_end())
    {}

    //! prefix increment
    void increment() {
      do {
        ++hostLeafIterator_;
      } while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, *hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (hostLeafIterator_ == mMesh_->getHostGrid().finite_edges_end())
        return false;
      return !mMesh_->partitionHelper().contains(pitype, dereference());
    }
    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
  };


  template<PartitionIteratorType pitype, class GridImp>
  class MMeshLeafIteratorImp<2, pitype, GridImp, std::enable_if_t<GridImp::dimension == 2>>
  {
  private:
    //! The type of the underlying entities
    using HostGridLeafIterator = typename GridImp::HostGridType::Finite_vertices_iterator;

  public:
    enum {codimension = GridImp::dimension};

    typedef typename GridImp::template Codim<codimension>::Entity Entity;

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_begin())
    {
      while( proceed() )
        increment();
    }

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_end())
    {}

    //! prefix increment
    void increment() {
      do {
        ++hostLeafIterator_;
      } while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (hostLeafIterator_ == mMesh_->getHostGrid().finite_vertices_end())
        return false;
      return !mMesh_->partitionHelper().contains(pitype, hostLeafIterator_->info().partition);
    }

    const GridImp* mMesh_;
    HostGridLeafIterator hostLeafIterator_;
  };


  /** \brief MMeshLeafIterator for 3D
   *  \ingroup MMesh
   */

  //! Cell iterator
  template<PartitionIteratorType pitype, class GridImp>
  class MMeshLeafIteratorImp<0, pitype, GridImp, std::enable_if_t<GridImp::dimension == 3>>
  {
  public:
    //! The type of the underlying entities
    using HostGridLeafIterator = typename GridImp::HostGridType::Finite_cells_iterator;

    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshLeafIteratorImp() :
      mMesh_(nullptr)
    {}

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(pitype == Interior_Partition
        ? mMesh->partitionHelper().leafInteriorBegin()
        : mMesh->getHostGrid().finite_cells_begin())
    {
      while( proceed() )
        increment();
    }

    /** \brief Constructor which creates the end iterator
     *  \param mMesh    Pointer to grid instance
     *  \param endDummy Here only to distinguish it from the other constructor
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(pitype == Interior_Partition
        ? mMesh->partitionHelper().leafInteriorEnd()
        : mMesh->getHostGrid().finite_cells_end())
    {}

    //! prefix increment
    void increment() {
      do {
        ++hostLeafIterator_;
      } while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      const auto endIterator = (pitype == Interior_Partition ? mMesh_->partitionHelper().leafInteriorEnd() : mMesh_->getHostGrid().finite_cells_end());
      if (hostLeafIterator_ == endIterator)
        return false;
      return !mMesh_->partitionHelper().contains(pitype, hostLeafIterator_->info().partition);
    }

    const GridImp* mMesh_;
    HostGridLeafIterator hostLeafIterator_;
  };


  //! Facet iterator
  template<PartitionIteratorType pitype, class GridImp>
  class MMeshLeafIteratorImp<1, pitype, GridImp, std::enable_if_t<GridImp::dimension == 3>>
  {
  private:
    //! The type of the underlying entities
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_facets_iterator;

  public:
    enum {codimension = 1};

    typedef typename GridImp::template Codim<1>::Entity Entity;

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_facets_begin())
    {
      while( proceed() )
        increment();
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_facets_end())
    {}

    //! prefix increment
    void increment() {
      do {
        ++hostLeafIterator_;
      } while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, *hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (hostLeafIterator_ == mMesh_->getHostGrid().finite_facets_end())
        return false;
      return !mMesh_->partitionHelper().contains(pitype, dereference());
    }

    const GridImp* mMesh_;
    HostGridLeafIterator hostLeafIterator_;
  };


  //! Edge iterator
  template<PartitionIteratorType pitype, class GridImp>
  class MMeshLeafIteratorImp<2, pitype, GridImp, std::enable_if_t<GridImp::dimension == 3>>
  {
  private:
    //! The type of the underlying entities
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_edges_iterator;

  public:
    enum {codimension = 2};

    typedef typename GridImp::template Codim<2>::Entity Entity;

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_edges_begin())
    {
      while( proceed() )
        increment();
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_edges_end())
    {}

    //! prefix increment
    void increment() {
      do {
        ++hostLeafIterator_;
      } while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, *hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (hostLeafIterator_ == mMesh_->getHostGrid().finite_edges_end())
        return false;
      return !mMesh_->partitionHelper().contains(pitype, dereference());
    }

    const GridImp* mMesh_;
    HostGridLeafIterator hostLeafIterator_;
  };


  //! Vertex iterator
  template<PartitionIteratorType pitype, class GridImp>
  class MMeshLeafIteratorImp<3, pitype, GridImp, std::enable_if_t<GridImp::dimension == 3>>
  {
  private:
    //! The type of the underlying entities
    typedef typename GridImp::HostGridType::Vertex_iterator HostGridLeafIterator;

  public:
    enum {codimension = GridImp::dimension};

    typedef typename GridImp::template Codim<codimension>::Entity Entity;

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_begin())
    {
      while( proceed() )
        increment();
    }

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_end())
    {}

    //! prefix increment
    void increment() {
      do {
        ++hostLeafIterator_;
      } while( proceed() );
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, hostLeafIterator_ }};
    }

    //! equality
    bool equals(const MMeshLeafIteratorImp& i) const {
      return hostLeafIterator_ == i.hostLeafIterator_;
    }

  private:
    //! return if this iterator should further be incremented
    bool proceed()
    {
      if (hostLeafIterator_ == HostGridLeafIterator(mMesh_->getHostGrid().finite_vertices_end()))
        return false;
      if (hostLeafIterator_ == mMesh_->getHostGrid().infinite_vertex())
        return true;
      return !mMesh_->partitionHelper().contains(pitype, hostLeafIterator_->info().partition);
    }

    const GridImp* mMesh_;
    HostGridLeafIterator hostLeafIterator_;
  };

}  // namespace Dune

#endif
