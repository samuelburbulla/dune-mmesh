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
    //! The type of the underlying entities
   using HostGridLeafIterator = typename GridImp::HostGridType::Edge_iterator;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshInterfaceGridLeafIteratorImp() :
      mMesh_(nullptr)
    {}

    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_edges_begin()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_edges_end())
    {
      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !mMesh_->isInterface( *hostLeafIterator_ ) )
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

      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !mMesh_->isInterface( *hostLeafIterator_ ) )
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
    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
    HostGridLeafIterator hostLeafIteratorEnd_;
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
      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !hostLeafIterator_->info().isInterface )
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

      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !hostLeafIterator_->info().isInterface )
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
    //! The type of the underlying entities
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_facets_iterator;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshInterfaceGridLeafIteratorImp() :
      mMesh_(nullptr)
    {}

    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_facets_begin()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_facets_end())
    {
      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !mMesh_->isInterface( *hostLeafIterator_ ) )
        ++hostLeafIterator_;
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshInterfaceGridLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_facets_end()),
      hostLeafIteratorEnd_(mMesh->getHostGrid().finite_facets_end())
    {}

    //! prefix increment
    void increment() {
      ++hostLeafIterator_;

      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !mMesh_->isInterface( *hostLeafIterator_ ) )
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
    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
    HostGridLeafIterator hostLeafIteratorEnd_;
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
      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !mMesh_->isInterface( *hostLeafIterator_ ) )
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

      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !mMesh_->isInterface( *hostLeafIterator_ ) )
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
      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !hostLeafIterator_->info().isInterface )
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

      while( hostLeafIterator_ != hostLeafIteratorEnd_ && !hostLeafIterator_->info().isInterface )
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
    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
    HostGridLeafIterator hostLeafIteratorEnd_;
  };

}  // namespace Dune

#endif
