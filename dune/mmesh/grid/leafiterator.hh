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
  private:
    //! The type of the underlying entities
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_faces_iterator;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshLeafIteratorImp() :
      mMesh_(nullptr)
    {}

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_faces_begin())
    {}

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_faces_end())
    {}

    //! prefix increment
    void increment() {
      ++hostLeafIterator_;
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
    {}

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
      ++hostLeafIterator_;
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
    {}

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
      ++hostLeafIterator_;
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
  private:
    //! The type of the underlying entities
   using HostGridLeafIterator = typename GridImp::HostGridType::Finite_cells_iterator;

  public:
    enum {codimension = 0};

    typedef typename GridImp::template Codim<0>::Entity Entity;

    explicit MMeshLeafIteratorImp() :
      mMesh_(nullptr)
    {}

    explicit MMeshLeafIteratorImp(const GridImp* mMesh) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_cells_begin())
    {}

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_cells_end())
    {}

    //! prefix increment
    void increment() {
      ++hostLeafIterator_;
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
    {}

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
      ++hostLeafIterator_;
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
    {}

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
      ++hostLeafIterator_;
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
    {}

    /** \brief Constructor which create the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshLeafIteratorImp(const GridImp* mMesh, bool endDummy) :
      mMesh_(mMesh),
      hostLeafIterator_(mMesh->getHostGrid().finite_vertices_end())
    {
      hostLeafIterator_--; // a bit surprisingly, we have to do this to match the number of vertices
    }

    //! prefix increment
    void increment() {
      ++hostLeafIterator_;
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
    const GridImp* mMesh_;

    HostGridLeafIterator hostLeafIterator_;
  };

}  // namespace Dune

#endif
