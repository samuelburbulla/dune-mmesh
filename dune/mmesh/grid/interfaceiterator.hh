// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_MMESHINTERFACEITERATOR_HH
#define DUNE_GRID_MMESHINTERFACEITERATOR_HH

/** \file
 * \brief The MMeshInterfaceIterator class
 */

// Dune includes
#include <dune/grid/common/gridenums.hh>

namespace Dune
{
  //! Forward declarations
  template<int codim, class GridImp, int dim>
  class MMeshInterfaceIteratorImp;

  //!  The Entity Iterator alias
  template<int codim, class Grid>
  using MMeshInterfaceIterator = EntityIterator<codim, Grid, MMeshInterfaceIteratorImp<codim, Grid, Grid::dimension>>;

  /** \brief Iterator over all interface entities of a given codimension.
   *  \ingroup MMesh
   */

  //! 2D codim 1
  template<class GridImp>
  class MMeshInterfaceIteratorImp<1, GridImp, 2>
  {
  private:
    //! The type of the edge iterator
   using EdgeIterator = typename GridImp::HostGridType::Finite_edges_iterator;

  public:
    typedef typename GridImp::template Codim<1>::Entity Entity;

    explicit MMeshInterfaceIteratorImp(const GridImp* mMesh, bool includeBoundary)
     : mMesh_(mMesh),
       edgeIterator_(mMesh_->getHostGrid().finite_edges_begin()),
       includeBoundary_(includeBoundary)
    {
      while( edgeIterator_ != mMesh_->getHostGrid().finite_edges_end() )
      {
        if (isInterface_())
          break;
        else
          edgeIterator_++;
      }
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshInterfaceIteratorImp(const GridImp* mMesh, bool endDummy, bool includeBoundary)
     : mMesh_(mMesh),
       edgeIterator_(mMesh_->getHostGrid().finite_edges_end()),
       includeBoundary_(includeBoundary)
    {}

    //! prefix increment
    void increment()
    {
      edgeIterator_++;
      while( edgeIterator_ != mMesh_->getHostGrid().finite_edges_end() )
      {
        if (isInterface_())
          break;
        else
          edgeIterator_++;
      }
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, *edgeIterator_ }};
    }

    //! equality
    bool equals(const MMeshInterfaceIteratorImp& iter) const {
      return edgeIterator_ == iter.edgeIterator_;
    }

  private:
    bool isInterface_() const
    {
      if ( includeBoundary_ )
        return mMesh_->isInterface( dereference() );
      else
      {
        return mMesh_->isInterface( dereference() )
          && !mMesh_->getHostGrid().is_infinite(edgeIterator_->first) // exclude boundary segments
          && !mMesh_->getHostGrid().is_infinite(edgeIterator_->first->neighbor(edgeIterator_->second));
      }
    }

    const GridImp* mMesh_;
    EdgeIterator edgeIterator_;
    bool includeBoundary_;
  };


  //! 3D codim 1
  template<class GridImp>
  class MMeshInterfaceIteratorImp<1, GridImp, 3>
  {
  private:
    //! The type of the edge iterator
   using FacetIterator = typename GridImp::HostGridType::Finite_facets_iterator;

  public:
    typedef typename GridImp::template Codim<1>::Entity Entity;

    explicit MMeshInterfaceIteratorImp(const GridImp* mMesh, bool includeBoundary)
     : mMesh_(mMesh),
       facetIterator_(mMesh_->getHostGrid().finite_facets_begin()),
       includeBoundary_(includeBoundary)
    {
      while( facetIterator_ != mMesh_->getHostGrid().finite_facets_end() )
      {
        if (isInterface_())
          break;
        else
          facetIterator_++;
      }
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshInterfaceIteratorImp(const GridImp* mMesh, bool endDummy, bool includeBoundary)
     : mMesh_(mMesh),
       facetIterator_(mMesh_->getHostGrid().finite_facets_end()),
       includeBoundary_(includeBoundary)
    {}

    //! prefix increment
    void increment()
    {
      facetIterator_++;
      while( facetIterator_ != mMesh_->getHostGrid().finite_facets_end() )
      {
        if (isInterface_())
          break;
        else
          facetIterator_++;
      }
    }

    //! dereferencing
    Entity dereference() const {
      return Entity {{ mMesh_, *facetIterator_ }};
    }

    //! equality
    bool equals(const MMeshInterfaceIteratorImp& iter) const {
      return facetIterator_ == iter.facetIterator_;
    }

  private:
    bool isInterface_() const
    {
      if ( includeBoundary_ )
        return mMesh_->isInterface( dereference() );
      else
      {
        return mMesh_->isInterface( dereference() )
          && !mMesh_->getHostGrid().is_infinite(facetIterator_->first) // exclude boundary segments
          && !mMesh_->getHostGrid().is_infinite(facetIterator_->first->neighbor(facetIterator_->second));
      }
    }

    const GridImp* mMesh_;
    FacetIterator facetIterator_;
    bool includeBoundary_;
  };



  //! Forward declarations
  template<class GridImp, int dim>
  class MMeshInterfaceVertexIteratorImp;

  //!  The Entity Iterator alias
  template<class Grid>
  using MMeshInterfaceVertexIterator = EntityIterator<Grid::dimension, Grid, MMeshInterfaceVertexIteratorImp<Grid, Grid::dimension>>;

  /** \brief Iterator over all interface vertices.
   *  \ingroup MMesh
   */

  //! 2D codim 1
  template<class GridImp, int dim>
  class MMeshInterfaceVertexIteratorImp
  {
  private:
    //! The type of the edge iterator
   using InterfaceIterator = MMeshInterfaceIterator<1, GridImp>;
   using InterfaceIteratorImpl = typename InterfaceIterator::Implementation;

  public:
    typedef typename GridImp::template Codim<dim>::Entity Vertex;

    explicit MMeshInterfaceVertexIteratorImp(const GridImp* mMesh, bool includeBoundary)
     : mMesh_(mMesh),
       interfaceIterator_( InterfaceIteratorImpl(mMesh_, includeBoundary) ),
       interfaceEnd_( InterfaceIteratorImpl(mMesh_, true, includeBoundary) ),
       i_(0)
    {
      if( interfaceIterator_ != interfaceEnd_ )
      {
        // add first vertex to visited
        std::size_t id = dereference().impl().hostEntity()->info().id;
        visited.insert( id );
      }
    }

    /** \brief Constructor which creates the end iterator
     *  \param endDummy      Here only to distinguish it from the other constructor
     *  \param mMesh  pointer to grid instance
     */
    explicit MMeshInterfaceVertexIteratorImp(const GridImp* mMesh, bool endDummy, bool includeBoundary)
     : mMesh_(mMesh),
       interfaceIterator_( InterfaceIteratorImpl(mMesh_, true, includeBoundary) ),
       interfaceEnd_( InterfaceIteratorImpl(mMesh_, true, includeBoundary) ),
       i_(0)
    {}

    //! prefix increment
    void increment()
    {
      i_++;

      if ( i_ == dim )
      {
        i_ = 0;
        interfaceIterator_++;
      }

      if ( interfaceIterator_ != interfaceEnd_ )
      {
        // check if already visited
        std::size_t id = dereference().impl().hostEntity()->info().id;
        int count = visited.count( id );
        assert( count <= 1 );
        if ( count > 0 )
          increment();
        else
          visited.insert( id );
      }
    }

    //! dereferencing
    Vertex dereference() const {
      const auto host = interfaceIterator_->impl().hostEntity();
      return Vertex {{ mMesh_, host.first->vertex((host.second+i_+1)%(dim+1)) }};
    }

    //! equality
    bool equals(const MMeshInterfaceVertexIteratorImp& iter) const {
      return interfaceIterator_ == iter.interfaceIterator_ && i_ == iter.i_;
    }

  private:
    const GridImp* mMesh_;
    InterfaceIterator interfaceIterator_, interfaceEnd_;
    std::unordered_set<std::size_t> visited;
    int i_;
  };

}  // namespace Dune

#endif
