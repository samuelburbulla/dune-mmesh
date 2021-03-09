// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_CONNECTEDCOMPONENT_HH
#define DUNE_MMESH_INTERFACE_CONNECTEDCOMPONENT_HH

/** \file
 * \brief The MMeshInterfaceConnectedComponent class
 */

#include <set>

// Dune includes
#include <dune/grid/common/grid.hh>

// MMesh includes
#include <dune/mmesh/interface/geometry.hh>
#include <dune/mmesh/interface/cachingentity.hh>

namespace Dune
{

  //**********************************************************************
  //
  // --MMeshInterfaceConnectedComponent
  //
  /** \brief The implementation of connected components in a MMeshInterfaceGrid
   *  \ingroup MMeshInterfaceGrid
   *  The connected component copies the vertex coordinates and ids.
   *
   */
  template<class GridImp>
  class MMeshInterfaceConnectedComponent
  {
    template <class GridImp_>
    friend class MMeshInterfaceGridLeafIndexSet;

    template <class GridImp_>
    friend class MMeshInterfaceGridLocalIdSet;

    template <class GridImp_>
    friend class MMeshInterfaceGridGlobalIdSet;

  private:
    // this type
    typedef MMeshInterfaceConnectedComponent<GridImp> ThisType;

    // type of scalars
    typedef typename GridImp::ctype ctype;

    // type of the host grid
    typedef typename GridImp::HostGridType HostGrid;

    // equivalent interface entity in the host grid as pointer
    typedef typename GridImp::template Codim<0>::Entity Element;

    // vertex in the host grid as object
    typedef typename GridImp::HostGridType::Vertex HostGridVertex;

    // type of caching entity
    using CachingEntity = MMeshInterfaceCachingEntity< 0, GridImp::dimension, const GridImp >;

    // id type
    using IdType = MMeshImpl::MultiId;

    // vertex storage
    using Vertices = std::array<HostGridVertex, GridImp::dimension+1>;

  public:
    typedef MMeshInterfaceGridGeometry<GridImp::dimension, GridImp::dimension+1, GridImp> Geometry;

    MMeshInterfaceConnectedComponent() {};

    //! Construct connected component with a single element
    explicit MMeshInterfaceConnectedComponent(const Element& element)
    {
      children_.emplace_back( element );
    }

    //! Add element to connected component
    void add (const Element& element)
    {
      children_.emplace_back( element );
      assert( children_.size() <= 2 ); // at the moment, more children are not supported
    }

    //! Return list of caching entities in this component
    const std::vector<CachingEntity>& children() const
    {
      return children_;
    }

    //! Return number of caching entities in this component
    const std::size_t size() const
    {
      return children_.size();
    }

  private:
    //! list of caching entities
    std::vector< CachingEntity > children_;

  };

} // namespace Dune

#endif
