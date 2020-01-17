// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_CONNECTEDCOMPONENT_HH
#define DUNE_MMESH_INTERFACE_CONNECTEDCOMPONENT_HH

/** \file
 * \brief The MMeshInterfaceFatherEntity class
 */

#include <set>

// Dune includes
#include <dune/grid/common/grid.hh>

// MMesh includes
#include "geometry.hh"

namespace Dune
{

  //**********************************************************************
  //
  // --MMeshInterfaceConnectedComponent
  // --Entity
  //
  /** \brief The implementation of connected components in a MMeshInterfaceGrid
   *   \ingroup MMeshInterfaceGrid
   *  The connected component copies the vertex coordinates.
   *
   */
  template<int dim, class GridImp>
  class MMeshInterfaceConnectedComponent<0,dim,GridImp> :
    public EntityDefaultImplementation <0,dim,GridImp,MMeshInterfaceConnectedComponent>
  {
    template <class GridImp_>
    friend class MMeshInterfaceGridLeafIndexSet;

    template <class GridImp_>
    friend class MMeshInterfaceGridLocalIdSet;

    template <class GridImp_>
    friend class MMeshInterfaceGridGlobalIdSet;

  private:
    // this type
    typedef MMeshInterfaceConnectedComponent<0,dim,GridImp> ThisType;

    // type of scalars
    typedef typename GridImp::ctype ctype;

    // type of the host grid
    typedef typename GridImp::HostGridType HostGrid;

    // equivalent interface entity in the host grid as pointer
    typedef typename GridImp::template Codim<0>::Entity Element;

    // vertex in the host grid as object
    typedef typename GridImp::HostGridType::Vertex HostGridVertex;

    // id type
    using IdType = Impl::MultiId;

  public:
    typedef MMeshInterfaceGridGeometry<dim, dim+1, GridImp> Geometry;

    MMeshInterfaceConnectedComponent() {};

    explicit MMeshInterfaceConnectedComponent(const Element& element)
    {
      // copy vertices
      const auto& host = element.impl().hostEntity();
      for ( int i = 0; i < dim+1; ++i )
        vertices_[i] = *(host.first->vertex( (host.second+i+1)%(dim+2) ));

      // store id
      id_ = element.impl().grid().globalIdSet().id( element );
    }

    //! returns true if host entities are equal
    bool equals(const MMeshInterfaceConnectedComponent& other) const
    {
      return id_ == other.id_;
    }

    //! returns true if host entities are equal
    bool operator==(const MMeshInterfaceConnectedComponent& other) const
    {
      return this->equals(other);
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return false;
    }

    //! returns true if this entity is new after adaptation
    const bool isNew () const
    {
      return false;
    }

    //! returns true if this entity will vanish after adaptation
    const bool mightVanish () const
    {
      return true;
    }

    //! Level of this element
    int level () const
    {
      // we only have one level
      return 0;
    }

    //! The partition type for parallel computing
    PartitionType partitionType () const {
      return PartitionType::InteriorEntity; /* dummy */
    }

    //! Geometry of this entity
    Geometry geometry () const
    {
      return Geometry( vertices_ );
    }

    //! Return vertex for given index
    const HostGridVertex& vertex ( int i ) const
    {
      assert( i < dim+1 );
      return vertices_[i];
    }

    //! Return the id of this entity
    IdType id() const
    {
      return id_;
    }

    //! Return the children of this connected component
    std::vector<ThisType> children() const
    {
      return { *this }; // TODO: this is different for coarsening
    }

    //! Return the intersection volume of this with given entity
    template< class Entity >
    ctype intersectionVolume( const Entity& entity ) const
    {
      return entity.geometry().volume(); // TODO: this is different for coarsening
    }

    //! Return the number of subEntities of codimension cc
    std::size_t subEntities (std::size_t cc) const
    {
      // we have a simplex grid
      int n = dim+1;
      int k = dim-cc+1;

      // binomial: n over k
      int binomial=1;
      for (int i=n-k+1; i<=n; i++)
        binomial *= i;
      for (long i=2; i<=k; i++)
        binomial /= i;

      return binomial;
    }

  private:
    //! list of entity vertices
    std::array<HostGridVertex, dim+1> vertices_;
    IdType id_;

  }; // end of MMeshInterfaceConnectedComponent codim = 0

} // namespace Dune

#endif
