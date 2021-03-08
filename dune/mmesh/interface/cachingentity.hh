// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_CACHINGENTITY_HH
#define DUNE_MMESH_INTERFACE_CACHINGENTITY_HH

/** \file
 * \brief The MMeshInterfaceCachingEntity class
 */

#include <set>

// Dune includes
#include <dune/grid/common/grid.hh>

// MMesh includes
#include <dune/mmesh/grid/multiid.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>

namespace Dune
{
  // Forward declarations
  template<int codim, int dim, class GridImp>
  class MMeshInterfaceCachingEntity;

  //**********************************************************************
  //
  // --MMeshInterfaceCachingEntity
  // --Entity
  //
  /** \brief The implementation of caching entities in a MMesh interface grid
   *  \ingroup MMeshInterfaceGrid
   *
   */
  template<int dim, class GridImp>
  class MMeshInterfaceCachingEntity<0,dim,GridImp> :
    public MMeshInterfaceGridEntity<0,dim,GridImp>
  {
    template <class GridImp_>
    friend class MMeshInterfaceGridLeafIndexSet;

    template <class GridImp_>
    friend class MMeshInterfaceGridLocalIdSet;

    template <class GridImp_>
    friend class MMeshInterfaceGridGlobalIdSet;

  private:
    // this type
    typedef MMeshInterfaceCachingEntity<0,dim,GridImp> ThisType;

    // base type
    typedef MMeshInterfaceGridEntity<0,dim,GridImp> BaseType;

    // type of scalars
    typedef typename GridImp::ctype ctype;

    // type of the host grid
    typedef typename GridImp::HostGridType HostGrid;

    // type of entity
    typedef typename GridImp::template Codim<0>::Entity Element;

    // standard MMesh entity implementation
    typedef typename GridImp::template Codim<0>::Entity MMeshEntity;

    // type of ids
    typedef MMeshImpl::MultiId IdType;

  public:
    // geometry type
    typedef AffineGeometry<ctype, dim, dim+1> Geometry;

    // type of global coordinate
    typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

    MMeshInterfaceCachingEntity() = delete;

    explicit MMeshInterfaceCachingEntity(const Element& element)
      : BaseType(
        &element.impl().grid(),
        element.impl().hostEntity(),
        element.impl().grid().globalIdSet().id( element )
      )
    {
      for ( int i = 0; i < dim+1; ++i )
        this->vertex_[i] = element.geometry().corner(i);
    }

    //! returns true if host entities are equal
    bool equals(const MMeshInterfaceCachingEntity& other) const
    {
      return this->id_ == other.id_;
    }

    //! returns true if host entities are equal
    bool operator==(const MMeshInterfaceCachingEntity& other) const
    {
      return this->equals(other);
    }

    //! returns true if caching entity has same id like mmesh entity
    bool operator==(const MMeshEntity& entity) const
    {
      return this->id_ == entity.impl().grid()->globalIdSet().id( entity );
    }

    //! returns true if id of other is greater
    bool operator<(const MMeshInterfaceCachingEntity& other) const
    {
      return this->id_ < other.id_;
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
      return Geometry( GeometryTypes::simplex(dim), this->vertex_ );
    }

    //! Return the number of subEntities of codimension cc
    unsigned int subEntities (unsigned int cc) const
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

    //! returns true if Entity has no children
    bool isLeaf() const {
      return false;
    }

    //calculates the intersection volume with another entity
    ctype intersectionVolume ( const MMeshEntity& entity ) const
    {
      // assuming simple components (of 2 entities)
      return std::min( geometry().volume(), entity.geometry().volume() );
    }

  }; // end of MMeshInterfaceCachingEntity codim = 0

} // namespace Dune

#endif
