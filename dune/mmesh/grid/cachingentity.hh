// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_CACHINGENTITY_HH
#define DUNE_MMESH_GRID_CACHINGENTITY_HH

/** \file
 * \brief The MMeshCachingEntity class
 */

#include <set>

// Dune includes
#include <dune/grid/common/grid.hh>

// MMesh includes
#include <dune/mmesh/grid/multiid.hh>
#include <dune/mmesh/grid/pointfieldvector.hh>
#include <dune/mmesh/grid/polygoncutting.hh>

namespace Dune
{

  //**********************************************************************
  //
  // --MMeshCachingEntity
  // --Entity
  //
  /** \brief The implementation of caching entities in a MMesh
   *   \ingroup MMesh
   *  The caching entity copys the CGAL face object instead of holding a Face_handle pointer.
   *
   *  A Grid is a container of grid entities. An entity is parametrized by the codimension.
   *  An entity of codimension c in dimension d is a d-c dimensional object.
   *
   */
  template<int dim, class GridImp>
  class MMeshCachingEntity<0,dim,GridImp> :
    public MMeshEntity<0,dim,GridImp>
  {
    template <class GridImp_>
    friend class MMeshLeafIndexSet;

    template <class GridImp_>
    friend class MMeshLocalIdSet;

    template <class GridImp_>
    friend class MMeshGlobalIdSet;

  private:
    // this type
    typedef MMeshCachingEntity<0,dim,GridImp> ThisType;

    // base type
    typedef MMeshEntity<0,dim,GridImp> BaseType;

    // type of scalars
    typedef typename GridImp::ctype ctype;

    // type of the host grid
    typedef typename GridImp::HostGridType HostGrid;

    // equivalent entity in the host grid as pointer
    typedef typename GridImp::template HostGridEntity<0> HostGridEntity;

    // standard MMesh entity implementation
    typedef typename GridImp::template Codim<0>::Entity MMeshEntityType;

    // type of ids
    typedef MMeshImpl::MultiId IdType;

  public:
    // geometry type
    typedef AffineGeometry<ctype, dim, dim> Geometry;

    // local geometry type
    typedef typename GridImp::template Codim<0>::LocalGeometry LocalGeometry;

    // type of global coordinate
    typedef typename Geometry::GlobalCoordinate GlobalCoordinate;

    MMeshCachingEntity() = delete;

    explicit MMeshCachingEntity(const GridImp* mMesh, const HostGridEntity& hostEntity)
      : BaseType(mMesh, hostEntity, mMesh->globalIdSet().id( mMesh->entity( hostEntity ) ))
    {
      const auto geo = mMesh->entity( hostEntity ).geometry();
      for( int i = 0; i < dim+1; ++i )
        this->vertex_[i] = geo.corner(i);
    }

    //! returns true if host entities are equal
    bool equals(const MMeshCachingEntity& other) const
    {
      return this->id_ == other.id_;
    }

    //! returns true if host entities are equal
    bool operator==(const MMeshCachingEntity& other) const
    {
      return this->equals(other);
    }

    //! returns true if caching entity has same id like mmesh entity
    bool operator==(const MMeshEntityType& entity) const
    {
      return this->id_ == this->mMesh_->globalIdSet().id( entity );
    }

    //! returns true if id of other is greater
    bool operator<(const MMeshCachingEntity& other) const
    {
      return this->id_ < other.id_;
    }

    //! returns true if father entity exists
    bool hasFather () const
    {
      return true;
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

    //! calculates the intersection volume with another MMesh entity
    template<int d = dim>
    std::enable_if_t<d == 2, ctype>
    intersectionVolume ( const MMeshEntityType& entity ) const
    {
      std::array<GlobalCoordinate, 3> entityPoints;

      for ( int i = 0; i < 3; ++i )
      {
        entityPoints[i] =
          makeFieldVector( entity.impl().hostEntity()->vertex(i)->point() );
      }

      using PC = Dune::PolygonCutting<ctype, GlobalCoordinate>;
      return std::abs( PC::polygonArea( PC::sutherlandHodgman(this->vertex_, entityPoints) ) );
    }

    template<int d = dim>
    std::enable_if_t<d == 3, ctype>
    intersectionVolume( const MMeshEntityType& entity ) const
    {
      DUNE_THROW( NotImplemented, "intersectionVolume in 3d" );
    }

  }; // end of MMeshCachingEntity codim = 0

} // namespace Dune

#endif
