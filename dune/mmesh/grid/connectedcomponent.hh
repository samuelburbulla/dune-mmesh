// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_CONNECTEDCOMPONENT_HH
#define DUNE_MMESH_GRID_CONNECTEDCOMPONENT_HH

/** \file
 * \brief The MMeshConnectedComponent class
 */

#include <memory>
#include <unordered_map>

// MMesh includes
#include <dune/mmesh/grid/multiid.hh>

namespace Dune
{
  // Forward declarations
  template<int codim, int dim, class GridImp>
  class MMeshCachingEntity;

  //**********************************************************************
  //
  // --MMeshConnectedComponent
  // --Entity
  //
  /** \brief The implementation of a connected component of entities in MMesh
   *   \ingroup MMesh
   *  The connected component stores a list of connected entities providing geometrical information for the remeshing step.
   *
   */
  template<class GridImp>
  class MMeshConnectedComponent
  {
    // dimension
    static constexpr int dim = GridImp::dimension;

    // type of scalars
    using ctype = typename GridImp::ctype;

    // type of the entity
    using Entity = typename GridImp::template Codim<0>::Entity;

    // type of the caching entity
    using CachingEntity = MMeshCachingEntity< 0, dim, const GridImp >;

    //type of ids
    using IdType = MMeshImpl::MultiId;

  public:
    MMeshConnectedComponent() : componentNumber_(0) {}

    MMeshConnectedComponent( const GridImp* mMesh, const Entity& entity )
     : mMesh_( mMesh ),
       componentNumber_( entity.impl().hostEntity()->info().componentNumber )
    {
      entities_.emplace_back( mMesh_, entity.impl().hostEntity() );
      CachingEntity& cachingEntity = entities_.back();

      const IdType id = mMesh_->globalIdSet().id( entity );
      entityIdToCachingPtr_.insert( std::make_pair( id, &cachingEntity ) );

      insertNeighbors_( entity, cachingEntity );
    }

    MMeshConnectedComponent& operator=(const MMeshConnectedComponent& other)
    {
      // check for self-assignment
      if(&other == this)
          return *this;

      entities_ = other.entities_;
      entityIdToCachingPtr_ = other.entityIdToCachingPtr_;
      mMesh_ = other.mMesh_;
      componentNumber_ = other.componentNumber_;

      return *this;
    }

    void update( const Entity& entity )
    {
      entities_.emplace_back( mMesh_, entity.impl().hostEntity() );
      CachingEntity& cachingEntity = entities_.back();

      const IdType id = mMesh_->globalIdSet().id( entity );
      entityIdToCachingPtr_.insert( std::make_pair( id, &cachingEntity ) );

      insertNeighbors_( entity, cachingEntity );
    }

    const std::list< CachingEntity >& entities() const
    {
      return entities_;
    }

    //! Return list of caching entities in this component
    const std::list< CachingEntity >& children() const
    {
      return entities();
    }

    //! Return number of caching entities in this component
    const std::size_t size() const
    {
      return entities_.size();
    }

    bool hasEntity( const Entity& entity ) const
    {
      const IdType& id = mMesh_->globalIdSet().id( entity );
      return ( entityIdToCachingPtr_.find( id ) != entityIdToCachingPtr_.end() );
    }

    std::size_t componentNumber() const
    {
      return componentNumber_;
    }

  private:
    /** \brief Insert all might vanishing neighbors of entity into entities_ recursively. */
    void insertNeighbors_(const Entity& entity, CachingEntity& cachingEntity)
    {
      for ( const auto& intersection : intersections( mMesh_->leafGridView(), entity ) )
        if ( intersection.neighbor() )
        {
          const Entity& neighbor = intersection.outside();
          if ( neighbor.impl().hostEntity()->info().componentNumber == componentNumber_ && !neighbor.isNew() )
          {
            const IdType id = mMesh_->globalIdSet().id( neighbor );

            const auto& it = entityIdToCachingPtr_.find( id );

            // if not found, add neighbor and start recursion
            if ( it == entityIdToCachingPtr_.end() )
            {
              entities_.emplace_back( mMesh_, neighbor.impl().hostEntity() );
              CachingEntity& cachingNeighbor = entities_.back();

              entityIdToCachingPtr_.insert( std::make_pair( id, &cachingNeighbor ) );

              insertNeighbors_( neighbor, cachingNeighbor );
            }
          }
        }
    }

    //! The entities of the connected component
    //! Note: we use a list here as we will share pointers to these entities
    std::list< CachingEntity > entities_;
    std::unordered_map< IdType, CachingEntity* > entityIdToCachingPtr_;

    //! pointer to the grid implementation
    const GridImp* mMesh_;
    std::size_t componentNumber_;

  }; // end of MMeshConnectedComponent

} // namespace Dune

#endif
