// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_ENTITYSEED_HH
#define DUNE_MMESH_GRID_ENTITYSEED_HH

/**
 * \file
 * \brief The MMeshEntitySeed class
 */

namespace Dune
{

  /**
   * \brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
   * \ingroup MMesh
   *  We use the host grid entity as seed.
   */
  template<int codim, class GridImp>
  class MMeshEntitySeed
  {
  protected:
    // Entity type of the hostgrid
    typedef typename GridImp::template HostGridEntity<codim> HostGridEntity;

  public:
    enum {codimension = codim};

    /**
     * \brief Construct empty seed
     */
    MMeshEntitySeed() {}

    /**
     * \brief Construct a seed from a host entity.
     */
    MMeshEntitySeed(const HostGridEntity& hostEntity)
     : hostEntity_( hostEntity )
    {}

    /**
     * \brief Construct a seed from another seed.
     */
    MMeshEntitySeed& operator= (const MMeshEntitySeed<codim, GridImp>& seed)
    {
      if( this != &seed )
        hostEntity_ = seed.hostEntity_;
      return *this;
    }

    /**
     * \brief Check whether it is safe to create an Entity from this Seed
     */
    template< int cd = codim >
    std::enable_if_t< cd != 2 || GridImp::dimension != 3, bool >
    isValid() const
    {
      return hostEntity_ != HostGridEntity();
    }

    //! Special handling for codim 2 entities in 3d
    template< int cd = codim >
    std::enable_if_t< cd == 2 && GridImp::dimension == 3, bool >
    isValid() const
    {
      return hostEntity_.first != HostGridEntity().first;
    }

    /**
     * \brief Return host entity
     */
    const HostGridEntity& hostEntity() const
    {
      return hostEntity_;
    }

  private:
    //! The host entity
    HostGridEntity hostEntity_;
  };

} // namespace Dune

#endif
