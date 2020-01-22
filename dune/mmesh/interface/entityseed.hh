// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_ENTITYSEED_HH
#define DUNE_MMESH_INTERFACE_ENTITYSEED_HH

/**
 * \file
 * \brief The MMeshInterfaceGridEntitySeed class
 */

namespace Dune
{

  /**
   * \brief The EntitySeed class provides the minimal information needed to restore an Entity using the grid.
   * \ingroup MMesh
   *  We use the host grid entity as seed.
   */
  template<int codim, class GridImp>
  class MMeshInterfaceGridEntitySeed
  {
  protected:
    // Entity type of the hostgrid
    typedef typename GridImp::template MMeshInterfaceEntity<codim> HostGridEntity;

  public:
    enum {codimension = codim};

    /**
     * \brief Construct empty seed
     */
    MMeshInterfaceGridEntitySeed() : valid_( false ) {}

    /**
     * \brief Construct a seed from a host entity.
     */
    MMeshInterfaceGridEntitySeed(const HostGridEntity& hostEntity)
     : hostEntity_( hostEntity ), valid_( true )
    {}

    /**
     * \brief Construct a seed from another seed.
     */
    MMeshInterfaceGridEntitySeed operator= (const MMeshInterfaceGridEntitySeed<codim, GridImp>& seed)
    {
      hostEntity_ = seed.hostEntity_;
      valid_ = seed.valid_;
      return seed;
    }

    /**
     * \brief Check whether it is safe to create an Entity from this Seed
     */
    bool isValid() const
    {
      return valid_;
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
    bool valid_;
  };

} // namespace Dune

#endif
