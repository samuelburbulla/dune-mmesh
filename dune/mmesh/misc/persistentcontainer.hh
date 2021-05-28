#ifndef DUNE_MMESH_MISC_PERSISTENTCONTAINER_HH
#define DUNE_MMESH_MISC_PERSISTENTCONTAINER_HH

#include <unordered_map>
#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/persistentcontainermap.hh>

namespace Dune
{

  template< class HostGrid, int dim, class T >
  class PersistentContainer< MMesh< HostGrid, dim >, T >
    : public Dune::PersistentContainerMap<
        MMesh< HostGrid, dim >,
        typename MMesh< HostGrid, dim >::LocalIdSet,
        std::unordered_map< typename MMesh< HostGrid, dim >::LocalIdSet::IdType, T > >
  {
    typedef MMesh< HostGrid, dim > G;
    typedef PersistentContainerMap< G, typename G::LocalIdSet, std::unordered_map< typename G::LocalIdSet::IdType, T > > Base;
    typedef typename MMesh< HostGrid, dim >::LocalIdSet IdSet;
    typedef std::unordered_map< typename MMesh< HostGrid, dim >::LocalIdSet::IdType, T > Map;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, grid.localIdSet(), value )
    {}

    using Base::codimension;
    using Base::migrateLevel;

    void resize ( const Value &value = Value() )
    {
      if( sequence_ != this->grid_->sequence() )
      {
        Hybrid::forEach( std::make_index_sequence< Grid::dimension+1 >{},
          [ & ]( auto i ){ if( i == this->codimension() ) this->template resize< i >( value ); } );
      }
    }

    template< int codim >
    void resize ( const Value &value )
    {
      std::integral_constant< bool, Capabilities::hasEntity< Grid, codim >::v > hasEntity;
      assert( codim == codimension() );

      // add one id for caching entity during adaptation
      MMeshImpl::MultiId id ( {42, 42, 42} );
      this->data_.insert( std::make_pair( id, value ) );

      // add one id for a codim 1 entity during adaptation
      MMeshImpl::MultiId eid ( {42, 42} );
      this->data_.insert( std::make_pair( eid, value ) );

      // add one id for caching vertices during adaptation
      auto vh = this->grid_->getHostGrid().infinite_vertex();
      MMeshImpl::MultiId vid ( vh->info().id );
      this->data_.insert( std::make_pair( vid, value ) );

      // create empty map, but keep old data
      Map data;

      // add new entries
      this->template migrateLevel< codim >( 0, value, data, hasEntity );
    }

  private:
    int sequence_ = -1;

  };

  template< class MMesh, class T >
  class PersistentContainer< MMeshInterfaceGrid< MMesh >, T >
    : public Dune::PersistentContainerMap<
        MMeshInterfaceGrid< MMesh >,
        typename MMeshInterfaceGrid< MMesh >::LocalIdSet,
        std::unordered_map< typename MMeshInterfaceGrid< MMesh >::LocalIdSet::IdType, T > >
  {
    typedef MMeshInterfaceGrid< MMesh > G;
    typedef PersistentContainerMap< G, typename G::LocalIdSet, std::unordered_map< typename G::LocalIdSet::IdType, T > > Base;
    typedef typename MMeshInterfaceGrid< MMesh >::LocalIdSet IdSet;
    typedef std::unordered_map< typename MMeshInterfaceGrid< MMesh >::LocalIdSet::IdType, T > Map;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
      : Base( grid, codim, grid.localIdSet(), value )
    {}

    using Base::codimension;
    using Base::migrateLevel;

    void resize ( const Value &value = Value() )
    {
      if( sequence_ != this->grid_->getMMesh().sequence() )
      {
        Hybrid::forEach( std::make_index_sequence< Grid::dimension+1 >{},
          [ & ]( auto i ){ if( i == this->codimension() ) this->template resize< i >( value ); } );
        sequence_ = this->grid_->getMMesh().sequence();
      }
    }

    template< int codim >
    void resize ( const Value &value )
    {
      std::integral_constant< bool, Capabilities::hasEntity< Grid, codim >::v > hasEntity;
      assert( codim == codimension() );

      // add one id for caching entity during adaptation
      MMeshImpl::MultiId id ( {42, 42} );
      this->data_.insert( std::make_pair( id, value ) );

      // add one id for caching vertices during adaptation
      auto vh = this->grid_->getHostGrid().infinite_vertex();
      MMeshImpl::MultiId vid ( vh->info().id );
      this->data_.insert( std::make_pair( vid, value ) );

      // create empty map, but keep old data
      Map data;

      // add new entries
      this->template migrateLevel< codim >( 0, value, data, hasEntity );
    }

  private:
    int sequence_ = -1;

  };

} // namespace Dune

#endif // #ifndef DUNE_MMESH_MISC_PERSISTENTCONTAINER_HH
