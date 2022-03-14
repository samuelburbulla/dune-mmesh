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
        sequence_ = this->grid_->sequence();
      }
    }

    template< int codim >
    void resize ( const Value &value )
    {
      std::integral_constant< bool, Capabilities::hasEntity< Grid, codim >::v > hasEntity;
      assert( codim == codimension() );

      // add one id for caching entity during adaptation
      if constexpr (codim == 0)
      {
        MMeshImpl::MultiId id ( { std::size_t(-4), std::size_t(-3), std::size_t(-2) } );
        this->data_.insert( std::make_pair( id, value ) );
      }

      // add one id for a every caching codim 1 entity during adaptation
      if constexpr (codim == 1)
      {
        MMeshImpl::MultiId eid0 ( { std::size_t(-4), std::size_t(-3) } );
        this->data_.insert( std::make_pair( eid0, value ) );
        MMeshImpl::MultiId eid1 ( { std::size_t(-4), std::size_t(-2) } );
        this->data_.insert( std::make_pair( eid1, value ) );
        MMeshImpl::MultiId eid2 ( { std::size_t(-3), std::size_t(-2) } );
        this->data_.insert( std::make_pair( eid2, value ) );
      }

      // add one id for every caching vertex during adaptation
      if constexpr (codim == dim)
      {
        MMeshImpl::MultiId vid0 ( std::size_t(-4) );
        this->data_.insert( std::make_pair( vid0, value ) );
        MMeshImpl::MultiId vid1 ( std::size_t(-3) );
        this->data_.insert( std::make_pair( vid1, value ) );
        MMeshImpl::MultiId vid2 ( std::size_t(-2) );
        this->data_.insert( std::make_pair( vid2, value ) );
      }

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
    static constexpr int dim = MMesh::dimensionworld-1;

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
        sequence_ = this->grid_->sequence();
      }
    }

    template< int codim >
    void resize ( const Value &value )
    {
      std::integral_constant< bool, Capabilities::hasEntity< Grid, codim >::v > hasEntity;
      assert( codim == codimension() );

      // add one id for caching entity during adaptation
      if constexpr (codim == 0)
      {
        MMeshImpl::MultiId id ( { std::size_t(-3), std::size_t(-2)} );
        this->data_.insert( std::make_pair( id, value ) );
      }

      // add one id for every caching vertex during adaptation
      if constexpr (codim == dim)
      {
        MMeshImpl::MultiId vid0 ( std::size_t(-3) );
        this->data_.insert( std::make_pair( vid0, value ) );
        MMeshImpl::MultiId vid1 ( std::size_t(-2) );
        this->data_.insert( std::make_pair( vid1, value ) );
      }

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
