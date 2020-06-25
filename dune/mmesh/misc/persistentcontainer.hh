#ifndef DUNE_MMESH_MISC_PERSISTENTCONTAINER_HH
#define DUNE_MMESH_MISC_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/persistentcontainermap.hh>

namespace Dune
{

  template< class HostGrid, int dim, class T >
  class PersistentContainer< MMesh< HostGrid, dim >, T >
    : public Dune::PersistentContainerMap<
        MMesh< HostGrid, dim >,
        typename MMesh< HostGrid, dim >::LocalIdSet,
        std::map< typename MMesh< HostGrid, dim >::LocalIdSet::IdType, T > >
  {
    typedef MMesh< HostGrid, dim > G;
    typedef PersistentContainerMap< G, typename G::LocalIdSet, std::map< typename G::LocalIdSet::IdType, T > > Base;
    typedef typename MMesh< HostGrid, dim >::LocalIdSet IdSet;
    typedef std::map< typename MMesh< HostGrid, dim >::LocalIdSet::IdType, T > Map;

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
      Hybrid::forEach( std::make_index_sequence< Grid::dimension+1 >{},
        [ & ]( auto i ){ if( i == codimension() ) resize< i >( value ); } );
    }

    template< int codim >
    inline void resize ( const Value &value )
    {
      std::integral_constant< bool, Capabilities::hasEntity< Grid, codim >::v > hasEntity;
      assert( codim == codimension() );

      // create empty map, but keep old data
      Map data;

      // add new entries
      this->template migrateLevel< codim >( 0, value, data, hasEntity );
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_MMESH_MISC_PERSISTENTCONTAINER_HH