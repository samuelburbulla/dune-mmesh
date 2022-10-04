#ifndef DUNE_MMESH_MISC_COMMUNICATION_HH
#define DUNE_MMESH_MISC_COMMUNICATION_HH

#include <vector>
#include <mpi.h>
#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/common/parallel/variablesizecommunicator.hh>

#include "objectstream.hh"

namespace Dune
{

  template< class Grid >
  class MMeshCommunication
  {
    typedef MMeshCommunication< Grid > This;
    typedef PartitionHelper< Grid > PartitionHelperType;
    typedef typename PartitionHelperType::LinksType Links;

    // prohibit copying and assignment
    MMeshCommunication( const This & );
    const This &operator= ( const This & );

    template< int codim >
    struct PackData;

    template< int codim >
    struct UnpackData;

  public:
    static const int dimension = Grid::dimension;

    MMeshCommunication ( const PartitionHelperType &partitionHelper )
    : partitionHelper_( partitionHelper ), tag_( 0 )
    {}

    template< class PackIterator, class UnpackIterator, class DataHandleImp >
    void operator() ( const PackIterator packBegin, const PackIterator packEnd,
                      const UnpackIterator unpackBegin, const UnpackIterator unpackEnd,
                      DataHandleImp &dataHandle,
                      const PartitionType sendType, const PartitionType recvType,
                      const bool packAll ) const
    {
      typedef MMeshImpl::ObjectStream BufferType;
      const Links& links = partitionHelper_.links();

      // vector of message buffers
      std::vector< BufferType > sendBuffers( links.size() ), recvBuffers( links.size() );
      for( int link = 0; link < links.size(); ++link )
      {
        sendBuffers[ link ].clear();
        recvBuffers[ link ].clear();
      }

      // pack data on send entities
      for( PackIterator it = packBegin; it != packEnd; ++it )
      {
        const typename PackIterator::Entity &entity = *it;
        if (entity.partitionType() == sendType && partitionHelper_.connectivity(entity).size() > 0)
        {
          Hybrid::forEach(std::make_index_sequence<dimension+1>{}, [&](auto codim){
            PackData<codim>::apply(links, partitionHelper_, dataHandle, sendBuffers, entity);
          });
        }
      }

      if( packAll )
      {
        // pack data on receive entities
        for( UnpackIterator it = unpackBegin; it != unpackEnd; ++it )
        {
          const typename UnpackIterator::Entity &entity = *it;
          if (entity.partitionType() == recvType && partitionHelper_.connectivity(entity).size() > 0)
          {
            Hybrid::forEach(std::make_index_sequence<dimension+1>{}, [&](auto codim){
              PackData<codim>::apply(links, partitionHelper_, dataHandle, sendBuffers, entity);
            });
          }
        }
      }

      // Invert links
      const auto& comm = partitionHelper_.comm();
      std::vector<int> invertedLinks (comm.size());
      for (int l = 0; l < links.size(); ++l)
        invertedLinks[ links[ l ] ] = l;

      // Exchange data
      for (int r = 0; r < comm.size(); ++r)
      {
        if (r == comm.rank())
        {
          // Send to all ranks
          for (int link = 0; link < links.size(); ++link)
          {
            BufferType& buf = sendBuffers[ link ];
            int dest = links[ link ];
            MPI_Ssend( buf._buf, buf._wb, MPI_BYTE, dest, tag_, comm );
          }
        }
        else
        {
          // Receive
          int link = invertedLinks[ r ];
          BufferType& buf = recvBuffers[ link ];

          MPI_Status status;
          MPI_Probe( r, tag_, comm, &status );

          int bufferSize;
          MPI_Get_count( &status, MPI_BYTE, &bufferSize );
          buf.reserve( bufferSize );
          buf.clear();

          MPI_Recv( buf._buf, bufferSize, MPI_BYTE, r, tag_, comm, &status );
          buf.seekp( bufferSize );

          buf.seekp( bufferSize );
        }
      }

      // unpack data on receive entities
      for( UnpackIterator it = unpackBegin; it != unpackEnd; ++it )
      {
        const typename UnpackIterator::Entity &entity = *it;
        if (entity.partitionType() == recvType && partitionHelper_.connectivity(entity).size() > 0)
        {
          Hybrid::forEach(std::make_index_sequence<dimension+1>{}, [&](auto codim){
              UnpackData<codim>::apply(links, partitionHelper_, dataHandle, recvBuffers, entity);
          });
        }
      }

      if( packAll )
      {
        // unpack data on send entities
        for( PackIterator it = packBegin; it != packEnd; ++it )
        {
          const typename PackIterator::Entity &entity = *it;
          if (entity.partitionType() == sendType && partitionHelper_.connectivity(entity).size() > 0)
          {
            Hybrid::forEach(std::make_index_sequence<dimension+1>{}, [&](auto codim){
                UnpackData<codim>::apply(links, partitionHelper_, dataHandle, recvBuffers, entity);
            });
          }
        }
      }

      tag_++;
      if (tag_ < 0) tag_ = 0;
    }

  private:
    const PartitionHelperType &partitionHelper_;
    mutable int tag_;
  };

  // MMeshCommunication::PackData
  // -----------------------------------

  template< class Grid >
  template< int codim >
  struct MMeshCommunication< Grid >::PackData
  {
    typedef typename Grid::template Codim< 0 >::Entity Element;

    typedef typename Grid::template Codim< codim >::Entity Entity;

    template< class DataHandleIF, class BufferType >
    static void apply ( const Links& links,
                        const PartitionHelperType &partitionHelper,
                        DataHandleIF &dataHandle,
                        std::vector< BufferType > &buffer,
                        const Element &element )
    {
      // if codim is not contained just go on
      if( !dataHandle.contains( dimension, codim ) )
        return;

      const auto& connectivity = partitionHelper.connectivity(element);

      const int numSubEntities = element.subEntities(codim);
      for( int subEntity = 0; subEntity < numSubEntities; ++subEntity )
      {
        // get subentity
        const Entity &entity = element.template subEntity< codim >( subEntity );

        for (int link = 0; link < links.size(); ++link)
        {
          // make sure entity belongs to the link
          if (connectivity.count(links[link]) == 0)
            continue;

          std::size_t size = dataHandle.size( entity );

          // write size into stream
          buffer[ link ].write( size );

          // write data to message buffer using data handle
          dataHandle.gather( buffer[ link ], entity );
        }
      }
    }
  };



  // MMeshCommunication::UnpackData
  // -------------------------------------

  template< class Grid >
  template< int codim >
  struct MMeshCommunication< Grid >::UnpackData
  {
    using Element = typename Grid::template Codim< 0 >::Entity;

    template< class DataHandleIF, class BufferType >
    static void apply ( const Links& links,
                        const PartitionHelperType &partitionHelper,
                        DataHandleIF &dataHandle,
                        std::vector< BufferType > &buffer,
                        const Element &element )
    {
      // if codim is not contained just go on
      if( !dataHandle.contains( dimension, codim ) )
        return;

      // get number of sub entities
      const int numSubEntities = element.subEntities(codim);
      for( int subEntity = 0; subEntity < numSubEntities; ++subEntity )
      {
        // get subentity
        const auto& entity = element.template subEntity< codim >( subEntity );

        for (int link = 0; link < links.size(); ++link)
        {
          // make sure entity belongs to the rank of the link
          if (links[link] != element.impl().hostEntity()->info().rank)
            continue;

          // read size from stream
          std::size_t size( 0 );
          buffer[ link ].read( size );

          // read data from message buffer using data handle
          dataHandle.scatter( buffer[ link ], entity, size );
        }
      }
    }
  };

} // namespace Dune

#endif
