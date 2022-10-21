#pragma GCC diagnostic ignored "-Wattributes"
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

  template< class Grid, class MMeshType >
  class MMeshCommunication
  {
    typedef MMeshCommunication< Grid, MMeshType > This;
    typedef PartitionHelper< MMeshType > PartitionHelperType;
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

      // Send to all links
      const auto& comm = partitionHelper_.comm();
      MPI_Request sendRequests[links.size()];
      for (int link = 0; link < links.size(); ++link)
      {
        BufferType& buf = sendBuffers[ link ];
        int dest = links[ link ];
        MPI_Request& request = sendRequests[ link ];
        MPI_Isend( buf._buf, buf._wb, MPI_BYTE, dest, tag_, comm, &request );
      }

      // Receive data
      int count = 0;
      std::vector<bool> received (links.size(), false);
      MPI_Request recvRequests[links.size()];
      while( count < links.size() )
      {
        for (int link = 0; link < links.size(); ++link)
        {
          if (received[ link ])
            continue;

          int source = links[ link ];

          int available = 0;
          MPI_Status status;
          MPI_Iprobe( source, tag_, comm, &available, &status );

          if (available)
          {
            int bufferSize;
            MPI_Get_count( &status, MPI_BYTE, &bufferSize );

            BufferType& buf = recvBuffers[ link ];
            buf.reserve( bufferSize );
            buf.clear();

            MPI_Request& request = recvRequests[ link ];
            MPI_Irecv( buf._buf, bufferSize, MPI_BYTE, source, tag_, comm, &request );
            buf.seekp( bufferSize );

            count++;
            received[ link ] = true;
          }
        }
      }

      MPI_Waitall( links.size(), sendRequests, MPI_STATUSES_IGNORE );
      MPI_Waitall( links.size(), recvRequests, MPI_STATUSES_IGNORE );

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

  template< class Grid, class MMeshType >
  template< int codim >
  struct MMeshCommunication< Grid, MMeshType >::PackData
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
          if (connectivity.count(links[link]) > 0 || links[link] == partitionHelper.rank(element))
          {
            std::size_t size = dataHandle.size( entity );

            // write size into stream
            buffer[ link ].write( size );

            // write data to message buffer using data handle
            dataHandle.gather( buffer[ link ], entity );
          }
        }
      }
    }
  };



  // MMeshCommunication::UnpackData
  // -------------------------------------

  template< class Grid, class MMeshType >
  template< int codim >
  struct MMeshCommunication< Grid, MMeshType >::UnpackData
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

      const auto& connectivity = partitionHelper.connectivity(element);

      // get number of sub entities
      const int numSubEntities = element.subEntities(codim);
      for( int subEntity = 0; subEntity < numSubEntities; ++subEntity )
      {
        // get subentity
        const auto& entity = element.template subEntity< codim >( subEntity );

        for (int link = 0; link < links.size(); ++link)
        {
          // make sure entity belongs to the rank of the link
          if (links[link] == partitionHelper.rank(element) || connectivity.count(links[link]) > 0)
          {
            // read size from stream
            std::size_t size( 0 );
            buffer[ link ].read( size );

            // read data from message buffer using data handle
            dataHandle.scatter( buffer[ link ], entity, size );
          }
        }
      }
    }
  };

} // namespace Dune

#endif
