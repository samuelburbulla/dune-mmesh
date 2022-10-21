// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <mpi.h>
#include <iostream>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/timer.hh>
#include <dune/common/parallel/mpicommunication.hh>
#include <dune/mmesh/mmesh.hh>

int main(int argc, char *argv[])
{
  Dune::MPIHelper::instance(argc, argv);

  using Comm = Dune::Communication<MPI_Comm>;
  auto comm = Comm( Dune::MPIHelper::getCommunicator() );

  // Create links
  const int rank = comm.rank();
  const int size = comm.size();
  std::vector<int> links;
  if (rank > 0)
    links.push_back( rank-1 );
  if (rank < size-1)
    links.push_back( rank+1 );

  // Create data (my rank for all links)
  using DataType = std::array<int, 1500>;
  std::vector<DataType> data (links.size());
  for (int link = 0; link < links.size(); ++link)
  {
    DataType d;
    std::fill(d.begin(), d.end(), rank);
    data[ link ] = d;
  }

  // Write buffers
  using BufferType = Dune::MMeshImpl::ObjectStream;
  std::vector<BufferType> sendBuffers (links.size());
  std::vector<BufferType> recvBuffers (links.size());

  for (int link = 0; link < links.size(); ++link)
    sendBuffers[ link ].write( data[ link ] );

  int tag = 0;
  Dune::Timer timer;
  timer.start();

  // Send
  MPI_Request sendRequests[links.size()];
  for (int link = 0; link < links.size(); ++link)
  {
    BufferType& buf = sendBuffers[ link ];
    int dest = links[ link ];
    MPI_Request& request = sendRequests[ link ];
    MPI_Isend(buf._buf, buf._wb, MPI_BYTE, dest, tag, comm, &request);
  }

  // Receive
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
      MPI_Iprobe( source, tag, comm, &available, &status );

      if (available)
      {
        int bufferSize;
        MPI_Get_count(&status, MPI_BYTE, &bufferSize);

        BufferType& buf = recvBuffers[ link ];
        buf.reserve(bufferSize);
        buf.clear();

        MPI_Request& request = recvRequests[ link ];
        MPI_Irecv(buf._buf, bufferSize, MPI_BYTE, source, tag, comm, &request);
        buf.seekp(bufferSize);

        count++;
        received[ link ] = true;
      }
    }
  }

  MPI_Waitall( links.size(), sendRequests, MPI_STATUSES_IGNORE );
  MPI_Waitall( links.size(), recvRequests, MPI_STATUSES_IGNORE );

  auto maxT = comm.max( timer.elapsed() );
  if (comm.rank() == 0)
    std::cout << "Comm took " << maxT << std::endl;

  for (int link = 0; link < links.size(); ++link)
  {
    DataType e;
    recvBuffers[ link ].read( e );
    data[ link ] = e;
  }

  // Check data
  for (int link = 0; link < links.size(); ++link)
    assert( data[ link ][0] == links[link]  );

  return EXIT_SUCCESS;
}
