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

  const int rank = comm.rank();
  const int size = comm.size();
  std::vector<int> links;
  for (int i = 0; i < size; i++)
    if (i != rank)
      links.push_back( i );

  // Create data
  std::vector<int> data (comm.size());
  for (int i = 0; i < data.size(); ++i)
    data[ i ] = rank;

  // Write buffers
  using BufferType = Dune::MMeshImpl::ObjectStream;
  std::vector<BufferType> sendBuffers (links.size());
  std::vector<BufferType> recvBuffers (links.size());

  for (int link = 0; link < links.size(); ++link)
  {
    auto r = links[ link ];
    int d = data[ r ];
    sendBuffers[ link ].write( d );
  }

  std::vector<int> invertedLinks (size);
  for (int l = 0; l < links.size(); ++l)
    invertedLinks[ links[ l ] ] = l;

  int tag = 0;

  // Send
  for (int r = 0; r < comm.size(); ++r)
  {
    if (r == rank)
      for (int link = 0; link < links.size(); ++link)
      {
        BufferType& buf = sendBuffers[ link ];
        int dest = links[ link ];
        MPI_Ssend(buf._buf, buf._wb, MPI_BYTE, dest, tag, comm);
      }
    else
    {
      int link = invertedLinks[ r ];
      BufferType& buf = recvBuffers[ link ];

      MPI_Status status;
      MPI_Probe(r, tag, MPI_COMM_WORLD, &status);

      int bufferSize;
      MPI_Get_count(&status, MPI_BYTE, &bufferSize);
      buf.reserve(bufferSize);

      MPI_Recv(buf._buf, bufferSize, MPI_BYTE, r, tag, comm, &status);
      buf.seekp(bufferSize);
    }
  }

  for (int link = 0; link < links.size(); ++link)
  {
    int e;
    auto r = links[ link ];
    recvBuffers[ link ].read( e );
    std::cout << rank << ": " << e << std::endl;
    assert( e == r );
  }

  return EXIT_SUCCESS;
}
