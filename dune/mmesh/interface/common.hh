// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_COMMON_HH
#define DUNE_MMESH_INTERFACE_COMMON_HH

/** \file
 * \brief Some common helper methods
 */

namespace Dune
{

  namespace MMeshInterfaceImpl
  {

    //! Return list of indices sorted by id
    template < typename HostEntity, int dim >
    static inline auto computeCGALIndices( const HostEntity& hostEntity )
    {
      std::array< std::pair< std::size_t, std::size_t >, dim+1 > pairs;
      for ( std::size_t i = 0; i < dim+1; ++i )
      {
        auto i0 = (hostEntity.second+i+1)%(dim+2);
        pairs[i].first = i0;
        pairs[i].second = hostEntity.first->vertex( i0 )->info().id;
      }

      std::sort(pairs.begin(), pairs.end(), [](auto& a, auto& b){ return a.second < b.second; });

      std::array< std::size_t, dim+1 > indices;
      for ( std::size_t i = 0; i < dim+1; ++i )
        indices[i] = pairs[i].first;

      assert( hostEntity.first->vertex(indices[0])->info().id < hostEntity.first->vertex(indices[1])->info().id );

      if constexpr (dim == 2)
        assert( hostEntity.first->vertex(indices[1])->info().id < hostEntity.first->vertex(indices[2])->info().id );

      return indices;
    }

  }

}

#endif
