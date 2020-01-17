// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_COMMON_HH
#define DUNE_MMESH_GRID_COMMON_HH

/** \file
 * \brief Some common helper methods
 */

namespace Dune
{

  //! Hash a UInt vector
  struct HashUIntVector {
    std::size_t operator() (const std::vector<std::size_t>& a) const
    {
      std::size_t hash = std::hash<std::size_t>{}(a[0]);
      for ( std::size_t i = 1; i < a.size(); ++i )
        hash = hash ^ (std::hash<std::size_t>{}(a[i]) << i);
      return hash;
    }
  };

  //! Hash a UInt array
  struct HashUIntArray {
    template< std::size_t dim >
    std::size_t operator() (const std::array<std::size_t, dim>& a) const
    {
      std::size_t hash = std::hash<std::size_t>{}(a[0]);
      for ( std::size_t i = 1; i < a.size(); ++i )
        hash = hash ^ (std::hash<std::size_t>{}(a[i]) << i);
      return hash;
    }
  };
}

#endif
