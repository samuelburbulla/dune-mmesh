// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_GRID_COMMON_HH
#define DUNE_MMESH_GRID_COMMON_HH

/** \file
 * \brief Some common helper methods
 */

namespace Dune
{

  struct HashUIntVector {
    std::size_t operator() (const std::vector<unsigned int>& a) const
    {
      std::size_t hash = std::hash<unsigned int>{}(a[0]);
      for ( int i = 1; i < a.size(); ++i )
        hash = hash ^ (std::hash<unsigned int>{}(a[i]) << i);
      return hash;
    }
  };

}

#endif
