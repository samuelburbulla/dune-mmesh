// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_MMESH_CC
#define DUNE_MMESH_MMESH_CC

/** \file
 * \brief The MMesh class
 */

#include "mmesh.hh"

namespace Dune
{
  // MMesh (2D)
  template <class HostGrid>
  class MMesh<HostGrid, 2>;

  // MMesh (3D)
  template <class HostGrid>
  class MMesh<HostGrid, 3>;

} // namespace Dune

#endif // DUNE_MMESH_MMESH_CC
