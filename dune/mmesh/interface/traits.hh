// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_INTERFACE_TRAITS_HH
#define DUNE_MMESH_INTERFACE_TRAITS_HH

/** \file
 * \brief The MMeshInterfaceGrid traits class
 */

#include <config.h>

namespace Dune
{
  // Forward declarations
  template<int dim, class MMesh>
  struct MMeshInterfaceGridFamily;

  template<class MMesh, int dim, class GridFamily = MMeshInterfaceGridFamily<dim, MMesh>>
  class MMeshInterfaceGrid;

  template<int dim, class MMesh>
  struct MMeshInterfaceGridFamily;

  template<int codim, int dim, class GridImp>
  class MMeshInterfaceGridEntity;

  template<class GridImp>
  class MMeshInterfaceGridLeafIntersectionIterator;

  template<class GridImp>
  class MMeshInterfaceGridHierarchicIterator;

  template<int cd, int dim, class GridImp>
  class MMeshInterfaceConnectedComponent;

} // namespace Dune

#endif // DUNE_MMESH_INTERFACE_TRAITS_HH
