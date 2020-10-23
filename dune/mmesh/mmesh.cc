// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MMESH_MMESH_CC
#define DUNE_MMESH_MMESH_CC

/** \file
 * \brief The explict MMesh template parameter instantiation
 */

#include "grid/declaration.hh"
#include "cgal/defaults.hh"
#include "cgal/triangulationwrapper.hh"
#include "grid/mmesh.hh"

namespace Dune
{

  // Instantiate explicit template parameters
  template class MMesh<TriangulationWrapper<2>, 2>;
  template class MMeshExplicitGridFactory<MMesh<TriangulationWrapper<2>, 2>>;
  template class MMeshLeafIndexSet<const MMesh<TriangulationWrapper<2>, 2>>;
  template class MMeshGeometry<2, 2, const MMesh<TriangulationWrapper<2>, 2>>;
  template class MMeshGeometry<1, 2, const MMesh<TriangulationWrapper<2>, 2>>;
  template class MMeshGeometry<0, 2, const MMesh<TriangulationWrapper<2>, 2>>;
  template class MMeshEntity<0, 2, const MMesh<TriangulationWrapper<2>, 2>>;
  template class MMeshEntity<1, 2, const MMesh<TriangulationWrapper<2>, 2>>;
  template class MMeshEntity<2, 2, const MMesh<TriangulationWrapper<2>, 2>>;

  template class MMeshInterfaceGrid<MMesh<TriangulationWrapper<2>, 2>>;
  template class MMeshInterfaceGridLeafIndexSet<const MMeshInterfaceGrid<MMesh<TriangulationWrapper<2>, 2>>>;
  template class MMeshInterfaceGridGeometry<1, 2, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<2>, 2>>>;
  template class MMeshInterfaceGridGeometry<0, 2, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<2>, 2>>>;
  template class MMeshInterfaceGridEntity<0, 1, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<2>, 2>>>;
  template class MMeshInterfaceGridEntity<1, 1, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<2>, 2>>>;

  template class MMesh<TriangulationWrapper<3>, 3>;
  template class MMeshExplicitGridFactory<MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshLeafIndexSet<const MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshGeometry<3, 3, const MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshGeometry<2, 3, const MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshGeometry<1, 3, const MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshGeometry<0, 3, const MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshEntity<0, 3, const MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshEntity<1, 3, const MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshEntity<2, 3, const MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshEntity<3, 3, const MMesh<TriangulationWrapper<3>, 3>>;

  template class MMeshInterfaceGrid<MMesh<TriangulationWrapper<3>, 3>>;
  template class MMeshInterfaceGridLeafIndexSet<const MMeshInterfaceGrid<MMesh<TriangulationWrapper<3>, 3>>>;
  template class MMeshInterfaceGridGeometry<2, 3, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<3>, 3>>>;
  template class MMeshInterfaceGridGeometry<1, 3, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<3>, 3>>>;
  template class MMeshInterfaceGridGeometry<0, 3, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<3>, 3>>>;
  template class MMeshInterfaceGridEntity<0, 2, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<3>, 3>>>;
  template class MMeshInterfaceGridEntity<1, 2, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<3>, 3>>>;
  template class MMeshInterfaceGridEntity<2, 2, const MMeshInterfaceGrid<MMesh<TriangulationWrapper<3>, 3>>>;

} // namespace Dune

#endif // DUNE_MMESH_MMESH_CC
