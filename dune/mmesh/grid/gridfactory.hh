// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_MMESH_GRID_GRIDFACTORY_HH
#define DUNE_MMESH_GRID_GRIDFACTORY_HH

// MMesh includes
#include "explicitgridfactory.hh"
#include "implicitgridfactory.hh"

namespace Dune
{

  /** \brief specialization of the GridFactory for MMesh
   *
   *  \ingroup GridFactory
   */

  //! Default grid factory for MMesh
  template< class HostGrid, int dim >
  class GridFactory< MMesh<HostGrid, dim> >
   : public MMeshExplicitGridFactory< MMesh<HostGrid, dim> > {};

} // end namespace Dune

#endif
