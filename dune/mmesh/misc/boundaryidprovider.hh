#ifndef DUNE_MMESH_MISC_BOUNDARYIDPROVIDER_HH
#define DUNE_MMESH_MISC_BOUNDARYIDPROVIDER_HH

#if HAVE_DUNE_FEM

#include <dune/fem/misc/boundaryidprovider.hh>

namespace Dune
{

  namespace Fem
  {

    template <class HostGrid, int dim>
    struct BoundaryIdProvider< MMesh<HostGrid,dim> >
    {
      typedef MMesh<HostGrid,dim> GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.impl().boundaryId()) : 0);
      }
    };

    template <class MMesh>
    struct BoundaryIdProvider< MMeshInterfaceGrid<MMesh> >
    {
      typedef MMeshInterfaceGrid<MMesh> GridType;

      template< class Intersection >
      static int boundaryId ( const Intersection &intersection )
      {
        return (intersection.boundary() ? (intersection.impl().boundaryId()) : 0);
      }
    };

  }  // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_FEM

#endif // #ifndef DUNE_MMESH_MISC_TWISTUTILITY_HH
