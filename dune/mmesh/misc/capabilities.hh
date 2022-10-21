#ifndef DUNE_MMESH_MISC_CAPABILITIES_HH
#define DUNE_MMESH_MISC_CAPABILITIES_HH

#if HAVE_DUNE_FEM

#include <dune/fem/misc/capabilities.hh>

namespace Dune
{

  namespace Capabilities
  {

    template< class HostGrid, int dim >
    struct hasHierarchicIndexSet< MMesh< HostGrid, dim > >
    {
      static const bool v = false;
    };

  } // namespace Capabilities


  namespace Fem
  {

    namespace Capabilities
    {

      template< class HostGrid, int dim >
      struct supportsCallbackAdaptation< MMesh< HostGrid, dim > >
      {
        static const bool v = true;
      };

      template< class HostGrid, int dim >
      struct isLocallyAdaptive< MMesh< HostGrid, dim > >
      {
        static const bool v = true;
      };

      template< class MMesh >
      struct supportsCallbackAdaptation< MMeshInterfaceGrid< MMesh > >
      {
        static const bool v = true;
      };

      template< class MMesh >
      struct isLocallyAdaptive< MMeshInterfaceGrid< MMesh > >
      {
        static const bool v = true;
      };

      template< class MMesh >
      struct isMMesh< MMeshInterfaceGrid< MMesh > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_DUNE_FEM

#endif // #if DUNE_MMESH_MISC_CAPABILITIES_HH
