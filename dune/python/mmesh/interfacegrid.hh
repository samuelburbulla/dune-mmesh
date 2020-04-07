// -*- tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PYTHON_MMESH_INTERFACEGRID_HH
#define DUNE_PYTHON_MMESH_INTERFACEGRID_HH

#include <array>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <type_traits>

#include <dune/common/hybridutilities.hh>
#include <dune/common/iteratorrange.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include <dune/python/common/mpihelper.hh>
#include <dune/python/common/typeregistry.hh>

#include <dune/python/grid/capabilities.hh>
#include <dune/python/grid/enums.hh>
#include <dune/python/grid/factory.hh>
#include <dune/python/grid/gridview.hh>
#include <dune/python/grid/idset.hh>

#include <dune/python/pybind11/functional.h>
#include <dune/python/pybind11/numpy.h>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

namespace Dune
{

  namespace Python
  {

    namespace MMGrid
    {

      template< class Grid, class... options >
      void registerInterfaceGrid ( pybind11::module module, pybind11::class_< Grid, options... > cls )
      {
        auto clsLeafView = insertClass< typename Grid::InterfaceGrid::LeafGridView >( module, "InterfaceGrid", GenerateTypeName( cls, "LeafGridView" ) );
        if( clsLeafView.second )
          registerGridView( module, clsLeafView.first );
        cls.def_property_readonly( "interfaceGrid", [] ( const Grid &grid ) {
            return grid.interfaceGrid().leafGridView();
          }, pybind11::return_value_policy::reference, pybind11::keep_alive< 0, 1 >(),
          R"doc(
            Obtain interface grid of the MMesh

            Returns:  interface grid
          )doc" );
      }

      template< int d, class... options >
      void registerHierarchicalGrid ( pybind11::module module, pybind11::class_<Dune::MovingMesh<d>, options...> cls )
      {
        typedef Dune::MovingMesh<d> Grid;
        Dune::Python::registerHierarchicalGrid( module, cls );
        registerInterfaceGrid ( module, cls );
      }

    } // namespace MMGrid

  } // namespace Python

} // namespace Dune

#endif
