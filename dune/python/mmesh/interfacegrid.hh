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
        auto clsLeafView = insertClass< typename Grid::InterfaceGrid::LeafGridView >( module, "InterfaceGrid", GenerateTypeName( cls, "InterfaceGrid::LeafGridView" ) );
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
        Dune::Python::registerHierarchicalGrid( module, cls );
        registerInterfaceGrid( module, cls );

        using Grid = Dune::MovingMesh< d >;
        using Element = typename Grid::template Codim< 0 >::Entity;
        using Intersection = typename Grid::Intersection;
        using InterfaceEntity = typename Grid::InterfaceEntity;
        using FieldVector = Dune::FieldVector< double, d >;

        cls.def( "preAdapt", [] ( Grid &self ) {
          self.preAdapt();
        },
        R"doc(
          Prepare grid for adaption
        )doc" );

        cls.def( "ensureInterfaceMovement", [] ( Grid &self, const std::vector< FieldVector >& shifts ) {
          self.ensureInterfaceMovement( shifts );
        },
        R"doc(
          Ensure the non-degeneration of the mesh after movement of the interface vertices
        )doc" );

        cls.def( "markElements", [] ( Grid &self ) {
          self.markElements();
        },
        R"doc(
          Mark all elements in accordance to the default indicator
        )doc" );

        cls.def( "adapt", [] ( Grid &self ) {
          self.adapt();
        },
        R"doc(
          Adapt the grid
        )doc" );

        cls.def( "getConnectedComponent", [] ( Grid &self, const Element& element ) {
          return self.getConnectedComponent( element );
        },
        R"doc(
          Return the connected component of an entity
        )doc" );

        cls.def( "moveInterface", [] ( Grid &self, const std::vector< FieldVector >& shifts ) {
          self.moveInterface( shifts );
        },
        R"doc(
          Move the interface by the given movement for each interface vertex
        )doc" );

        cls.def( "isInterface", [] ( Grid &self, const Intersection& intersection ) {
          return self.isInterface( intersection );
        },
        R"doc(
          Return if intersection is part of the interface
        )doc" );

        cls.def( "asInterfaceEntity", [] ( Grid &self, const Intersection& intersection ) {
          return self.asInterfaceEntity( intersection );
        },
        R"doc(
          Return intersection as entity of the interface grid
        )doc" );

        cls.def( "asIntersection", [] ( Grid &self, const InterfaceEntity& interfaceEntity ) {
          return self.asIntersection( interfaceEntity );
        },
        R"doc(
          Return entity of the interface grid as (some) intersection of the MMesh
        )doc" );

        cls.def( "postAdapt", [] ( Grid &self ) {
          self.postAdapt();
        },
        R"doc(
          Remove adaption markers and connected components
        )doc" );
      }

    } // namespace MMGrid

  } // namespace Python

} // namespace Dune

#endif
