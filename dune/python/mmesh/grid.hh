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

    namespace MMIFGrid
    {
      //! register interface grid
      template< class Grid, class... options >
      void registerHierarchicalGrid ( pybind11::module module, pybind11::class_<Grid, options... > cls )
      {
        Dune::Python::registerHierarchicalGrid( module, cls );

        cls.def_property_readonly( "bulkGrid", [] ( const Grid &grid ) {
            return grid.getMMesh().leafGridView();
          },
          R"doc(
            Obtain bulk grid of the MMesh

            Returns:  bulk grid
          )doc" );
      }

    } // end namespace MMIFGrid

    namespace MMGrid
    {
      //! register bulk grid
      template< int d, class... options >
      void registerHierarchicalGrid ( pybind11::module module, pybind11::class_<Dune::MovingMesh<d>, options...> cls )
      {
        Dune::Python::registerHierarchicalGrid( module, cls );

        using Grid = Dune::MovingMesh< d >;

        cls.def_property_readonly( "interfaceHierarchicalGrid", [] ( const Grid &grid ) -> const auto & {
            return grid.interfaceGrid();
          }, pybind11::return_value_policy::reference, pybind11::keep_alive< 0, 1 >(),
          R"doc(
            Obtain interface hierarchical grid of the MMesh

            Returns:  interface hierarchical grid
          )doc" );

        using Element = typename Grid::template Codim< 0 >::Entity;
        using Vertex = typename Grid::template Codim< d >::Entity;
        using Intersection = typename Grid::Intersection;
        using InterfaceEntity = typename Grid::InterfaceEntity;
        using InterfaceGrid = typename Grid::InterfaceGrid;
        using InterfaceVertex = typename InterfaceGrid::template Codim< InterfaceGrid::dimension >::Entity;
        using FieldVector = Dune::FieldVector< double, d >;

        cls.def( "preAdapt", [] ( Grid &self ) {
          self.preAdapt();
        },
        R"doc(
          Prepare grid for adaption
        )doc" );

        cls.def( "ensureInterfaceMovement", [] ( Grid &self, const std::vector< FieldVector >& shifts ) {
          return self.ensureInterfaceMovement( shifts );
        },
        R"doc(
          Ensure the non-degeneration of the mesh after movement of the interface vertices
        )doc" );

        cls.def( "markElements", [] ( Grid &self ) {
          return self.markElements();
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

        cls.def( "moveVertices", [] ( Grid &self, const std::vector< FieldVector >& shifts ) {
          self.moveVertices( shifts );
        },
        R"doc(
          Move all vertices of the triangulation by the given movement
        )doc" );

        cls.def( "isInterface", [] ( Grid &self, const Intersection& intersection ) {
          return self.isInterface( intersection );
        },
        R"doc(
          Return if intersection is part of the interface
        )doc" );

        cls.def( "isInterface", [] ( Grid &self, const Vertex& vertex ) {
          return self.isInterface( vertex );
        },
        R"doc(
          Return if vertex is part of the interface
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

        cls.def( "addInterface", [] ( Grid &self, const Intersection& intersection ) {
          self.addInterface( intersection );
        },
        R"doc(
          Add the intersection to the set of interface edges
        )doc" );

        cls.def( "addInterface", [] ( Grid &self, const Intersection& intersection, const std::size_t marker ) {
          self.addInterface( intersection, marker );
        },
        R"doc(
          Add the intersection to the set of interface edges and mark it with marker
        )doc" );

        cls.def( "postAdapt", [] ( Grid &self ) {
          self.postAdapt();
        },
        R"doc(
          Remove adaption markers and connected components
        )doc" );

        cls.def( "isTip", [] ( Grid &self, const InterfaceVertex& interfaceVertex ) {
          return interfaceVertex.impl().isTip();
        },
        R"doc(
          Return if interface vertex is a tip
        )doc" );

        cls.def( "refineEdge", [] ( Grid &self, const Element& element, const std::size_t edgeIndex, const double where ) {
          return self.refineEdge(element, edgeIndex, where);
        },
        R"doc(
          Refine edge manually
        )doc" );

        cls.def( "refineEdge", [] ( Grid &self, const Element& element, const std::size_t edgeIndex ) {
          return self.refineEdge(element, edgeIndex, 0.5);
        },
        R"doc(
          Refine edge manually
        )doc" );

        cls.def( "removeVertex", [] ( Grid &self, const Vertex& vertex ) {
          return self.removeVertex(vertex);
        },
        R"doc(
          Remove vertex manually
        )doc" );

        cls.def( "removeVertex", [] ( Grid &self, const InterfaceVertex& vertex ) {
          return self.removeVertex(vertex);
        },
        R"doc(
          Remove interface vertex manually
        )doc" );

        cls.def( "insertVertexInCell", [] ( Grid &self, const FieldVector& position ) {
          return self.insertVertexInCell(position);
        },
        R"doc(
          Insert vertex in cell manually
        )doc" );

      }

    } // namespace MMGrid

  } // namespace Python

} // namespace Dune

#endif
