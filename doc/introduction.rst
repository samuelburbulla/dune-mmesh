.. _introduction:

************
Introduction
************
In several physical and environmental processes that concern multiphase flows, biological systems, and geophysical phenomena, important physical processes occur along thin physical interfaces. These processes include effects that may alter the interface's position or topology over time creating a moving interface, which complicates traditional modelling techniques. Moving interface problems thus require advanced numerical tools with specific treatment of the interface and the simultaneous ability to implement complex physical effects.

Dune-MMesh is tailored for numerical applications with moving physical interfaces. It is an implementation of the well-developed `Dune <https://www.dune-project.org>`_ grid interface and is well-suited for the numerical discretization of partial differential equations. The package wraps two and three dimensional `CGAL <https://www.cgal.org>`_ triangulations in high-level objects like intersections of grid entities, index and id sets and geometry transformations and exports a predefined set of facets as a separate interface grid.
In two dimensions, the arbitrary movement of vertices is enhanced with a re-meshing algorithm that implements non-hierarchical adaptation procedures. Besides the adaptation of the triangulation, Dune-MMesh provides the necessary data structures to adapt discrete functions defined on the bulk grid or the interface. This adaptation approach complements existing grid implementations within the Dune framework that strictly rely on hierarchical adaptation.
Various examples in Python have been implemented based on the discretization module `Dune-Fem <https://www.dune-project.org/sphinx/dune-fem/>`_ that demonstrate the versatile applicability of Dune-MMesh. Due to the ability to handle custom PDEs in their weak from written in Unified Form Language (UFL) and the mesh adaptation capabilities, we believe Dune-MMesh provides a useful tool for solving mixed-dimensional PDEs on moving interfaces that arise from various fields of modelling.


More details about the concepts behind Dune-MMesh are described in :ref:`concepts`.
The procedure to install and use Dune-MMesh is specified in :ref:`installation`.
You can find a collection of examples of what can be done with Dune-MMesh based on the discretization module Dune-Fem in :ref:`examples`.
The programming interface is described in the sections :ref:`pyreference` and :ref:`cppreference`.
