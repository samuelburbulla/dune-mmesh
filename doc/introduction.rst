.. _introduction:

************
Introduction
************

In many technical applications, in particular in the field of fluid dynamics, comparably thin physical interfaces can have a large impact on the overall behavior of a modeled system.
For instance, interfaces occur as separating layer between fluid phases in multiphase flows, in fluid- structure interaction and fluid-solid phase change.
Even fractures in porous media can be modeled by lower-dimensional surfaces.
Oftentimes, these interfaces move over time and the processes become a kind of free-boundary value problems.

The grid implementation Dune-MMesh aims at providing numerical capabilities for grid based methods to model interface-driven processes within the DUNE framework.
Essentially, it consists of two things:
 1. A triangulation based on CGAL where a set of facets is considered as interface and
 2. the possibility to re-mesh the triangulation when necessary.

These two ingredients enable many new possibilities within the DUNE framework.
First, the representation of some grid facets as an interface makes Dune-MMesh a useful tool for the implementation of mixed-dimensional models.
Second, the inevitable non-hierarchical adaptation complements the existing grid implementations within the DUNE framework and allows for unprecedent flexibility of grid adaptation.

More details about the concepts behind Dune-MMesh are described in :ref:`concepts`.
The procedure to install and use Dune-MMesh is specified in :ref:`installation`.
You can find a collection of examples of what can be done with Dune-MMesh based on the discretization module Dune-Fem in :ref:`examples`.
The programming interface is described in the sections :ref:`pyreference` and :ref:`cppreference`.
