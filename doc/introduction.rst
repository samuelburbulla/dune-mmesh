.. _introduction:

************
Introduction
************

In many technical applications - especially in fluid dynamics - thin physical interfaces have a large impact on
the behavior of a system. For instance, interfaces occur as separating layer between fluid phases in multiphase flows,
in fluid-structure interaction and fluid-solid phase change, and even fractures in porous media can be modeled
by lower-dimensional surfaces. Oftentimes, these interfaces move over time and the processes become free-boundary value problems.

The grid implementation Dune-MMesh aims at providing numerical capabilities for grid based methods to model interface-driven processes
within the DUNE framework. Essentially, it consists of a triangulation with a set of facets that are considered as interface
and the possibility to remesh.

More details about the concepts behind Dune-MMesh are described in :ref:`concepts`.
A few examples of what you can do with Dune-MMesh are presented in :ref:`examples`.
The whole programming interface is described in the sections :ref:`pyreference` and :ref:`cppreference`.
