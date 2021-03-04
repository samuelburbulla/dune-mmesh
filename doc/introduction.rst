************
Introduction
************

In many technical applications - especially in fluid dynamics - thin physical interfaces have a large impact on
the behavior of a system. For instance, interfaces occur as separating layer between fluid phases in multiphase flows,
in fluid-structure interaction and fluid-solid phase change, and even fractures in porous media can be modeled
by lower-dimensional surfaces.
Oftentimes, these interfaces move over time and the processes become free-boundary value problems.

The grid implementation dune-mmesh aims at providing numerical capabilities for grid based methods to model interface-driven processes.
Basically, it consists of a triangulation that is able to remesh and a set of facets that are considered as interface.

More details about the background of dune-mmesh can be found in `about`_, the general concepts are described in `concepts`_.
A few examples are shown in `examples_`.

.. toctree::
  :maxdepth: 1

  introduction/about
