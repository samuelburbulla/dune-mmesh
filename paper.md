---
title: 'Dune-MMesh: The Dune Grid Module for Moving Interfaces'
tags:
  - Python
  - C++
  - Dune
  - partial differential equations
  - mixed-dimensional
  - moving mesh
authors:
  - name: Samuel Burbulla^[corresponding author]
    orcid: 0000-0002-2566-9777
    affiliation: 1
  - name: Andreas Dedner
    affiliation: 2
  - name: Maximilian Hörl
    affiliation: 1
  - name: Christian Rohde
    affiliation: 1
affiliations:
 - name: University of Stuttgart
   index: 1
 - name: University of Warwick
   index: 2
date: 05 October 2021
bibliography: paper.bib

---

# Summary

Dune-MMesh is an implementation of the Dune [@BBD+21] grid interface that is tailored for numerical applications with possibly moving physical interfaces. The implementation based on CGAL triangulations [@CGAL] supports two and three dimensional meshes and can export a predefined set of facets as a separate interface grid. In spatial dimension two, arbitrary movement of vertices is enhanced with a remeshing algorithm that implements non-hierarchical adaptation procedures. Various examples based on the python bindings of the discretization module dune-fem [@DNK20] have been implemented that demonstrate the versatile applicability of Dune-MMesh.

# Statement of need

In many technical applications, in particular in the field of fluid dynamics, comparably thin physical interfaces can have a large impact on the overall behavior of a modeled system. Interfaces occur as separating layer between fluid phases in multiphase flows, in fluid-structure interaction, fluid-solid phase change and even fractures are modeled by lower-dimensional surfaces.

The grid implementation Dune-MMesh aims at providing numerical capabilities for grid based methods to model interface-driven processes within the Dune framework. Essentially, it consists of two things: A triangulation based on CGAL where a set of facets is considered as interface and the possibility to re-mesh the triangulation when necessary.
The representation of some grid facets as an interface makes Dune-MMesh a useful tool for the implementation of mixed-dimensional models.
The inevitable non-hierarchical adaptation complements the existing grid implementations within the Dune framework and allows for unprecedented flexibility of grid adaptation.

# CGAL Wrapper

In its core, Dune-MMesh is a wrapper of CGAL Triangulations in $\mathbb{R}^d, d = 2, 3,$ that implements the Dune grid interface.
A CGAL triangulation is a set of simplicial cells and vertices where each cell gives access to its $d+1$ incident vertices and cells.
Facets are not explicitly represented: a facet is given by the pair of a cell $c$ and an index $i$ and has two implicit representations.
For $d=3$, edges are represented by triples of a cell $c$ and two indices $i$ and $j$ that indicate the two vertices of the edge.

![CGAL representation of cells and differing Dune numbering in brackets.\label{fig:wrapper}](img/wrapper.png)

In order to match the Dune grid reference cell numbering we apply an index mapping, cf. Figure \ref{fig:wrapper}.
Dune intersections can directly be represented by CGAL's cell-index representations of facets which are already equipped with an orientation.
The index and id sets of the Dune grid interface are realized by consecutive numbering of cells and vertices.
Various iterators of CGAL triangulations can directly be used to construct the Dune grid range generators.
Additional (non-standard Dune) iterators have been added, e.g. iterating over incident cells of a vertex.

# Interface Grid

Consider a domain $\Omega \subset \mathbb{R}^d, d \in \{2,3\},$ that includes a
$(d-1)$-dimensional interface $\Gamma \subset \Omega$, as depicted in Figure \ref{fig:triangulation}.
We assume the domain is triangulated conforming to the interface $\Gamma$.

![A domain with a T-shaped interface and an example for a conforming triangulation.\label{fig:triangulation}](img/triangulation.png){ width=80% }

Dune-MMesh features a second implementation of the Dune grid interface that represents the interface triangulation.
Here, a codim-0 entity of the interface grid is represented by a CGAL cell-index pair.
The interface grid also supports networks, cf. Figure \ref{fig:junction}, and it is possible to convert bulk intersections to interface grid cells and vice versa.

![Outer normals at junctions.\label{fig:junction}](img/junction.png){ width=30% }

# Moving Mesh

Most interface driven-problems have time-dependent interfaces $\Gamma = \Gamma(t)$.
Therefore, Dune-MMesh features capabilities of moving and remeshing in spatial dimension two.

### Moving Vertices

We assume that movement is given by a shift of interface vertices (or all grid vertices).

![Moving the interface.\label{fig:moving}](img/moving.png){ width=30% } ![Marking cells. Green: refine. Red: coarsen.\label{fig:mark}](img/mark.png){ width=30% }

To prevent degeneration of the triangulation, i.e. cells have non-positive volume, Dune-MMesh is equipped with remeshing routines.

### Adaptation

Adaption in Dune is usually hierarchical by definition and the adaptation procedure is performed in two stages:

1. Mark: Grid cells are marked for coarsening or refinement.
2. Adapt: The cells are modified due to their markers and discrete functions are restricted or prolongated.

In Dune-MMesh, due to the moving mesh, non-hierarchic adaptation is inavoidable.
However, we will try to follow the general Dune approach and separate the adaptation into two stages.

__1. Mark__

Dune-MMesh provides utility functions to mark cells either in expectation of a movement of vertices or regarding to their current geometrical properties. However, one can also use a proprietary procedure marking cells manually, or one can insert and remove vertices directly.

__2. Adapt__

After marking cells the `adapt` routine performs the actual adaptation process.
The adaptation is performed by insertion and removal of points.

![Inserting and removing points.\label{fig:adapt}](img/adapt.png){ width=30% } ![Connected components.\label{fig:conncomp}](img/conncomp.png){ width=30% }

In each cell that is marked for refinement we bisect the longest edge, cf. Figure \ref{fig:adapt}.
In all cells marked for coarsening, the least important vertex is removed.
When a vertex is removed, the resulting star-shaped hole is re-triangulated with respect to the interface.

For the purpose of projection, we introduce \emph{connected components}, see Figure \ref{fig:conncomp},
and implement a generalized callback adaptation in dune-fem.

![Non-hierarchic projection with cut-set triangulation.\label{fig:projection}](img/projection.png)

A conservative projection of discrete functions can be performed by intermediate
prolongation and restriction on the cut-set cells, cf. Figure \ref{fig:projection}.
We use a similar concept on the interface grid that enables projection of discrete functions on the interface.


# Trace and skeleton

Dune-MMesh exports both traces of bulk discrete functions on the interface and skeleton representations of interface discrete functions on bulk edges.

The trace is a discrete function on the interface grid that evaluates a given bulk discrete function.
It can be restricted to both sides of the interface and might be used in UFL forms.

Analogously, the skeleton function is a discrete function that returns the interface's discrete function values on interface bulk facets.

Both `trace` and `skeleton` can be used to couple bulk and interface problems.
Such couplings occur, e.g., in mixed-dimensional PDEs.

# Coupled solve

We provide two helper functions to solve bulk and interface schemes in a coupled way.

The first method `iterativeSolve` uses an iterative solution strategy which alternately solves both schemes until the two norm between two iterates is below an objective tolerance.

The second helper function `monolithicSolve` solves bulk and interface scheme coupled monolithically.
A newton method is implemented assembling the underlying jacobian matrix where the coupling jacobian blocks are evaluated by finite differences.

# Examples

We implemented a few examples to display how Dune-MMesh can be used in different contexts.
All examples can be found in [`dune-mmesh/doc/examples`](https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh/-/tree/master/doc/examples) as IPython notebooks.
Some numerical results of these examples are visualized in Figure \ref{fig:fvmm}, Figure \ref{fig:poro} and Figure \ref{fig:navierstokes}.

![Finite volume moving mesh method to track a discontinuity [@CMR+18]\label{fig:fvmm}](img/fvmm.png)

![Mixed-dimensional model of poro-elasticity with a t-shaped fracture.\label{fig:poro}](img/poro.png)

![Two-phase Navier-Stokes equation [@GBK20].\label{fig:navierstokes}](img/navierstokes.png)

# Acknowledgements

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project Number 327154368 - SFB 1313.

We thank all contributors that improved Dune-MMesh via the GitLab repository, especially Timo Koch.

# References
