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
 - name: University of Stuttgart, Germany
   index: 1
 - name: University of Warwick, UK
   index: 2
date: 05 October 2021
bibliography: paper.bib

---

# Summary

In several physical and environmental processes that concern multiphase flows, biological systems, and geophysical phenomena, important physical processes occur along thin physical interfaces. These processes include effects that may alter the interface's position or topology over time creating a moving interface, which complicates traditional modeling techniques. Moving interface problems thus require advanced numerical tools with specific treatment of the interface and the simultaneous ability to implement complex physical effects, which this work seeks to create solutions for.

# Statement of Need

In this work, we present Dune-MMesh that is tailored for numerical applications with moving physical interfaces. Dune-MMesh is an implementation of the well-developed Dune [@BBD+21] grid interface and is well-suited for the numerical discretization of partial differential equations (PDEs). The package wraps two and three dimensional CGAL triangulations [@CGAL] in high-level objects like intersections of grid entities, index and id sets and geometry transformations and exports a predefined set of facets as a separate interface grid.
In two dimensions, the arbitrary movement of vertices is enhanced with a re-meshing algorithm that implements non-hierarchical adaptation procedures. Besides the adaptation of the triangulation, Dune-MMesh provides the necessary data structures to adapt discrete functions defined on the bulk grid or the interface. This adaptation approach complements existing grid implementations within the Dune framework that strictly rely on hierarchical adaptation.
Various examples in Python have been implemented based on the discretization module Dune-Fem [@DNK20] that demonstrate the versatile applicability of Dune-MMesh. Due to the ability to handle custom PDEs in their weak from written in Unified Form Language (UFL) and the mesh adaptation capabilities, we believe Dune-MMesh provides a useful tool for solving mixed-dimensional PDEs on moving interfaces that arise from various fields of modelling.

# CGAL Wrapper

In its core, Dune-MMesh is a wrapper of CGAL Triangulations in $\mathbb{R}^d, d = 2, 3,$ that implements the Dune grid interface.
A CGAL triangulation is a set of simplicial cells and vertices where each cell gives access to its $d+1$ incident vertices and cells.
Facets are not explicitly represented: a facet is given by the pair of a cell $c$ and an index $i$ and has two implicit representations.
For $d=3$, edges are represented by triples of a cell $c$ and two indices $i$ and $j$ that indicate the two vertices of the edge.

![CGAL representation of cells and differing Dune numbering in brackets. The vertex numbering is maintained, facets are renumbered, and the edges of tetrahedrons are equipped with indices according to the Dune reference element numbering. \label{fig:wrapper}](img/wrapper.png)

In order to match the Dune grid reference cell numbering we apply an index mapping, cf. Figure \ref{fig:wrapper}.
Here, the edges of tetrahedrons are equipped with indices according to the Dune reference element numbering.
Dune intersections, i.e., intersections of mesh entities of codimension 0 with a neighboring element or with the domain boundary, can directly be represented by CGAL's cell-index representations of facets which are already equipped with an orientation.
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
Therefore, Dune-MMesh features capabilities of moving and re-meshing in spatial dimension two.
Here we follow the approach of moving the interface edges and adapt the mesh next to the interface.

### Moving Vertices

We assume that movement is given by a shift of interface vertices (or all grid vertices), cf. Figure \ref{fig:movmark} (left).

![Left: Moving the interface is performed by shifting vertices. The blue shift vectors tranform the gray two-dimensional triangulation into the black one. Right: Marking cells for refinement (green) or coarsening (red).\label{fig:movmark}](img/movmark.png){ width=80% }

To prevent degeneration of the triangulation, i.e. cells have non-positive volume, Dune-MMesh is equipped with re-meshing routines that will be described in the subsequent.

### Adaptation

Adaptation in Dune is usually hierarchical by definition and the adaptation procedure is performed in two stages:

1. Mark: Grid cells are marked for coarsening or refinement.
2. Adapt: The cells are modified due to their markers and discrete functions are restricted or prolongated.

In Dune-MMesh, due to the moving mesh, non-hierarchic adaptation is unavoidable.
However, we will try to follow the general Dune approach and separate the adaptation into two stages.

__Stage 1: Mark__

Dune-MMesh provides utility functions to mark cells either in expectation of a movement of vertices or regarding to their current geometrical properties, cf. Figure \ref{fig:movmark} (right).
For instance, when moving the interface would cause a cell to get a negative volume, we mark this cell for coarsening (marked red in Figure \ref{fig:movmark}). Similarly, we use the edge length as indicator for coarsening or refinement (marked green).
However, one can also use a proprietary procedure marking cells manually, or one can insert and remove vertices directly.

__Stage 2: Adapt__

After marking cells an adapt routine performs the actual adaptation process.
The adaptation is performed by insertion and removal of points.

![Left: Inserting and removing points. Right: Connected components.\label{fig:adaptconn}](img/adaptconn.png){ width=80% }

In each cell that is marked for refinement we bisect the longest edge, cf. Figure \ref{fig:adaptconn} (left).
In all cells marked for coarsening, the least important vertex is removed.
When a vertex is removed, the resulting star-shaped hole is re-triangulated with respect to the interface.

For the purpose of projection, we introduce \emph{connected components}, see Figure \ref{fig:adaptconn} (right),
and implement a generalized callback adaptation in Dune-Fem.

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
Such couplings occur, for example, in mixed-dimensional PDEs.

# Coupled solve

We provide two helper functions to solve bulk and interface schemes in a coupled way.

The first method `iterativeSolve` uses an iterative solution strategy which alternately solves both schemes until the two norm between two iterates is below an objective tolerance.

The second helper function `monolithicSolve` solves bulk and interface scheme coupled monolithically.
A newton method is implemented assembling the underlying jacobian matrix where the coupling jacobian blocks are evaluated by finite differences.

# Examples

We implemented a few examples to display how Dune-MMesh can be used in different contexts.
All examples can be found in [`dune-mmesh/doc/examples`](https://github.com/samuelburbulla/dune-mmesh/tree/master/doc/examples) as IPython notebooks.
Some numerical results of these examples are visualized in Figure \ref{fig:fvmm}, Figure \ref{fig:poro} and Figure \ref{fig:navierstokes}.

![Finite volume moving mesh method to track a discontinuity [@CMR+18]\label{fig:fvmm}](img/fvmm.png)

![Mixed-dimensional model of poro-elasticity with a T-shaped fracture.\label{fig:poro}](img/poro.png)

![Two-phase Navier-Stokes equation [@GBK20].\label{fig:navierstokes}](img/navierstokes.png)

# Acknowledgements

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project Number 327154368 - SFB 1313.

We thank all contributors that improved Dune-MMesh via the GitLab repository, especially Timo Koch.

# References
