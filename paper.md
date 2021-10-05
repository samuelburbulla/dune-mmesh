---
title: 'Dune-MMesh: The DUNE Grid Module for Moving Interfaces'
tags:
  - Python
  - C++
  - DUNE
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

Dune-MMesh is an implementation of the DUNE `@BBD+21` grid interface that is tailored for numerical applications with possibly moving physical interfaces. The implementation based on CGAL `@CGAL` triangulations supports two and three dimensional meshes and can export a predefined set of facets as a separate interface grid. In spatial dimension two, arbitrary movement of vertices is enhanced with a remeshing algorithm that implements non-hierarchical adaptation procedures. We present a collection of examples based on the python bindings of the discretization module dune-fem `@DNK20` that demonstrate the versatile applicability of Dune-MMesh.

# Statement of need

In many technical applications, in particular in the field of fluid dynamics, comparably thin physical interfaces can have a large impact on the overall behaviour of a modeled system. For instance, interfaces occur as separating layer between fluid phases in multiphase flows, in fluid- structure interaction and fluid-solid phase change. Even fractures in porous media can be modeled by lower-dimensional surfaces. Oftentimes, these interfaces move over time and the processes become free-boundary value problems.
The grid implementation Dune-MMesh aims at providing numerical capabilities for grid based methods to model interface-driven processes within the DUNE framework. Essentially, it consists of two things:
1. A triangulation based on CGAL where a set of facets is considered as interface and
2. the possibility to re-mesh the triangulation when necessary.
These two ingredients enable many new possibilities within the DUNE framework. First, the representation of some grid facets as an interface makes Dune-MMesh a useful tool for the implementation of mixed-dimensional models. Second, the inevitable non-hierarchical adaptation
complements the existing grid implementations within the DUNE framework and allows for unprecedent flexibility of grid adaptation.

# Concepts

The concepts behind Dune-MMesh can be split into three main parts: the CGAL wrapper, the interface grid and the moving mesh capability.

## CGAL Wrapper

In its core, Dune-MMesh is a wrapper of CGAL Triangulations in $\mathbb{R}^d, d = 2, 3,$
that implements the Dune grid interface.
Therefore, it is essential to understand how CGAL triangulation objects are translated into Dune entities.

First of all, a CGAL triangulation is a set of simplicial cells and vertices.
Each cell gives access to its $d+1$ incident vertices and its $d+1$ adjacent cells.
Each vertex gives access to one of its incident cells.
The $d+1$ vertices are indexed with $0, 1, \dots, d$ in positive orientation being defined by the orientation of
the underlying Euclidian space $\mathbb{R}^d$.
The neighbors of a cell are also indexed with $0, 1, \dots, d$ in such a way
that the neighbor is opposite to the vertex with the same index.
Facets are not explicitly represented: a facet is given by the pair of a cell $c$
and an index $i$. Here, the facet $i$ of cell $c$ is the facet of $c$ that is
opposite to the vertex with index $i$. Remark that a facet has two implicit representations.
For $d=3$, edges are represented by triples of a cell $c$ and
two indices $i$ and $j$ that indicate the two vertices of the edge.

![CGAL representation of cells and differing Dune numbering in brackets.\label{fig:wrapper}](img/wrapper.png)

In order to match the Dune grid interface we have to follow the reference element numbering, cf. Firgure \ref{fig:wrapper}.
Fortunately, the vertex numbering of cells can be retained.
However, each facet $i$ of the CGAL representation corresponds to the codim-1 subentity $d-i$ in the Dune reference element.
For the representation of Dune intersections we can directly use CGAL's cell-index representation of facets
which is already equipped with an orientation.
With this reference mapping all geometry and sub-entity objects of the Dune grid interface can be specified.

Various iterators of CGAL triangulations can directly be used to construct the Dune grid range generators.
For instance, the element iterator coincides with the \texttt{{finite\_faces\_iterator}} or \texttt{{finite\_cells\_iterator}}.
Additional (non-standard Dune) iterators could be added easily, e.g. \texttt{{incidentElements}} or \texttt{{incidentVertices}} of a vertex.

The main large objects that have to be implemented are the index and id sets.
For this purpose, we define ids of entities as follows. At creation, each vertex is equipped with a unique integer id.
Each higher dimensional entity's id is defined by the sorted tuple of corresponding vertex ids.

As CGAL vertices and cells allow to append data (called: \emph{info}) to the objects, we can store and access the vertex ids directly within the vertex objects.
Entity indices are consecutively distributed at grid creation (or after adaptation) and also can be stored in the corresponding cell or vertex info.
For entities of codimensions different than $0$ and $d$, an id-index mapping is used.

The geometrical representation of entities that are not intrinsically CGAL entities (i.e., codimensions $1,...,d-1$) is made unique
by an ascending order of vertex ids. In additon, this prevents twists of intersections and we obtain a twist free grid implementation.

We extend the above described concepts of wrapping the CGAL triangulation to export a set of facets as interface grid.


## Interface Grid

Consider a domain $\Omega \subset \mathbb{R}^d, d \in \{2,3\},$ that includes a
$(d-1)$-dimensional interface $\Gamma \subset \Omega$, as depicted in Figure \ref{fig:triangulation}.
We assume the domain is triangulated conforming to the interface $\Gamma$.

![A domain with a T-shaped interface and an example for a conforming triangulation.\label{fig:triangulation}](img/triangulation.png)

Let us denote this triangulation by $\mathcal{T}$ and the set of facets by $\mathcal{F}$.
Due to conforming meshing, there exists a subset of facets $\mathcal{F}_\Gamma \subset \mathcal{F}$
that belong to the interface $\Gamma$.
Therefore, these facets in $\mathcal{F}_\Gamma$ can also be interpreted as a triangulation of a surface.
We call this surface triangulation the \emph{interface grid} and denote it by $\mathcal{T}_\Gamma$.

Dune-MMesh features a second implementation of the Dune grid interface that represents the interface triangulation $\mathcal{T}_\Gamma$.
Therefore, facets have to be marked as belonging to the interface - usually this is done when parsing a .msh file.

The interface grid can be used like any other Dune grid as it implements all necessary functionality.

A codim-0 entity of the interface grid is represented by a CGAL cell-index pair as used for the codim-1 entities of the wrapper implementation.
This representation is made unique by taking the representation where the cell has the lower index - which
is also considered to be the positive side of the facet.

All subentity objects can be generated by this representation using the right indexing of vertices.
The geometry representations and element ids are made unique by ascending order of vertex ids as it is done
in the full-dimensional wrapper implementation.

For iteration, CGAL's \texttt{{finite\_edges\_iterator}} or \texttt{{finite\_facets\_iterator}} is used skipping all facets not belonging to the interface. Intersections and neighbor relationships are obtained by CGAL's \texttt{{incident\_edges}} or \texttt{{incident\_facets}} iterators.
Index sets are implemented by mappings of vertex ids.

The interface grid also supports networks. For this purpose, the intersection iterator returns all common intersections with
adjacent cells. Note that this can be more than one for a single codim-1 subentity.
However, the intersection's outer normal is always independent of the neighbor entity, cf. Figure \ref{fig:junction}.

![Outer normals at junctions.\label{fig:junction}](img/junction.png)

Each bulk grid intersection can be identified belonging to the interface or not.
It is also possible to convert bulk intersections to interface grid elements and vice versa as the underlying representation is the same.
When converting an interface grid entity to a bulk intersection, Dune-MMesh returns the intersection as seen from the cell with the lower index.


## Moving Mesh

Most interface driven-problems have time-dependent interfaces $\Gamma = \Gamma(t)$.
Therefore, Dune-MMesh features capabilities of moving and remeshing in spatial dimension two.

### Moving Vertices

Dune-MMesh allows the movement of interface vertices (or all grid vertices) by a predescribed movement.

For this, we assume that movement is given by the shift of vertices.
This movement can be performed by simply changing the coordinates of the vertices.
Dune-MMesh provides the method \texttt{{moveInterface(shifts)}} that takes a vector of shift coordinates indexed by interface vertex indices as in Figure \ref{fig:moving}.

![Moving the interface.\label{fig:moving}](img/moving.png)
![Marking elements. Green: refine. Red: coarsen.\label{fig:mark}](img/mark.png)

A second method \texttt{{moveVertices(shifts)}} is available for moving all vertices of the triangulation that is indexed by bulk vertex indices.

Remark that moving vertices might lead to degeneration of the triangulation, i.e. cells can have non-positive volume.
To prevent that, Dune-MMesh is equipped with remeshing routines we describe in the following.


### Adaptation

Adaption in DUNE is hierarchical by definition. Whenever a grid element is supposed to be refined,
it is split into smaller cells belonging to a higher level of the grid hierarchy.
If all children in the highest refinement level of a grid element are supposed to be coarsened,
the children cells are put together to form a parent cell one level lower.

Hereby, the adaptation procedure is performed in two stages:
1. Mark: Grid elements are marked for coarsening or refinement.
2. Adapt: The elements are adapted due to their markers and discrete functions are restricted or prolongated.

In Dune-MMesh, due to the moving mesh, non-hierarchic adaptation is inavoidable.
However, we will try to follow the general DUNE approach of adaptation as good as possible.
For this reason, we similarily separate the adaptation into two stages.

__1. Mark__

In advance of moving, two methods are provided for marking elements in a convenient way.
First, the method \texttt{{ensureInterfaceMovement(shifts)}} (respectively \\ \texttt{{ensureVertexMovement(shifts)}}) can be called
to prepare Dune-MMesh for moving the vertices. The routine takes the vertex shifts as argument
and marks presumbly degenerate cells for coarsening, cf. Figure \ref{fig:mark}. Hence, they will be somehow removed during adaptation.


The second method available for marking elements is \texttt{{markElements()}}.
This method uses a default indicator that marks elements depending on their current geometrical properties.

This indicator considers primarily maximal and minimal edge length and aims at an objective edge length between $h_{max}$ and $h_{min}$.

- If an edge is longer than the maximum edge length $h_{max}$, the cell will be marked for refine.
- If an edge is shorter than the minimum edge length $h_{min}$, the cell will be marked for coarsening.

Additionally, if the ratio of longest to shortest edge is larger than 4, the cell is marked for coarsening.
The number 4 occurs from the fact that we we will use bisection and a triangle where two edges are longer then $h_{max}$
should not be splitted into smaller triangles where an edge is shorter than $h_{min}$.

Finally, a maximal radius ratio is taken into account to remove very ugly cells.
Always coarsening has priority before refinement because refinement would not remove ugly cells.

The minimal and maximal edge lengths $h_{max}$ and $h_{min}$ are
initialized automatically when constructing a mesh by determining the range of edge lengths occuring the grid.

Remark that \texttt{{markElements()}} also checks the elements of the interface grid.
Therefore, the interface will be refined and coarsened as well if edges of the interface get too long or too short.

The methods \texttt{{ensureInterfaceMovement(shifts)}} and \texttt{{markElements()}} are just convenience methods.
Instead, one can also use a proprietary procedure marking elements manually, or one can insert and remove vertices directly
using \texttt{{removeVertex(vertex)}} and \texttt{{refineEdge(element, edgeIndex)}}.


__2. Adapt__

After marking elements the \texttt{{adapt()}} routine performs the actual adaptation process.
The adaptation is performed by insertion and removal of points.

![Inserting and removing points.\label{fig:adapt}](img/adapt.png)
![Connected components.\label{fig:conncomp}](img/conncomp.png)

- In each element that is marked for refinement the center of the longest edge is inserted, i.e. refinement is done via bisection, cf. Figure \ref{fig:adapt}.

- In all elements marked for coarsening, one vertex is removed. Here, the vertex incident to the shortest edges of the cell is chosen, but we give priority on non-interface and non-boundary vertices.

When a vertex is removed, the resulting star-shaped hole is re-triangulated with respect to the interface.
Here, for the purpose of projection, we introduce \emph{connected components}, see Figure \ref{fig:conncomp}.
These are the minimal sets of cells from the triangulation before adaptation that
cover the same area as a set of cells in the triangulation afterwards.
The easiest representatives of these connected components are the incident cells when bisecting an edge and
the incident cells to a vertex that is removed.
Though, we have to combine overlapping sets of these representatives.

![Non-hierarchic projection with cut-set triangulation.\label{fig:projection}](img/projection.png)


For a conservative projection of discrete functions we compute a cut-set triangulation
which enables evalutation with agglomerated quadrature rules on triangles.
Here, we prolong from an old cell onto such a cut triangle and prolong onto the new cell, cf. Figure \ref{fig:projection}.
This whole projection is performed under the hood and just assumes that you use the callback adaptation in dune-fem.
We use a similar concept on the interface grid that enables projection of discrete functions on the interface.


## Trace and skeleton

Dune-MMesh exports both traces of bulk discrete functions on the interface and skeleton representations of interface discrete functions on bulk edges.

````
from dune.mmesh import trace, skeleton
tr = trace(bulkFunction)
sk = avg(skeleton(interfaceFunction, grid=geometryGridView))
````

The trace is a discrete function on the interface grid that evaluates the given bulk discrete function. It can be restricted to positive ('+') or negative side ('-') and it might be used in UFL forms.

Analogously, the skeleton function is a discrete function that returns the interface's discrete function values on interface bulk facets. Hereby, it has to be restricted formally though the value is independent of the side. The bulk grid can be passed optionally, but it is necessary if the grid view has been wrapped.


Please note that if a trace is used in surface integrals of forms on the interface grid, one must use already restricted trace objects to prevent confusion with the interface facet's restriction.
````
trp = trace(bulkFunction, restrictTo='+')
trm = trace(bulkFunction, restrictTo='-')
````


## Coupled solve

We provide two helper functions to solve bulk and interface schemes in a coupled way. Remark that for both methods the target discrete functions must be used in the coupling forms, too.

The first method uses an iterative solution strategy with a vector formulation of Aitken's fix point acceleration.

````
from dune.mmesh import iterativeSolve
iterativeSolve(schemes=(scheme, ischeme), targets=(sol, isol), callback=None, iter=100, tol=1e-8, f_tol=None, verbose=False)
````

The callback function is called every time before solving a scheme.
The iteration is performed until the residuum (measured in two norm between two iterates) is below the objective tolerance.

The second helper function solves bulk and interface scheme coupled monolithically.
A newton method is implemented assembling the underlying jacobian matrix.
Hereby, the coupling jacobian blocks are evaluated by finite difference on demand.
There are two implementations: A fast version with a C++ backend using UMFPACK and a slow python version relying on scipy.

````
from dune.mmesh import monolithicSolve
monolithicSolve(schemes=(scheme, ischeme), targets=(sol, isol), callback=None, iter=30, tol=1e-8, f_tol=1e-5, eps=1e-8, verbose=0, python=False)
````



# Examples

We implemented a few examples to display how Dune-MMesh can be used in different contexts. All examples can be found in `dune-mmesh/doc/examples` as IPython notebooks. Some examples for the creation of grid files can be found in `doc/examples/grids` which rely on gmsh `@GR09`.

![Finite volume moving mesh method to track a discontinuity `@CMR+18`\label{fig:fvmm}](img/fvmm.png)

![Mixed-dimensional model of poro-elasticity.\label{fig:poro}](img/poro.png)

![Two-phase Navier-Stokes equation `@GBK20`.\label{fig:navierstokes}](img/navierstokes.png)


# Acknowledgements

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project Number 327154368 - SFB 1313.

We thank all contributors that improved Dune-MMesh via the GitLab repository, especially Timo Koch.

# References
