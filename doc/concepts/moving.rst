.. _moving:

***********
Moving Mesh
***********

Most interface driven-problems have time-dependent interfaces :math:`\Gamma = \Gamma(t)`.
Therefore, Dune-MMesh features capabilities of moving and remeshing in spatial dimension two.

.. note::
  The remeshing feature is not (yet) supported in spatial dimension three because the removal of a vertex is not
  offered by the underlying CGAL Triangulation_3 class. In fact, it could appear that the region formed by its
  adjacent tetrahedrons is an instance of the untetrahedralizable Schönhardt's polyhedron. In this case, the
  removal of the vertex might be impossible without rebuilding the whole triangulation.


Moving Vertices
***************

Dune-MMesh allows the movement of interface vertices (or all grid vertices) by a predescribed movement.

For this, we assume that movement is given by the shift of vertices.
This movement can be performed by simply changing the coordinates of the vertices.
Dune-MMesh provides the method :code:`moveInterface(shifts)` that takes a vector of shift coordinates indexed by interface vertex indices.

.. tikz:: Moving the interface.
  :xscale: 20

  \coordinate (A) at (-0.7,-0.5);
  \coordinate (B) at (1,1);
  \coordinate (C) at (1.5,1.5);
  \coordinate (D) at (0.2,1.2);
  \coordinate (E) at (0.1,3);
  \coordinate (F) at (1.1,-0.6);
  \coordinate (G) at (1.8,0.4);
  \coordinate (H) at (0.5,0.1);
  \coordinate (A2) at ($(A)+(0.1,-0.4)$);
  \coordinate (B2) at ($(B)+(0.2,-0.4)$);
  \coordinate (C2) at ($(C)+(0.5,-0.3)$);
  \coordinate (I) at ($(A2)!0.5!(B2)$);
  \coordinate (J) at ($(B2)!0.5!(E)$);
  \coordinate (K) at ($(C2)!0.5!(E)$);

  \draw[white!80!black] (A) -- (D) -- (B);
  \draw[white!80!black] (D) -- (E) -- (B);
  \draw[white!80!black] (B) -- (E) -- (C);
  \draw[white!80!black] (A) -- (F) -- (B);
  \draw[white!80!black] (C) -- (J) -- (K);
  \draw[white!80!black] (C) -- (G) -- (B);
  \draw[white!80!black, very thick] (A) -- (B) -- (C);

  \draw[->, black!20!blue, thick] (A) -- +(0.1,-0.4);
  \draw[->, black!20!blue, thick] (B) -- +(0.2,-0.4);
  \draw[->, black!20!blue, thick] (C) -- +(0.5,-0.3);

  \draw (A2) -- (D) -- (B2);
  \draw (D) -- (E) -- (B2);
  \draw (B2) -- (E) -- (C2);
  \draw (A2) -- (F) -- (B2);
  \draw (I) -- (F) -- (G);
  \draw (I) -- (D);
  \draw (J) -- (D);
  \draw (C2) -- (J) -- (K);
  \draw (C2) -- (G) -- (B2);
  \draw[very thick] (A2) -- (B2) -- (C2);


A second method :code:`moveVertices(shifts)` is available for moving all vertices of the triangulation that is indexed by bulk vertex indices.


Remark that moving vertices might lead to degeneration of the triangulation, i.e. cells can have non-positive volume.
To prevent that, Dune-MMesh is equipped with remeshing routines we describe in the following.


Adaptation
**********

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

1. Mark
-------

In advance of moving, two methods are provided for marking elements in a convenient way.

.. tikz:: Marking elements. Green: refine. Red: coarsen.
  :xscale: 20

  \coordinate (A) at (-0.7,-0.5);
  \coordinate (B) at (1,1);
  \coordinate (C) at (1.5,1.5);
  \coordinate (D) at (0.2,1.2);
  \coordinate (E) at (0.1,3);
  \coordinate (F) at (1.1,-0.6);
  \coordinate (G) at (1.8,0.4);
  \coordinate (H) at (0.5,0.1);
  \draw[fill={rgb:red,50;green,205;blue,50}] (A) -- (D) -- (B);
  \draw (D) -- (E) -- (B);
  \draw[fill={rgb:red,50;green,205;blue,50}] (B) -- (E) -- (C);
  \draw[fill=black!10!red] (A) -- (H) -- (B);
  \draw (A) -- (F) -- (B);
  \draw (H) -- (F) -- (G);
  \draw (C) -- (G) -- (B);
  \draw[very thick] (A) -- (B) -- (C);

  \draw[->, black, thick] (A) -- +(0.1,-0.4);
  \draw[->, black, thick] (B) -- +(0.2,-0.4);
  \draw[->, black, thick] (C) -- +(0.5,-0.3);


First, the method :code:`ensureInterfaceMovement(shifts)` (respectively :code:`ensureVertexMovement(shifts)`) can be called
to prepare Dune-MMesh for moving the vertices. The routine takes the vertex shifts as argument
and marks presumbly degenerate cells for coarsening. Hence, they will be somehow removed during adaptation.

The second method available for marking elements is :code:`markElements()`.
This method uses a default indicator that marks elements depending on their current geometrical properties.

This indicator considers primarily maximal and minimal edge length and aims at an objective edge length between :math:`h_{max}` and :math:`h_{min}`.

 - If an edge is longer than the maximum edge length :math:`h_{max}`, the cell will be marked for refine.
 - If an edge is shorter than the minimum edge length :math:`h_{min}`, the cell will be marked for coarsening.

Additionally, if the ratio of longest to shortest edge is larger than 4, the cell is marked for coarsening.
The number 4 occurs from the fact that we we will use bisection and a triangle where two edges are longer then :math:`h_{max}`
should not be splitted into smaller triangles where an edge is shorter than :math:`h_{min}`.

Finally, a maximal radius ratio is taken into account to remove very ugly cells.
Always coarsening has priority before refinement because refinement would not remove ugly cells.

The minimal and maximal edge lengths :math:`h_{max}` and :math:`h_{min}` are
initialized automatically when constructing a mesh by determining the range of edge lengths occuring the grid.

Remark that :code:`markElements()` also checks the elements of the interface grid.
Therefore, the interface will be refined and coarsened as well if edges of the interface get too long or too short.

.. note:: The methods :code:`ensureInterfaceMovement(shifts)` and :code:`markElements()` are just convenience methods.
  Instead, one can also use a proprietary procedure marking elements manually, or one can insert and remove vertices directly
  using :code:`removeVertex(vertex)` and :code:`refineEdge(element, edgeIndex)`.


2. Adapt
--------

After marking elements the :code:`adapt()` routine performs the actual adaptation process.
The adaptation is performed by insertion and removal of points.

.. tikz:: Inserting and removing points.
  :xscale: 20

  \coordinate (A) at (-0.7,-0.5);
  \coordinate (B) at (1,1);
  \coordinate (C) at (1.5,1.5);
  \coordinate (D) at (0.2,1.2);
  \coordinate (E) at (0.1,3);
  \coordinate (F) at (1.1,-0.6);
  \coordinate (G) at (1.8,0.4);
  \coordinate (H) at (0.5,0.1);
  \coordinate (I) at ($(A)!0.5!(B)$);
  \coordinate (J) at ($(B)!0.5!(E)$);
  \coordinate (K) at ($(C)!0.5!(E)$);

  \draw (A) -- (D) -- (B);
  \draw (D) -- (E) -- (B);
  \draw (B) -- (E) -- (C);
  \draw (A) -- (F) -- (B);
  \draw (I) -- (F) -- (G);
  \draw (I) -- (D);
  \draw (J) -- (D);
  \draw (C) -- (J) -- (K);
  \draw (C) -- (G) -- (B);
  \draw[very thick] (A) -- (B) -- (C);

  \draw[->, black, thick] (A) -- +(0.1,-0.4);
  \draw[->, black, thick] (B) -- +(0.2,-0.4);
  \draw[->, black, thick] (C) -- +(0.5,-0.3);

  \draw[fill={rgb:red,50;green,205;blue,50}] (I) circle (2pt);
  \draw[fill={rgb:red,50;green,205;blue,50}] (J) circle (2pt);
  \draw[fill={rgb:red,50;green,205;blue,50}] (K) circle (2pt);
  \draw[fill=black!10!red] (H) circle (2pt);


- In each element that is marked for refinement the center of the longest edge is interserted, i.e. refinement is done via bisection.
- In all elements marked for coarsening, one vertex is removed. Here, the vertex incident to the shortest edges of the cell is chosen, but we give priority on non-interface and non-boundary vertices.

.. tikz:: Connected components.
  :xscale: 20

  \coordinate (A) at (-0.7,-0.5);
  \coordinate (B) at (1,1);
  \coordinate (C) at (1.5,1.5);
  \coordinate (D) at (0.2,1.2);
  \coordinate (E) at (0.1,3);
  \coordinate (F) at (1.1,-0.6);
  \coordinate (G) at (1.8,0.4);
  \coordinate (H) at (0.5,0.1);
  \draw[fill=yellow] (D) -- (E) -- (B);
  \draw[fill=yellow] (B) -- (E) -- (C);
  \draw[fill=white!70!blue] (A) -- (D) -- (B);
  \draw[fill=white!70!blue] (A) -- (F) -- (B);
  \draw (A) -- (H) -- (B);
  \draw (G) -- (F) -- (H);
  \draw (A) -- (F) -- (B);
  \draw (C) -- (G) -- (B);
  \draw[very thick] (A) -- (B) -- (C);

  \draw[->, black, thick] (A) -- +(0.1,-0.4);
  \draw[->, black, thick] (B) -- +(0.2,-0.4);
  \draw[->, black, thick] (C) -- +(0.5,-0.3);


When a vertex is removed, the resulting star-shaped hole is re-triangulated with respect to the interface.
Here, for the purpose of projection, we introduce *connected components*.
These are the minimal sets of cells from the triangulation before adaptation that
cover the same area as a set of cells in the triangulation afterwards.
The easiest representatives of these connected components are the incident cells when bisecting an edge and
the incident cells to a vertex that is removed.
Though, we have to combine overlapping sets of these representatives.

.. tikz:: Non-hierarchic projection with cut-set triangulation.
  :xscale: 75

  \coordinate (A) at (-0.7,-0.5);
  \coordinate (B) at (1,1);
  \coordinate (D) at (0.2,1.2);
  \coordinate (F) at (1.1,-0.6);
  \coordinate (H) at (0.5,0.1);
  \coordinate (I) at ($(A)!0.5!(B)$);
  \coordinate (P) at (0.383,0.0383);
  \draw[fill=white!70!red] (A) -- (H) -- (B);
  \fill[white!10!red] (I) -- (H) -- (P);
  \draw (A) -- (D) -- (B);
  \draw (A) -- (F) -- (B);
  \draw (A) -- (H) -- (B);
  \draw (F) -- (H);
  \draw (A) -- (F) -- (B);
  \draw[very thick] (A) -- (B);

  \draw[->] (1.6,0.3) -- (2,0.3) node[midway, above] {\tiny prolong};

  \coordinate (A2) at (-0.7+3,-0.5);
  \coordinate (B2) at (1+3,1);
  \coordinate (C2) at (1.5+3,1.5);
  \coordinate (D2) at (0.2+3,1.2);
  \coordinate (F2) at (1.1+3,-0.6);
  \coordinate (H2) at (0.5+3,0.1);
  \coordinate (I2) at ($(A2)!0.5!(B2)$);
  \coordinate (P2) at (0.383+3,0.0383);
  \draw (A2) -- (D2) -- (B2);
  \draw (A2) -- (F2) -- (B2);
  \fill[white!10!red] (I2) -- (H2) -- (P2);
  \draw[dashed] (I2) -- (F2);
  \draw[dashed] (I2) -- (D2);
  \draw[dashed] (A2) -- (B2);
  \draw[dashed] (A2) -- (H2) -- (B2);
  \draw[dashed] (F2) -- (H2);
  \draw[dashed] (I2) -- (H2);


  \draw[->] (4.7,0.3) -- (5.1,0.3) node[midway, above] {\tiny restrict};

  \coordinate (A3) at (-0.7+6,-0.5);
  \coordinate (B3) at (1+6,1);
  \coordinate (C3) at (1.5+6,1.5);
  \coordinate (D3) at (0.2+6,1.2);
  \coordinate (F3) at (1.1+6,-0.6);
  \coordinate (H3) at (0.5+6,0.1);
  \coordinate (I3) at ($(A3)!0.5!(B3)$);
  \coordinate (P3) at (0.383+6,0.0383);
  \fill[white!70!red] (I3) -- (F3) -- (B3);
  \fill[white!10!red] (I3) -- (H3) -- (P3);
  \draw (A3) -- (D3) -- (B3);
  \draw (A3) -- (F3) -- (B3);
  \draw (I3) -- (F3);
  \draw (I3) -- (D3);
  \draw[very thick] (A3) -- (B3);


For a conservative projection of discrete functions we compute a cut-set triangulation
which enables evalutation with agglomerated quadrature rules on triangles.
Here, we prolong from an old cell onto such a cut triangle and prolong onto the new cell.
This whole projection is performed under the hood and just assumes that you use the callback adaptation in dune-fem.
We use a similar concept on the interface grid that enables projection of discrete functions on the interface.
