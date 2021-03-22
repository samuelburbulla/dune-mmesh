.. _moving:

***********
Moving Mesh
***********

Most interface driven-problems have time-dependent interfaces :math:`\Gamma = \Gamma(t)`.
Therefore, Dune-MMesh features capabilities of moving and remeshing in spatial dimension two.

We assume that movement is given by shifts of vertices.
Then, this movement can be performed by just changing the coordinates of vertices.
Therefore, Dune-MMesh implements the method

.. code-block:: cpp

  moveInterface(shifts)

that takes a vector of shift coordinates indexed by interface vertex indices.
The shifts are simply used to update the vertex position coordinates.


.. tikz:: Moving the interface.

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


There is a second method :code:`moveVertices` that can be used to move all vertices of the triangulation.
This method is indexed by bulk vertex indices.

Remark that moving might lead to degeneration of the triangulation, i.e. cells have non-positive volume.
To prevent that, Dune-MMesh is equipped with remeshing routines.
Remeshing is performed in two stages.
First, elements are marked for coarsening or refinement.
Second, the mesh is adapted.

Mark
----

In advance of moving, the method

.. code-block:: cpp

  ensureInterfaceMovement(shifts)

(respectively :code:`ensureVertexMovement`) can be called to prepare Dune-MMesh for moving the vertices.
The routine checks for presumbly degenerate cells and marks them for coarsening. Hence, they will be removed during adaptation.

The second method available for marking elements is

.. code-block:: cpp

  markElements()

which uses a default indicator that marks elements
for coarsening or refinement depending on their current geometrical properties.

.. tikz:: Marking elements. Here, green for refinement, red for coarsening.

  \coordinate (A) at (-0.7,-0.5);
  \coordinate (B) at (1,1);
  \coordinate (C) at (1.5,1.5);
  \coordinate (D) at (0.2,1.2);
  \coordinate (E) at (0.1,3);
  \coordinate (F) at (1.1,-0.6);
  \coordinate (G) at (1.8,0.4);
  \coordinate (H) at (0.5,0.1);
  \draw[fill=black!40!green] (A) -- (D) -- (B);
  \draw (D) -- (E) -- (B);
  \draw[fill=black!40!green] (B) -- (E) -- (C);
  \draw[fill=black!10!red] (A) -- (H) -- (B);
  \draw (A) -- (F) -- (B);
  \draw (H) -- (F) -- (G);
  \draw (C) -- (G) -- (B);
  \draw[very thick] (A) -- (B) -- (C);

  \draw[->, black!20!blue, thick] (A) -- +(0.1,-0.4);
  \draw[->, black!20!blue, thick] (B) -- +(0.2,-0.4);
  \draw[->, black!20!blue, thick] (C) -- +(0.5,-0.3);



This indicator considers primarily maximal and minimal edge length.
If an edge is longer (shorter) than a maximum (minimum) edge length :math:`h_{max}` (:math:`h_{min}`),
the cell will be marked for refine (coarsening).

Additionally, if the ratio of longest to shortest edge is larger than 4, the cell is marked for coarsening.
The number 4 occurs from the fact that we we will use bisection and a triangle where two edges are longer then :math:`h_{max}`
should not be splitted into smaller triangles where an edge is shorter than :math:`h_{min}`.

Finally, a maximal radius ratio is taken into account to remove very ugly cells.
Always coarsening has priority before refinement.

The minimal and maximal edge lengths :math:`h_{max}` and :math:`h_{min}` are
initialized automatically when constructing a mesh by determining the range of edge lengths occuring the grid.

The `markElements()` routine also checks all elements of the interface grid.
Therefore, the interface will be refined and coarsened as well if edges of the interface get too long or too short.

Adapt
-----

After marking elements the

.. code-block:: cpp

 adapt()

routine performs the actual adaptation process.

.. tikz:: Inserting and removing points.

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

  \draw[->, black!20!blue, thick] (A) -- +(0.1,-0.4);
  \draw[->, black!20!blue, thick] (B) -- +(0.2,-0.4);
  \draw[->, black!20!blue, thick] (C) -- +(0.5,-0.3);

  \draw[fill, black!40!green] (I) circle (2pt);
  \draw[fill, black!40!green] (J) circle (2pt);
  \draw[fill, black!40!green] (K) circle (2pt);
  \draw[fill, black!10!red] (H) circle (2pt);


The adaptation is performed by insertion and removal of points.
In each element that is marked for refinement the center of the longest edge is interserted,
i.e. refinement is done via bisection.
In all elements marked for coarsening, one vertex is removed. Here, the vertex incident to the
shortest edges of the cell is chosen, but we give priority on non-interface and non-boundary vertices.
When a vertex is removed, the resulting whole is retriangulated with respect to the interface.

For the purpose of projection we introduce *connected components*.
These are defined as sets of cells from the triangulation before adaptation that
cover the same area as a set of cells in the triangulation afterwards.
The easiest representatives of these connected components are the incident cells when bisecting an edge and
the incident cells to a vertex that is to be removed.
Though, we have to combine overlapping sets of these representatives.

.. tikz:: Connected components.

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

  \draw[->, black!20!blue, thick] (A) -- +(0.1,-0.4);
  \draw[->, black!20!blue, thick] (B) -- +(0.2,-0.4);
  \draw[->, black!20!blue, thick] (C) -- +(0.5,-0.3);


For a conservative projection of discrete functions we compute a cut-set triangulation
which enables evalutation with agglomerated quadrature rules on triangles.
Here, we prolong from an old cell onto such a cut triangle and prolong onto the new cell.
This whole projection is performed under the hood and just assumes that you use the callback adaptation in dune-fem.
We use a similar concept on the interface grid that enables projection of discrete functions on the interface.
