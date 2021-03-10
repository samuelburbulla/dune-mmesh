.. _wrapper:

************
CGAL Wrapper
************

In its core, Dune-MMesh is a wrapper of CGAL (Delaunay) Triangulations in :math:`\mathbb{R}^d, d = 2, 3,`
that implements the Dune grid interface.
Therefore, it is essential to understand how CGAL triangulation objects are translated into Dune entities.

First of all, a CGAL triangulation is a set of simplicial cells and vertices.
Each cell gives access to its :math:`d+1` incident vertices and its :math:`d+1` adjacent cells.
Each vertex gives access to one of its incident cells.
The :math:`d+1` vertices are indexed with :math:`0, 1, \dots, d` in positive orientation being defined by the orientation of
the underlying Euclidian space :math:`\mathbb{R}^d`.
The neighbors of a cell are also indexed with :math:`0, 1, \dots, d` in such a way
that the neighbor is opposite to the vertex with the same index.
Facets are not explicitly represented: a facet is given by the pair of a cell :math:`c`
and an index :math:`i`. Here, the facet :math:`i` of cell :math:`c` is the facet of :math:`c` that is
opposite to the vertex with index :math:`i`. Remark that a facet has two implicit representations.
For :math:`d=3`, edges are represented by triples of a cell :math:`c` and
two indices :math:`i` and :math:`j` that indicate the two vertices of the edge.

.. tikz:: CGAL representation of cells.

  \tikzset{vertex/.style={circle, fill=white, draw, inner sep=1pt}}
  \tikzset{edge/.style={midway, sloped, circle, fill=white, inner sep=1pt}}

  \draw[thick] (0,0) node[vertex] {\tiny 0}
    -- (2,0) node[vertex] {\tiny 1} node[edge] {\tiny 2}
    -- (0,2) node[vertex] {\tiny 2} node[edge] {\tiny 0}
    -- (0,0) node[edge] {\tiny 1};

  \draw[thick] (4.5,0) node[vertex] {\tiny 0}
    -- (7,0) node[vertex] {\tiny 1}
    -- (5.5,1) node[vertex] {\tiny 2}
    -- (4.5,0)
    -- (5,2.5) node[vertex] {\tiny 3}
    -- (7,0) -- (5,2.5) -- (5.5,1);
  \node at (5.5,0.4) {\tiny 3};
  \node at (5.05,1.2) {\tiny 1};
  \node at (5.8,1.2) {\tiny 0};
  \node at (4.9,0.15) {\tiny 2};

In order to match the Dune grid interface we have to follow the reference element numbering.
Fortunately, the vertex numbering of cells can be retained.
However, each facet :math:`i` of the CGAL representation corresponds to the codim-1 subentity :math:`d-i` in the Dune reference element.
For the representation of Dune intersections we can directly use CGAL's cell-index representation of facets
which is already equipped with an orientation.
With this reference mapping all geometry and sub-entity objects of the Dune grid interface can be specified.

..  tikz:: Dune representation of codim-0-entities.

   \tikzset{vertex/.style={circle, fill=white, draw, inner sep=1pt}}
   \tikzset{edge/.style={midway, sloped, circle, fill=white, inner sep=1pt}}

   \draw[thick] (0,0) node[vertex] {\tiny 0}
     -- (2,0) node[vertex] {\tiny 1} node[edge] {\tiny 0}
     -- (0,2) node[vertex] {\tiny 2} node[edge] {\tiny 2}
     -- (0,0) node[edge] {\tiny 1};

   \draw[thick] (4.5,0) node[vertex] {\tiny 0}
     -- (7,0) node[vertex] {\tiny 1}
     -- (5.5,1) node[vertex] {\tiny 2}
     -- (4.5,0)
     -- (5,2.5) node[vertex] {\tiny 3}
     -- (7,0) -- (5,2.5) -- (5.5,1);
   \node at (5.5,0.4) {\tiny 0};
   \node at (5.05,1.2) {\tiny 2};
   \node at (5.8,1.2) {\tiny 3};
   \node at (4.9,0.15) {\tiny 1};


Various iterators of CGAL triangulations can directly by used to construct the Dune grid range generators.

The remaining large objects that have to be implemented are the index and id sets.
For this purpose, we define ids of entities as follows. At creation, each vertex is equipped with a unique integer id.
Each higher dimensional entity's id is defined by the sorted tuple of corresponding vertex ids.

As CGAL vertices and cells allow to append data (called *info*) to the objects, we can store and access the vertex ids directly within the vertex objects.
Indices are consecutively distributed at grid creation and stored in the corresponding info, too.

The geometrical representation of each entity is made unique by using the vertices sorted after their id.
This prevents twists of intersections and we obtain a twist free grid implementation.
