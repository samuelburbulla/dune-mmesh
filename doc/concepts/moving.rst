***********
Moving Mesh
***********

The remeshing feature of dune-mmesh is able to insert and remove arbitrary vertices
at any time and assists in projecting unknowns.

A default adaptation strategy is provided that adapts the triangulation
in preparation of vertex movement. It uses an indicator that is
defined over several criteria as edge length ratio, radius ratio,
edge length and distance to the interface.

Per default, remeshing is performed by retriangulation of holes after
removal of vertices and the bisection of edges.

For conservative projection of discrete functions, connected components of
cells of the old triangulation are constructed that cover the same area as a
set of cells in the new one. The weighted average using the exact intersection
volumes of cells can be used for conservative projection of cell-wise defined values.

If some edges that belong to the interface have to be refined or coarsened,
a similar concept defined on the lower-dimensional triangulation is used.
