************
Known Issues
************

There are a few known issues that either just haven't been implemented so far
or need further clarification:

* Adaptation in spatial dimension three and for parallel computations
* Topology change of interface (e.g. merging)


The remeshing feature is not (yet) supported in spatial dimension three because the removal of a vertex is not offered by the underlying CGAL triangulation class.
In fact, it could appear that the region formed by its adjacent tetrahedrons is an instance of the untetrahedralizable Schönhardt’s polyhedron.
In this case, the removal of the vertex might be impossible without rebuilding the whole triangulation.
The Dune-MMesh grid implementation supports solvers that are parallelized with MPI, but the adaptation is not yet generalized to the parallel case.
For parallelization one can also rely on multi-threading features like the threaded Galerkin operator in Dune-Fem.

We do not handle topology changes of interface (e.g. merging of two interfaces), yet.
A simple intersection algorithm for a tip crossing a second interface has been implemented, but other kinds of interaction of two interfaces need further clarification on how to handle this interaction.
