The `dune-mmesh` module
=======================

The `dune-mmesh` module is an implementation of the `dune-grid` interface that wraps CGAL Delaunay triangulations in 2D and 3D.


Installation
------------

`dune-mmesh` requires the DUNE core modules, version 2.4 or later.

Additionally, you'll need `CGAL >= 4.12` installed on your system.


Examples
--------

See the `examples` folder for some grid creations.
You can either use a `.msh` file, a `.dgf` file or the `CGAL` triangulation interface directly.

Short example using the GmshReader:
````
    using Grid = Dune::MovingMesh< dim >;

    using GridFactory = Dune::GmshGridFactory< Grid >;
    GridFactory gridFactory( "grid.msh" );

    Grid& grid = *gridFactory.grid();
````
