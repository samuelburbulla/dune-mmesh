******************************************
The DUNE grid module for moving interfaces
******************************************

Dune-MMesh is an implementation of the `DUNE <https://www.dune-project.org>`_ [BBD+21]_ grid interface that
is tailored for numerical applications with moving physical interfaces.
The implementation based on `CGAL <https://www.cgal.org>`_ [TCP20]_ triangulations supports two and three dimensional meshes
and can export a predefined set of facets as a separate interface grid.
In spatial dimension two, arbitrary movement of vertices is enhanced with a remeshing algorithm
that implements non-hierarchical adaptation procedures.
We present a collection of :ref:`examples` based on the python bindings of the discretization module `dune-fem <https://www.dune-project.org/sphinx/dune-fem/>`_ [DNK+20]_.

.. toctree::
  :hidden:
  :maxdepth: 1
  :caption: First steps

  introduction
  installation

.. toctree::
  :hidden:
  :maxdepth: 1
  :caption: Next steps

  concepts
  examples

.. toctree::
  :hidden:
  :maxdepth: 1
  :caption: API Reference

  pyreference
  cppreference

.. toctree::
  :hidden:
  :maxdepth: 1
  :caption: Appendix

  issues
  acknowledgements
  contribute


.. [BBD+21] P. Bastian, M. Blatt, A. Dedner, N.-A. Dreier, C. Engwer, R. Fritze, C. Gräser, C. Grüninger, D. Kempf, R. Klöfkorn, M. Ohlberger, O. Sander. The DUNE framework: Basic concepts and recent developments. Computers & Mathematics with Applications 81, 2021, pp. 75-112.

.. [TCP20] The CGAL Project. CGAL User and Reference Manual. CGAL Editorial Board, 5.2 edition, 2020.

.. [DNK+20] A. Dedner, M. Nolte, and R. Klöfkorn. Python Bindings for the DUNE-FEM module. Zenodoo, 2020, DOI 10.5281/zenodo.3706994.
