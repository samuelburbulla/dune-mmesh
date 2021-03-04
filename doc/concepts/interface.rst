**************
Interface Grid
**************

Consider a domain :math:`\Omega \subset \mathbb{R}^n, n \in \{2,3\},` that includes a
:math:`(n-1)`-dimensional interface :math:`\Gamma \subset \Omega`. The interface might be curved and have junctions.

.. tikz:: A domain with a T-shaped interface.

  \draw (0,0) node[anchor=south west]{$\Omega$} rectangle (3,3);
  \draw[thick] (1,1) -- (1.5,1.5);
  \draw[thick] (1.5,1.5) -- (2,2) node[anchor=south west]{$\Gamma$};
  \draw[thick] (1.5,1.5) -- (2,1);

We assume the domain is triangulated conforming to the interface :math:`\Gamma`.

.. include:: grid.msh.rst

Let us denote this triangulation by :math:`\mathcal{T}`
and the set of facets - edges in 2d and faces in 3d - by :math:`\mathcal{F}`.

Due to conforming meshing, there exists a subset of facets :math:`\mathcal{F}_\Gamma \subset \mathcal{F}`
that belong to the interface :math:`\Gamma`.

These facets in :math:`\mathcal{F}_\Gamma` can also be interpreted as a triangulation of a surface.
We call this surface triangulation the `interface grid` and denote it by :math:`\mathcal{T}_\Gamma`.
