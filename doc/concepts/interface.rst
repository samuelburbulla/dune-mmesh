**************
Interface Grid
**************

Consider a domain :math:`\Omega \subset \mathbb{R}^n, n \in \{2,3\},` that includes a
:math:`(n-1)-dimensional` interface :math:`\Gamma \subset \Omega`.

.. tikz:: Domain with interface

  \draw (0,0) node[anchor=south west]{$\Omega$} rectangle (3,3);
  \draw[thick] (1,1) -- (1.5,1.5);
  \draw[thick] (1.5,1.5) -- (2,2) node[anchor=south west]{$\Gamma$};
  \draw[thick] (1.5,1.5) -- (2,1);
