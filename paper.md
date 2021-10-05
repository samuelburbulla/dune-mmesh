---
title: 'Gala: A Python package for galactic dynamics'
tags:
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Adrian M. Price-Whelan^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-0872-7098
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID^[co-first author] # note this makes a footnote saying 'co-first author'
    affiliation: 2
  - name: Author with no affiliation^[corresponding author]
    affiliation: 3
affiliations:
 - name: Lyman Spitzer, Jr. Fellow, Princeton University
   index: 1
 - name: Institution Name
   index: 2
 - name: Independent Researcher
   index: 3
date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

Dune-MMesh is an implementation of the DUNE grid interface that is tailored for numerical applications with possibly moving physical interfaces. The implementation based on CGAL triangulations supports two and three dimensional meshes and can export a predefined set of facets as a separate interface grid. In spatial dimension two, arbitrary movement of vertices is enhanced with a remeshing algorithm that implements non-hierarchical adaptation procedures. We present a collection of examples based on the python bindings of the discretization module dune-fem that demonstrate the versatile applicability of Dune-MMesh.

# Statement of need

In many technical applications, in particular in the field of fluid dynamics, comparably thin physical interfaces can have a large impact on the overall behaviour of a modeled system. For instance, interfaces occur as separating layer between fluid phases in multiphase flows, in fluid- structure interaction and fluid-solid phase change. Even fractures in porous media can be modeled by lower-dimensional surfaces. Oftentimes, these interfaces move over time and the processes become free-boundary value problems.
The grid implementation Dune-MMesh aims at providing numerical capabilities for grid based methods to model interface-driven processes within the DUNE framework. Essentially, it consists of two things:
1. A triangulation based on CGAL where a set of facets is considered as interface and
2. the possibility to re-mesh the triangulation when necessary.
These two ingredients enable many new possibilities within the DUNE framework. First, the representation of some grid facets as an interface makes Dune-MMesh a useful tool for the implementation of mixed-dimensional models. Second, the inevitable non-hierarchical adaptation
complements the existing grid implementations within the DUNE framework and allows for unprecedent flexibility of grid adaptation.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

Funded by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) - Project Number 327154368 - SFB 1313.

We thank all contributors that improved Dune-MMesh via the GitLab repository, especially Timo Koch.

# References
