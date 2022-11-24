[![Release](https://img.shields.io/github/v/release/samuelburbulla/dune-mmesh.svg)](https://github.com/samuelburbulla/dune-mmesh/releases)
[![Test](https://github.com/samuelburbulla/dune-mmesh/actions/workflows/test.yml/badge.svg)](https://github.com/samuelburbulla/dune-mmesh/actions/)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03959/status.svg)](https://doi.org/10.21105/joss.03959)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# The Dune-MMesh Module

In several physical and environmental processes that concern multiphase flows, biological systems, and geophysical phenomena, important physical processes occur along thin physical interfaces. These processes include effects that may alter the interface's position or topology over time creating a moving interface, which complicates traditional modeling techniques. Moving interface problems thus require advanced numerical tools with specific treatment of the interface and the simultaneous ability to implement complex physical effects.

Dune-MMesh is tailored for numerical applications with moving physical interfaces. It is an implementation of the well-developed [Dune](https://www.dune-project.org) grid interface and is well-suited for the numerical discretization of partial differential equations. The package wraps two and three dimensional [CGAL](https://www.cgal.org) triangulations in high-level objects like intersections of grid entities, index and id sets and geometry transformations and exports a predefined set of facets as a separate interface grid.
In two dimensions, the arbitrary movement of vertices is enhanced with a re-meshing algorithm that implements non-hierarchical adaptation procedures. Besides the adaptation of the triangulation, Dune-MMesh provides the necessary data structures to adapt discrete functions defined on the bulk grid or the interface. This adaptation approach complements existing grid implementations within the Dune framework that strictly rely on hierarchical adaptation.
Various examples in Python have been implemented based on the discretization module [dune-fem](https://www.dune-project.org/sphinx/dune-fem/) that demonstrate the versatile applicability of Dune-MMesh. Due to the ability to handle custom PDEs in their weak from written in Unified Form Language (UFL) and the mesh adaptation capabilities, we believe Dune-MMesh provides a useful tool for solving mixed-dimensional PDEs on moving interfaces that arise from various fields of modelling.

You can find the full documentation of Dune-MMesh at [dune-mmesh.readthedocs.io](https://dune-mmesh.readthedocs.io).

## Installation

Note that Dune-MMesh has a list of dependencies: C++ compiler, CMake, Python3 + pip (+ venv), pkg-config, Boost, OpenMPI, SuiteSparse, Gmsh.

We strongly recommend using a virtual environment:
````
python3 -m venv dune-env
source dune-env/bin/activate
````

Install the Dune-MMesh package using pip:
````
pip install dune-mmesh
````
This will take some time in order to compile all dependent Dune modules.

Now, you should be able to execute Dune-MMesh's python code. For instance:
````
git clone https://github.com/samuelburbulla/dune-mmesh.git
cd dune-mmesh/doc/examples
python coupling.py
````

For more details on the installation procedure we refer to [Installation](https://dune-mmesh.readthedocs.io/en/latest/installation.html).


### Docker image

The easiest starting point is to use Docker with a preconfigured setup.
Using the pre-built Docker container you can simply run:

````
docker run -it ghcr.io/samuelburbulla/dune-mmesh:master
````

This will open an interactive shell in the Dune-MMesh's examples directory.


## Examples

You can find a collection of examples of how to use Dune-MMesh on our [Examples](https://dune-mmesh.readthedocs.io/en/latest/examples.html) page.

## Testing
You can test your installation of Dune-MMesh by running the python tests
````
python -m dune.mmesh test
````
Further tests of the C++ backend can be performed with a [source build](https://dune-mmesh.readthedocs.io/en/latest/installation.html#from-source) executing `make build_test` and `make test` in the build directory.

## Contribution

Contributions are highly welcome. If you want to contribute, please use [GitHub](https://github.com/samuelburbulla/dune-mmesh/)
or our [GitLab](https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh) repository to report an issue or open a merge/pull request.

## License
Dune-MMesh is licensed under the terms and conditions of the GNU General Public License (GPL) version 3 or - at your option - any later version.
