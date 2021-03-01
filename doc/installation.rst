************
Installation
************

The following dependencies are needed for dune-mmesh:

* C++ compiler (at least C++17 compatible)
* Python (3.4 or later)
* `CGAL <https://www.cgal.org>`_ (5.1 or later)


Using Pip
=========

The easiest way to install dune-mmesh is using pip.

We strongly recommend the usage of a virtual environment.
In some folder setup a virtual environment and activate it

.. code-block:: bash

  python3 -m venv dune-env
  source dune-env/bin/activate

Then download and build dune-mmesh and its dependecies

.. code-block:: bash

  pip install dune-mmesh

Note that this takes some time in order to compile all dependent Dune libraries.


From Source
===========

You can install dune-mmesh from source to get full access to the source code.
It also enables git support if you want to contribute.

* Required Dune modules (2.8 or later)

  * `dune-common <https://gitlab.dune-project.org/core/dune-common.git>`_
  * `dune-geometry <https://gitlab.dune-project.org/core/dune-geometry.git>`_
  * `dune-grid <https://gitlab.dune-project.org/core/dune-grid.git>`_

* Recommended Dune modules (2.8 or later)

  * `dune-fem <https://gitlab.dune-project.org/dune-fem/dune-fem.git>`_

Read the instructions on how to `build Dune with Python support`_ which also
links to general instructions on how to `build Dune modules`_. Note that
the Python bindings require some additional CMake flags to be set as
described on the pages linked above.

.. _build Dune modules: https://dune-project.org/doc/installation
.. _build Dune with Python support: https://dune-project.org/doc/pythonbindings
