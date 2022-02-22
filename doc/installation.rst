.. _installation:

************
Installation
************

In order to install and use Dune-MMesh you need:
- C++ compiler (at least C++17 compatible, e.g. clang >= 5 or g++ >= 7)
- CMake (3.13 or later)
- pkg-config
- Boost (1.66 or later)
- OpenMPI (this will become optional)
- Python (3.7 or later)

There are two ways of installing Dune-MMesh.

Using Pip
---------

The easiest way to install Dune-MMesh is using pip and the package uploaded to `PyPI <https://pypi.org/project/dune-mmesh/>`_.

1. Activate a virtual environment (strongly recommended).

.. code-block:: bash

  python3 -m venv dune-env
  source dune-env/bin/activate

2. Download and build Dune-MMesh and its dependencies.

.. code-block:: bash

  pip install dune-mmesh

Note that this takes some time in order to compile all dependent Dune modules.

If you encounter problems with, e.g., Boost headers missing (e.g. on an M1 Mac),
make sure that the include paths can be found. For instance, use the export

.. code-block:: bash

  export CXXFLAGS="-I/opt/homebrew/Cellar/boost/1.66.0/include/"

before installing Dune-MMesh.


From Source
-----------

You can install Dune-MMesh from source to get full access to the source code.
It also enables git support if you want to contribute.

1. Clone the Dune modules `dune-common <https://gitlab.dune-project.org/core/dune-common.git>`_,
`dune-geometry <https://gitlab.dune-project.org/core/dune-geometry.git>`_,
`dune-grid <https://gitlab.dune-project.org/core/dune-grid.git>`_,
`dune-istl <https://gitlab.dune-project.org/core/dune-istl.git>`_,
`dune-localfunctions <https://gitlab.dune-project.org/core/dune-localfunctions.git>`_,
`dune-alugrid <https://gitlab.dune-project.org/extensions/dune-alugrid.git>`_
`dune-fem <https://gitlab.dune-project.org/dune-fem/dune-fem.git>`_
and `dune-mmesh <https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh.git>`_.

.. code-block:: bash

  git clone https://gitlab.dune-project.org/core/dune-common.git
  git clone https://gitlab.dune-project.org/core/dune-geometry.git
  git clone https://gitlab.dune-project.org/core/dune-grid.git
  git clone https://gitlab.dune-project.org/core/dune-istl.git
  git clone https://gitlab.dune-project.org/core/dune-localfunctions.git
  git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git
  git clone https://gitlab.dune-project.org/dune-fem/dune-fem.git
  git clone https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh.git

2. Build the modules. This will create an internal virtual environment and install the python bindings.

.. code-block:: bash

  ./dune-common/bin/dunecontrol --opts=dune-mmesh/cmake/config.opts all

3. Activate the DUNE internal virtual environment.

.. code-block:: bash

  source ./dune-common/build-cmake/dune-env/bin/activate

Remark that a `dune-py` module will be generated automatically that is necessary to perform the just-in-time compilation of DUNE python modules.
