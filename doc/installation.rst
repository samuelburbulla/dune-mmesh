.. _installation:

************
Installation
************

In order to install and use Dune-MMesh you need a recent C++ compiler (at least C++17 compatible), Python (3.7 or later), CMake (3.13 or later) and Boost (1.66 or later).

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


From Source
-----------

You can install Dune-MMesh from source to get full access to the source code.
It also enables git support if you want to contribute.

1. Clone the Dune modules `dune-common <https://gitlab.dune-project.org/core/dune-common.git>`_,
`dune-geometry <https://gitlab.dune-project.org/core/dune-geometry.git>`_,
`dune-grid <https://gitlab.dune-project.org/core/dune-grid.git>`_,
`dune-istl <https://gitlab.dune-project.org/core/dune-istl.git>`_,
`dune-localfunctions <https://gitlab.dune-project.org/core/dune-localfunctions.git>`_,
`dune-fem <https://gitlab.dune-project.org/dune-fem/dune-fem.git>`_
and `dune-mmesh <https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh.git>`_.

.. code-block:: bash

  git clone https://gitlab.dune-project.org/core/dune-common.git
  git clone https://gitlab.dune-project.org/core/dune-geometry.git
  git clone https://gitlab.dune-project.org/core/dune-grid.git
  git clone https://gitlab.dune-project.org/core/dune-istl.git
  git clone https://gitlab.dune-project.org/core/dune-localfunctions.git
  git clone https://gitlab.dune-project.org/dune-fem/dune-fem.git
  git clone https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh.git

2. In a virtual environment (strongly recommended, see above) build the modules and install the python bindings.

.. code-block:: bash

  ./dune-common/bin/dunecontrol --opts=dune-mmesh/cmake/config.opts all
  ./dune-common/bin/dunecontrol --opts=dune-mmesh/cmake/config.opts make install_python

3. Configure the automatically generated dune-py module by calling

.. code-block:: bash

  ./dune-common/bin/setup-dunepy.py --opts=dune-mmesh/cmake/config.opts

Remark that this generated dune-py module is used to perform the just-in-time compilation that is used for
the python bindings of DUNE.
