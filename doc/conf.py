# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
import subprocess
sys.path.insert(0, os.path.abspath('../python/'))

# -- Project information -----------------------------------------------------

project = 'dune-mmesh'
copyright = '2021, Samuel Burbulla'
author = 'Samuel Burbulla, Andreas Dedner, Maximilian Hörl, Christian Rohde'

# The full version, including alpha/beta/rc tags
release = 'master'


# -- General configuration ---------------------------------------------------

autodoc_mock_imports = ["dune.generator", "dune.fem"]

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
  "sphinx.ext.intersphinx",
  "sphinx.ext.napoleon",
  "sphinx.ext.autodoc",
  "sphinx.ext.mathjax",
  "sphinx.ext.viewcode",
  "sphinx_rtd_theme",
  "breathe",
  "nbsphinx",
  "sphinxcontrib.tikz",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

breathe_projects = { "dune-mmesh": "xml" }
breathe_default_project = 'dune-mmesh'

tikz_proc_suite = 'GhostScript'
tikz_tikzlibraries = 'calc'
tikz_resolution = 1000

nbsphinx_execute = 'never'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "sphinx_rtd_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# -- Options for Latex output -------------------------------------------------
latex_documents = [
    ('index', 'dune-mmesh.tex', 'Dune-MMesh: The DUNE Grid Module for Moving Interfaces', author, 'article')
]

latex_toplevel_sectioning = 'section'
latex_elements = {
    'tableofcontents': '',
    'preamble': r'''
\newcommand{\jump}[1]{[\mskip-5mu[ #1 ]\mskip-5mu]}
\newcommand{\avg}[1]{\{\mskip-5mu\{ #1 \}\mskip-5mu\}}
''',
    'printindex': '',
}

def generate_doxygen_xml(app):
    subprocess.call(["doxygen", "doxygen/Doxyfile"], cwd=app.confdir)

def setup(app):
    app.connect("builder-inited", generate_doxygen_xml)
