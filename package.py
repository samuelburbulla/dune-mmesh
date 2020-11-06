# Meta data
name="dune-mmesh"
version="2.8"
author="Samuel Burbulla"
author_email="samuel.burbulla@mathematik.uni-stuttgart.de"
description="MMesh is a grid implementation based on CGAL triangulations."
url="https://gitlab.dune-project.org/samuel.burbulla/dune-mmesh"

# DUNE dependencies
dune_dependencies=['dune-grid']

# Package dependencies
install_requires=[]

# Module libaries that have to be compiled (without the _ prefix)
modules=['mmesh']

# Files to include in the source package
manifest='''\
graft cmake
graft doc
graft dune
graft examples
graft lib
graft python
graft scripts
include CMakeLists.txt
include dune-mmesh.pc.in
include dune.module
include LICENSE.md
include README.md
'''
