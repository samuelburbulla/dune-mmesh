"""Solve a simple laplace problem."""

from ufl import TrialFunction, TestFunction, SpatialCoordinate, FacetNormal, FacetArea, inner, dot, grad, div, dx, dS, ds, jump, avg, sqrt, sin, pi
from dune.grid import reader
from dune.mmesh import mmesh
from dune.fem.space import dglagrange
from dune.fem.function import integrate
from dune.ufl import Constant
from dune.fem.scheme import galerkin

from dune.mmesh.test.grids import line, vertical

import sys
from time import time
from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()

def algorithm(grid, name):
  print(grid.size(0), "cells,", grid.ghostSize(0), "ghost")
  space = dglagrange(grid, order=1, storage="istl")
  x = SpatialCoordinate(space)
  uh = space.interpolate(0, name="uh")

  u = TrialFunction(space)
  v = TestFunction(space)

  exact = sin(pi*x[0]) * sin(pi*x[1])
  beta = Constant(3.01, name="beta")
  n = FacetNormal(space)
  h = avg(FacetArea(space))
  hBnd = FacetArea(space)

  a = inner(grad(u), grad(v)) * dx
  a += beta / h * inner(jump(u), jump(v)) * dS
  a -= dot(dot(avg(grad(u)), n("+")), jump(v)) * dS
  a -= dot(dot(avg(grad(v)), n("+")), jump(u)) * dS
  a += beta / hBnd * inner(u - exact, v) * ds
  a -= dot(dot(grad(u), n), v) * ds
  a -= dot(dot(grad(v), n), u - exact) * ds

  if name == "interface":
    b = -exact.dx(0).dx(0) * v * dx
  else:
    b = -div( grad(exact) ) * v * dx

  scheme = galerkin([a == b], solver="cg", parameters={"newton.linear.verbose": "false"})

  took = -time()
  scheme.solve(uh)
  took += time()

  if rank == 0:
    print(f"Took {took:.6f}")

  grid.writeVTK("laplace-"+name, pointdata={"uh": uh, "exact": exact}, nonconforming=True)

  l2 = integrate(grid, sqrt(dot(uh-exact, uh-exact)), order=5)
  if l2 > 1e-2:
    print(l2)
    sys.exit(1)


########
# MAIN #
########

for file in [line.filename, vertical.filename]:
  gridView = mmesh((reader.gmsh, file), 2)
  hgrid = gridView.hierarchicalGrid
  igrid = hgrid.interfaceGrid

  # Run bulk
  algorithm(gridView, "bulk")

  # Run interface
  algorithm(igrid, "interface")
