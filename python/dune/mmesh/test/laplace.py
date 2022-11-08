from dune.grid import reader
from dune.mmesh import mmesh
from dune.fem.space import dglagrange
from dune.fem.function import integrate
from ufl import *
from dune.ufl import Constant
from dune.fem.scheme import galerkin
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

  A = inner(grad(u), grad(v)) * dx
  A += beta / h * inner(jump(u), jump(v)) * dS
  A -= dot(dot(avg(grad(u)), n('+')), jump(v)) * dS
  A -= dot(dot(avg(grad(v)), n('+')), jump(u)) * dS
  A += beta / hBnd * inner(u - exact, v) * ds
  A -= dot(dot(grad(u), n), v) * ds
  A -= dot(dot(grad(v), n), u - exact) * ds

  if name == "interface":
    b = -exact.dx(0).dx(0) * v * dx
  else:
    b = -div( grad(exact) ) * v * dx

  scheme = galerkin([A == b], solver='cg', parameters={"newton.linear.verbose": "false"})

  took = -time()
  scheme.solve(uh)
  took += time()

  if rank == 0:
    print(f"Took {took:.6f}")

  grid.writeVTK("laplace-"+name, pointdata={"uh": uh, "exact": exact}, nonconforming=True)

  l2 = integrate(grid, sqrt(dot(uh-exact, uh-exact)), order=5)
  if(l2 > 1e-2):
    print(l2)
    raise


########
# MAIN #
########


from dune.mmesh.test.grids import line, vertical

for file in [line.filename, vertical.filename]:
  grid = mmesh((reader.gmsh, file), 2)
  hgrid = grid.hierarchicalGrid
  igrid = hgrid.interfaceGrid

  # Run bulk
  algorithm(grid, "bulk")

  # Run interface
  algorithm(igrid, "interface")
