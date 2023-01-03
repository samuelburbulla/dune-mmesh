"""Test a specific coupled problem."""

from dune.grid import reader
from dune.mmesh import mmesh, monolithicSolve, iterativeSolve, trace
from dune.mmesh.test.grids import line, plane
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.function import integrate
from dune.ufl import DirichletBC
from ufl import TrialFunction, TestFunction, SpatialCoordinate, dot, grad, dx, sin, pi

# Let   Ω := (0, 1)^d,
#       γ := {x in Ω | x[d-1] = 0.5},
#   ∂ΩD := {x in ∂Ω | x[d-1] = 0 or x[d-1] = 1},
#   ∂ΩN := ∂Ω - ∂ΩD,
#   ∂γD := {x in ∂Ω | x[0] = 0 or x[0] = 1},
#   ∂γN := ∂γ - ∂γD,
#
# We use a finite element method to solve the following problem:
#        -Δu = 0  in Ω,
#  -Δl + u*l = g  in γ,
#          u = 0  on  ∂ΩD,
#  grad(u)*n = 0  on ∂ΩN,
#          l = 0  on  ∂γD,
#  grad(l)*n = 0  on ∂γN.

def coupledproblem(grid):
  dim = grid.dimension
  igrid = grid.hierarchicalGrid.interfaceGrid

  def uExact(x):
    return x[dim - 1]

  def lExact(x):
    return sin(pi * x[0])

  def ig(x):
    return (pi * pi + 0.5) * sin(pi * x[0])

  space = lagrange(grid, order=1)
  ispace = lagrange(igrid, order=1)

  x = SpatialCoordinate(space)
  u = TrialFunction(space)
  v = TestFunction(space)

  ix = SpatialCoordinate(ispace)
  l = TrialFunction(ispace)
  m = TestFunction(ispace)

  uh = space.interpolate(0., name="uh")
  lh = ispace.interpolate(0., name="lambdah")

  b = dot(grad(u), grad(v)) * dx
  dbc_top = DirichletBC(space, 0., x[dim - 1] <= 1e-8)
  dbc_btm = DirichletBC(space, 1., x[dim - 1] >= 1. - 1e-8)

  ib = dot(grad(l), grad(m)) * dx
  ib += trace(uh, igrid)("+") * l * m * dx
  il = ig(ix) * m * dx
  idbc_lft = DirichletBC(ispace, 0., ix[0] <= 1e-8)
  idbc_rgt = DirichletBC(ispace, 0., ix[0] >= 1. - 1e-8)

  umfp = ("suitesparse", "umfpack")
  scheme = galerkin([b == 0., dbc_top, dbc_btm], solver=umfp)
  ischeme = galerkin([ib == il, idbc_lft, idbc_rgt], solver=umfp)

  for coupledSolve in [iterativeSolve, monolithicSolve]:
    coupledSolve(schemes=(scheme, ischeme), targets=(uh, lh), verbose=1)

    err = integrate(grid, abs(uh - uExact(x)), order=5)
    ierr = integrate(igrid, abs(lh - lExact(ix)), order=5)

    assert err < 1e-10
    assert ierr < 5e-3

  grid.writeVTK(f"coupledsolve-{dim}d", pointdata={"uh": uh, "uExact": uExact(x)})
  igrid.writeVTK(f"coupledsolve-{dim}d-interface", pointdata={"lh": lh, "lExact": lExact(ix)})


# 2D
gridView = mmesh((reader.gmsh, line.filename), 2)
coupledproblem(gridView)

# 3D
gridView = mmesh((reader.gmsh, plane.filename), 3)
coupledproblem(gridView)
