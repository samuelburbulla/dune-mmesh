from dune.grid import reader
from dune.mmesh import mmesh, monolithicSolve, skeleton, trace
from dune.mmesh.test.grids import line, plane
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.function import integrate
from dune.ufl import DirichletBC, Constant
from ufl import *

# Let   Ω := (0, 1)^d,
#       γ := {x in Ω | x[d-1] = 0.5},
#     ∂ΩD := {x in ∂Ω | x[d-1] = 0 or x[d-1] = 1},
#     ∂ΩN := ∂Ω - ∂ΩD,
#     ∂γD := {x in ∂Ω | x[0] = 0 or x[0] = 1},
#     ∂γN := ∂γ - ∂γD,
#
# We use a finite element method to solve the following problem:
#          -Δu = 0  in Ω,
#    -Δλ + u*λ = g  in γ,
#            u = 0  on ∂ΩD,
#    grad(u)*n = 0  on ∂ΩN,
#            λ = 0  on ∂γD,
#    grad(λ)*n = 0  on ∂γN.

def coupledproblem(grid):
    dim = grid.dimension
    igrid = grid.hierarchicalGrid.interfaceGrid

    u_exact = lambda x : x[dim - 1]
    λ_exact = lambda x : sin(pi * x[0])
    ig = lambda x : (pi * pi + 0.5) * sin(pi * x[0])

    space = lagrange(grid, order=1)
    ispace = lagrange(igrid, order=1)

    x = SpatialCoordinate(space)
    u = TrialFunction(space)
    v = TestFunction(space)

    ix = SpatialCoordinate(ispace)
    λ = TrialFunction(ispace)
    μ = TestFunction(ispace)

    uh = space.interpolate(0., name="uh")
    λh = ispace.interpolate(0., name="lambdah")

    b = dot(grad(u), grad(v)) * dx
    dbc_top = DirichletBC(space, 0., x[dim - 1] <= 1e-8)
    dbc_btm = DirichletBC(space, 1., x[dim - 1] >= 1. - 1e-8)

    ib = dot(grad(λ), grad(μ)) * dx
    ib += trace(uh, igrid)('+') * λ * μ * dx
    il = ig(ix) * μ * dx
    idbc_lft = DirichletBC(ispace, 0., ix[0] <= 1e-8)
    idbc_rgt = DirichletBC(ispace, 0., ix[0] >= 1. - 1e-8)

    umfp = ("suitesparse", "umfpack")
    scheme = galerkin([b == 0., dbc_top, dbc_btm], solver=umfp)
    ischeme = galerkin([ib == il, idbc_lft, idbc_rgt], solver=umfp)
    monolithicSolve(schemes=(scheme, ischeme), targets=(uh, λh), verbose=1)

    grid.writeVTK(f"coupledsolve-{dim}d", pointdata={"uh": uh, "u_exact": u_exact(x)})
    igrid.writeVTK(f"coupledsolve-{dim}d-interface", pointdata={"λh": λh, "λ_exact": λ_exact(ix)})

    err = integrate(grid, abs(uh - u_exact(x)), order=5)
    ierr = integrate(igrid, abs(λh - λ_exact(ix)), order=5)

    assert(err < 1e-10)
    assert(ierr < 5e-3)

# 2D
grid = mmesh((reader.gmsh, line.filename), 2)
coupledproblem(grid)

# 3D
grid  = mmesh((reader.gmsh, plane.filename), 3)
coupledproblem(grid)
