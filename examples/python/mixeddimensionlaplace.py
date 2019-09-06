#!/usr/bin/python3
import ufl
from ufl import grad, div, dot, dx, ds, inner, sin, cos, pi, exp, sqrt
import dune.ufl
import dune.grid
import dune.mmesh
import dune.fem

endTime  = 0.05
saveInterval = 0.05

# YaspGrid
# gridView = dune.grid.structuredGrid([0, 0], [1, 1], [140, 140])

# MMesh
domain = dune.grid.cartesianDomain([0, 0], [1, 1], [100, 100])
gridView = dune.mmesh.mmesh(domain)

interfaceGridView = dune.mmesh.mmeshInterfaceGrid(domain)


space = dune.fem.space.lagrange(gridView, order=1)
u     = ufl.TrialFunction(space)
phi   = ufl.TestFunction(space)
x     = ufl.SpatialCoordinate(space)
dt    = dune.ufl.Constant(5e-3, "timeStep")
t     = dune.ufl.Constant(0.0, "time")

# define storage for discrete solutions
uh     = space.interpolate(0, name="uh")
uh_old = uh.copy()

# initial and exact solution
initial = 1/2*(x[0]**2+x[1]**2) - 1/3*(x[0]**3 - x[1]**3) + 1
exact   = exp(-2*t)*(initial - 1) + 1

# problem definition

# space form
abs_du = lambda v: sqrt(inner(grad(v), grad(v)))
K = lambda v: 2/(1 + sqrt(1 + 4*abs_du(v)))
f = -2*exp(-2*t)*(initial - 1) - div( K(exact)*grad(exact) )
g = K(exact)*grad(exact)
n = ufl.FacetNormal(space)
xForm = (dot(K(u)*grad(u), grad(phi)) - f*phi) * dx - dot(g, n) * phi * ds

# add time discretization
form = dot(u - uh_old, phi) * dx + dt * xForm

scheme = dune.fem.scheme.galerkin(form == 0, solver="cg")

nextSaveTime = saveInterval

uh.interpolate(initial)
vtk = gridView.sequencedVTK("parabolic", pointdata=[uh])
vtk()

while t.value < endTime:
    uh_old.assign(uh)
    info = scheme.solve(target=uh)
    t.value += dt.value

    print("Computed solution at time", t.value,
              "iterations: ", info["linear_iterations"],
              "#Ent: ", gridView.size(0) )
    if t.value >= nextSaveTime or t.value >= endTime:
        vtk()
        nextSaveTime += saveInterval
