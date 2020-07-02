import io
from dune.grid import reader
from dune.mmesh import mmesh, trace, domainMarker, moveInterface, interfaceCurvature
from dune.ufl import Space, Constant
import dune.create as create
import numpy as np
import timeit

dim = 2
file = "../grids/circle.msh"

gridView = mmesh((reader.gmsh, file), dim)
hgrid = gridView.hierarchicalGrid

import ufl
from dune.fem.function import integrate
from dune.fem.space import *
from dune.fem.scheme import galerkin
from dune.ufl import DirichletBC, latex

dt = 0.25
T = 100.0

start = timeit.default_timer()

# Start the time loop
for t in range(0, (int)(T / dt)):
    gridView = hgrid.leafView
    igridView = hgrid.interfaceGrid

    fvspace = finiteVolume(gridView)
    marker = fvspace.interpolate(domainMarker(gridView), name="marker")

    mus = [1, 10, 10]
    mu = ufl.conditional( ufl.eq(marker, 3), mus[2], ufl.conditional( ufl.eq(marker, 2), mus[1], mus[0] ) )

    rhos = [1, 5, 1]
    rho = ufl.conditional( ufl.eq(marker, 3), rhos[2], ufl.conditional( ufl.eq(marker, 2), rhos[1], rhos[0] ) )

    grav = ufl.as_vector([0.0, -9.81])

    sigma = 1e-2

    space = lagrange(gridView, dimRange=3, order=2)
    trial = ufl.TrialFunction(space)
    test  = ufl.TestFunction(space)
    u = ufl.as_vector([trial[0], trial[1]])
    uu = ufl.as_vector([test[0], test[1]])
    p = trial[2]
    pp = test[2]

    a  = ( mu * ufl.inner(ufl.grad(u), ufl.grad(uu)) ) * ufl.dx
    a += ( ufl.inner(ufl.grad(p), uu) ) * ufl.dx
    a += ( ufl.div(u)*pp ) * ufl.dx
    b = rho*ufl.inner(grav, uu) * ufl.dx

    # TODO add surface tension
    # b += sigma * ufl.inner(curvature, uu) * ufl.dx

    x = ufl.SpatialCoordinate(space)
    top = DirichletBC(space, [ 0.0, 0.0, 0.0 ], x[1] > 3.0-1e-6)
    noslip = DirichletBC(space, [ 0.0, 0.0, None ], x[1] < 3.0-1e-6)

    model = create.model("integrands", gridView, a == b, top, noslip)

    pspace = lagrange(gridView, dimRange=1, order=1)

    space3 = lagrange(igridView, dimRange=3, order=1)
    space2 = lagrange(igridView, dimRange=2, order=1)

    print('setup ', timeit.default_timer() - start)

    start = timeit.default_timer()
    scheme = galerkin(model, space=space, solver=("suitesparse","umfpack"))
    vh = space.interpolate([0,0,0], name="solution")
    scheme.solve(target=vh)
    print('solve ', timeit.default_timer() - start)

    start = timeit.default_timer()
    vel = space.interpolate([vh[0], vh[1], 0], name="velocity")
    p = pspace.interpolate([vh[2]], name="pressure")
    gridView.writeVTK("twophase"+str(t), pointdata=[vel, p], celldata=[marker])

    trace = trace(igridView, vel)
    movement = space3.interpolate(trace, name="movement")
    curvature = interfaceCurvature(igridView)
    igridView.writeVTK("interface"+str(t), pointdata=[movement, curvature])
    ivel = space2.interpolate([dt*movement[0], dt*movement[1]], name="movement")
    print('write ', timeit.default_timer() - start)

    start = timeit.default_timer()
    moveInterface(hgrid, ivel)
    print('move ', timeit.default_timer() - start)

    print("t =", t*dt)
