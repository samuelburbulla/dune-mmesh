import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

import io
from dune.grid import reader
from dune.mmesh import mmesh, trace, skeleton, domainMarker
from dune.ufl import Space, Constant
import dune.create as create
import numpy as np
import timeit

dim = 2
file = "../grids/circle.msh"

gridView = mmesh((reader.gmsh, file), dim)
hgrid = gridView.hierarchicalGrid

import ufl
from dune.fem import parameter, adapt
from dune.fem.function import integrate
from dune.fem.space import *
from dune.fem.scheme import galerkin
from dune.fem.view import adaptiveLeafGridView as adaptive
from dune.ufl import DirichletBC, latex

parameter.append( { "fem.verboserank": 0,
                    "fem.adaptation.method": "callback" } )

dt = 0.25
T = 100.0

start = timeit.default_timer()

gridView = adaptive(hgrid.leafView)
igridView = hgrid.interfaceGrid

fvspace = finiteVolume(gridView)
marker = fvspace.interpolate(domainMarker(gridView), name="marker")

mus = [1, 10, 10]
mu = ufl.conditional( ufl.eq(marker, 3), mus[2], ufl.conditional( ufl.eq(marker, 2), mus[1], mus[0] ) )

rhos = [1, 5, 1]
rho = ufl.conditional( ufl.eq(marker, 3), rhos[2], ufl.conditional( ufl.eq(marker, 2), rhos[1], rhos[0] ) )

grav = ufl.as_vector([0.0, -9.81])

sigma = 1e-2

space = dglagrange(gridView, dimRange=3, order=2)
trial = ufl.TrialFunction(space)
test  = ufl.TestFunction(space)
u = ufl.as_vector([trial[0], trial[1]])
uu = ufl.as_vector([test[0], test[1]])
p = trial[2]
pp = test[2]
beta = 1e6

x = ufl.SpatialCoordinate(space)
bcu = ufl.as_vector([ 0.0, 0.0 ])
bcp = ufl.conditional(x[1] > 3.0-1e-6, 1.0, 0.0)

a  = ( mu * ufl.inner(ufl.grad(u), ufl.grad(uu)) ) * ufl.dx
a += ( ufl.inner(ufl.grad(p), uu) ) * ufl.dx
a += ( ufl.div(u)*pp ) * ufl.dx
# continuity
a += beta * ufl.inner( ufl.jump(u), ufl.jump(uu) ) * ufl.dS
a += beta * ( ufl.jump(p) * ufl.jump(pp) ) * ufl.dS
# boundary conditions
a += beta * ufl.inner( (u - bcu), uu ) * ufl.ds
a += beta * ( (p - 0) * pp ) *  bcp * ufl.ds
# gravity source term
b = rho*ufl.inner(grav, uu) * ufl.dx

# surface tension
# kspace = lagrange(igridView, dimRange=dim, order=1)
# k = ufl.TrialFunction(kspace)
# kk = ufl.TestFunction(kspace)
# x = ufl.SpatialCoordinate(kspace)
# ka = ufl.inner(k, kk) * ufl.dx - ufl.inner(ufl.grad(x), ufl.grad(kk)) * ufl.dx
# curvature = kspace.interpolate([0]*dim, name="curvature")
# kscheme = galerkin([ka == 0], kspace)
#
# b += sigma * ufl.inner(skeleton(curvature)('+'), uu) * ufl.dS

scheme = galerkin([a == b], solver=("suitesparse","umfpack"),
    parameters = {"newton.verbose": True, "fem.solver.preconditioning.method": "ilu"})
vh = space.interpolate([0,0,0], name="solution")
print('setup ', timeit.default_timer() - start)

# Start the time loop
for t in range(0, (int)(T / dt)):

    start = timeit.default_timer()
    # kscheme.solve(target=curvature)
    scheme.solve(target=vh)
    print('solve ', timeit.default_timer() - start)

    start = timeit.default_timer()
    gridView.writeVTK("twophase"+str(t),
        pointdata={"velocity": [vh[0], vh[1], 0], "pressure": vh[2]},
        celldata=[marker])

    ivel = trace(vh)
    igridView.writeVTK("interface"+str(t),
        pointdata={"velocity": [ivel[0], ivel[1], 0]})#, "curvatureTimesNormal": curvature})
    print('write ', timeit.default_timer() - start)

    start = timeit.default_timer()

    # move interface
    def getShifts():
      mapper = igridView.mapper({dune.geometry.vertex: 1})
      shifts = np.zeros((mapper.size, dim))
      for e in igridView.elements:
          for v in e.subEntities(igridView.dimension):
              m = ivel(e, v.geometry.center)
              shifts[ mapper.index(v) ][0] += 0.5*dt*m[0]
              shifts[ mapper.index(v) ][1] += 0.5*dt*m[1]
      return shifts

    hgrid.markElements()
    for i in range(3):
      hgrid.ensureInterfaceMovement( getShifts() )
      adapt([vh, marker])

    hgrid.moveInterface( getShifts() )
    print('move ', timeit.default_timer() - start)

    print("t =", t*dt)
