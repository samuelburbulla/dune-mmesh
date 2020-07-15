from dune.grid import reader, cartesianDomain
from dune.mmesh import mmesh, skeleton, edgemovement, cellVolumes
from dune.ufl import Constant
import ufl
import numpy as np
from dune.fem import parameter, adapt
from dune.fem.space import *
from dune.fem.scheme import galerkin
from dune.fem.view import adaptiveLeafGridView as adaptive
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

parameter.append( { "fem.verboserank": -1,
                    "fem.adaptation.method": "callback" } )

dim = 2
file = "../grids/horizontal2d.msh"
gridView = adaptive( grid=mmesh((reader.gmsh, file), dim) )
hgrid = gridView.hierarchicalGrid
igridView = hgrid.interfaceGrid

t = 0
dt = 0.01

# Movement
time = Constant(t, name="time")
def speed():
    return ufl.as_vector([0.0, -1.0])

def movement(x):
  return ufl.as_vector([0.0, -1.0])

def getShifts():
  mapper = igridView.mapper({dune.geometry.vertex: 1})
  shifts = np.zeros((mapper.size, dim))
  for v in igridView.vertices:
    shifts[ mapper.index(v) ] = ufl.as_vector(movement( v.geometry.center ))
  return shifts

shifts = getShifts()
em = edgemovement(gridView, shifts)

space = dglagrange(gridView, dimRange=1, order=0)
trial = ufl.TrialFunction(space)
test = ufl.TestFunction(space)
u = trial[0]
v = test[0]

x = ufl.SpatialCoordinate(space.cell())
n = ufl.FacetNormal(space.cell())
top = ufl.conditional( x[1] < 1e-6, 1.0, 0.0 )
bottom = ufl.conditional( x[1] > 1 - 1e-6, 1.0, 0.0 )

def f(u):
    return speed() * u

# numerical flux
def g(u, n):
    sgn = ufl.inner(speed(), n('+'))
    return ufl.inner( ufl.conditional( sgn > 0, f( u('+') ), f( u('-') ) ), n('+') )

# geometrical flux
def h(u, n):
    sgn = ufl.inner(em('+'), n('+'))
    return ufl.conditional( sgn > 0, sgn * u('+'), sgn * u('-') )

def u0(x):
    return ufl.conditional( x[1] > 0.5, 1.0, -1.0 )

uh_old = space.interpolate(0, name="uh_old")
tau = Constant(dt, name="timeStep")

cellVolume = ufl.CellVolume(space)

a = (u - uh_old[0] / cellVolume) / tau * v * ufl.dx
a -= ufl.inner( f(u), ufl.grad(v) ) * ufl.dx
a += g(u, n) * ufl.jump(v) * ufl.dS
a -= h(u, n) * ufl.jump(v) * ufl.dS
a += ufl.inner(f(u), n) * v * ufl.ds

scheme = galerkin([a == 0], solver=("suitesparse","umfpack"))
uh = space.interpolate(u0(x), name="uh")

def writeVTK(step):
    gridView.writeVTK("adaptation-"+str(step), pointdata={"u": uh, "em" : em}, nonconforming=True, subsampling=max(0,space.order-1))

writeVTK(0)

for step in range(1, 31):
  print("step =", step)

  hgrid.ensureInterfaceMovement(shifts*dt)
  hgrid.markElements()
  adapt([uh])

  # multiply cell values by (old) cell volumes
  vol_old = cellVolumes(gridView)
  uh_old.assign(uh)
  uh_old.as_numpy[:] *= vol_old.as_numpy[:]

  shifts = getShifts()
  hgrid.moveInterface(shifts*dt)

  em = edgemovement(gridView, shifts)

  scheme.solve(target=uh)

  t += dt
  writeVTK(step)
