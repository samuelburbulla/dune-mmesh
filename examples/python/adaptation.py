from dune.grid import reader, cartesianDomain
from dune.mmesh import mmesh
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
file = "../grids/tip2d.msh"
gridView = adaptive( grid=mmesh((reader.gmsh, file), dim) )
hgrid = gridView.hierarchicalGrid
igridView = hgrid.interfaceGrid

t = 0
dt = 0.1

# Movement
time = Constant(t, name="time")
def speed():
    return ufl.as_vector([1.0]+[0.0]*(dim-1))

def movement(x):
  return ufl.conditional((x[0] - (0.25 + 0.5 * time))**2 + (x[1] - 0.5)**2 < 1e-3, \
    ufl.as_vector([0.5, 0.0]), \
    ufl.as_vector([0.0, 0.0]) )

def getShifts():
  mapper = igridView.mapper({dune.geometry.vertex: 1})
  shifts = np.zeros((mapper.size, dim))
  for v in igridView.vertices:
    shifts[ mapper.index(v) ] = ufl.as_vector(movement( v.geometry.center )) * dt
  return shifts

space = dglagrange(gridView, dimRange=1, order=1)
trial = ufl.TrialFunction(space)
test = ufl.TestFunction(space)
u = trial[0]
v = test[0]

x = ufl.SpatialCoordinate(space.cell())
n = ufl.FacetNormal(space.cell())
left = ufl.conditional( x[0] < 1e-6, 1.0, 0.0 )
right = ufl.conditional( x[0] > 1 - 1e-6, 1.0, 0.0 )

def f(u):
    return speed() * u

# numerical flux
def g(u, n):
    return ufl.inner( ufl.conditional( ufl.inner(speed(), n('+')) > 0, f( u('+') ), f( u('-') ) ), n('+') )

# geometrical flux
def h(u, n):
    sgn = ufl.inner(movement(x), n('+')) # TODO: use effective edge movement here
    return ufl.conditional( sgn > 0, -sgn * u('+'), -sgn * u('-') )

def u0(x):
    return 1.0 - x[0]

uh_old = space.interpolate(0, name="uh_old")
tau = Constant(dt, name="timeStep")

a = (u - uh_old[0]) / tau * v * ufl.dx
a -= ufl.inner( f(u), ufl.grad(v) ) * ufl.dx
a += g(u, n) * ufl.jump(v) * ufl.dS
# a += h(u, n) * ufl.jump(v) * ufl.dS # TODO add this if h(u, n) is correct
a += ufl.inner(f(1.0), n) * v * left * ufl.ds
a += ufl.inner(f(u), n) * v * right * ufl.ds

scheme = galerkin([a == 0], solver=("suitesparse","umfpack"))
uh = space.interpolate(u0, name="uh")

gridView.writeVTK("transport-mesh")
gridView.writeVTK("transport-"+str(0), pointdata={"u": uh}, nonconforming=True, subsampling=max(0,space.order-1))

for step in range(1, 11):
  print("step =", step)
  hgrid.ensureInterfaceMovement( getShifts() )
  hgrid.markElements()
  adapt([uh])

  gridView.writeVTK("adaptation-"+str(2*step-1), pointdata={"u": uh}, nonconforming=True, subsampling=max(0,space.order-1))

  hgrid.moveInterface( getShifts() )

  t += dt
  time = Constant(t, name="time")

  uh_old.assign(uh)
  scheme.solve(target=uh)
  gridView.writeVTK("adaptation-"+str(2*step), pointdata={"u": uh}, nonconforming=True, subsampling=max(0,space.order-1))
