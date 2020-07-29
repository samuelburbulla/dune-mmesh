from dune.grid import reader, cartesianDomain
from dune.mmesh import mmesh, skeleton, edgemovement, cellVolumes
from dune.ufl import Constant
import ufl
import numpy as np
from dune.fem import parameter, adapt
from dune.fem.space import *
from dune.fem.scheme import galerkin
from dune.fem.function import integrate
from dune.fem.view import adaptiveLeafGridView as adaptive
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

parameter.append( { "fem.verboserank": -1,
                    "fem.adaptation.method": "callback" } )

dim = 2
file = "../grids/quad.msh"
gridView = adaptive( grid=mmesh((reader.gmsh, file), dim) )
hgrid = gridView.hierarchicalGrid
igridView = hgrid.interfaceGrid

t = 0
tEnd = 1
dt = 0.01

# Movement
time = Constant(t, name="time")
def speed():
    return ufl.as_vector([0.0, 1.0])

def movement(x):
  return ufl.as_vector([0.0, 1.0])

def getShifts():
  mapper = igridView.mapper({dune.geometry.vertex: 1})
  shifts = np.zeros((mapper.size, dim))
  for v in igridView.vertices:
    shifts[ mapper.index(v) ] = ufl.as_vector(movement( v.geometry.center ))
  return shifts

shifts = getShifts()
em = edgemovement(gridView, shifts)

space = dglagrange(gridView, dimRange=1, order=1)
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
def gBnd(u, n):
    sgn = ufl.inner(speed(), n)
    return ufl.inner( ufl.conditional( sgn > 0, f(u), f(uexact(x)) ), n )

# geometrical flux
def h(u, n):
    sgn = ufl.inner(em('+'), n('+'))
    return ufl.conditional( sgn > 0, sgn * u('+'), sgn * u('-') )

def u0(x):
    return 3.0*ufl.conditional(abs(x[0]-0.5)<0.2, 1.0, 0.0)*ufl.conditional(abs(x[1]-0.5)<0.2, 1.0, 0.0) + ufl.sin(5*3.14*x[0]) + ufl.cos(6*3.14*x[1])

def uexact(x):
    return u0(x - time*speed())

uh_old = space.interpolate(0, name="uh_old")
tau = Constant(dt, name="timeStep")

cellVolume = ufl.CellVolume(space)
vol_old = cellVolumes(gridView)

a = (u - uh_old[0] * vol_old / cellVolume) / tau * v * ufl.dx
a -= ufl.inner( f(u), ufl.grad(v) ) * ufl.dx
a += g(u, n) * ufl.jump(v) * ufl.dS
a += ufl.inner( em * u, ufl.grad(v) ) * ufl.dx
a -= h(u, n) * ufl.jump(v) * ufl.dS
a += gBnd(u, n) * v * ufl.ds

scheme = galerkin([a == 0], solver=("suitesparse","umfpack"))
uh = space.interpolate(u0(x), name="uh")

def writeVTK(step):
    gridView.writeVTK("adaptation-subsampled-"+str(step), pointdata={"u": uh}, nonconforming=True, subsampling=space.order-1)
    gridView.writeVTK("adaptation-"+str(step), pointdata={"u": uh, "em" : em}, nonconforming=True)

writeVTK(0)

for step in range(1, int(tEnd/dt)+1):
  print("step =", step)

  while hgrid.ensureInterfaceMovement(getShifts()*dt):
    print("ensure1")
    adapt([uh])

  print("mark")
  hgrid.markElements()
  adapt([uh])

  while hgrid.ensureInterfaceMovement(getShifts()*dt):
    print("ensure2")
    adapt([uh])

  vol_old.assign( cellVolumes(gridView) )

  hgrid.moveInterface(shifts*dt)

  em = edgemovement(gridView, shifts)

  t += dt
  time.assign(t)

  uh_old.assign(uh)
  scheme.solve(target=uh)

  writeVTK(step)

print("error =", integrate(gridView, ufl.sqrt( ufl.dot(uh[0] - uexact(x), uh[0] - uexact(x)) ), order=5 ) )
