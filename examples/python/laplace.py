## @example laplace.py
#  This is an example of how to use MMesh with dune-python

import io
from dune.grid import reader
from dune.mmesh import mmesh, trace, skeleton
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

dim = 2

if dim == 3:
    file = "../grids/horizontal3d.msh"
else:
    file = "../grids/horizontal2d.msh"

# MMesh
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid


import ufl
import dune.ufl
from dune.fem import parameter
from dune.fem.function import integrate, uflFunction
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.ufl import DirichletBC
parameter.append({"fem.verboserank": 0})
solverParameters =\
   {"newton.tolerance": 1e-8,
    "newton.linear.tolerance": 1e-12,
    "newton.linear.preconditioning.method": "jacobi",
    "newton.linear.maxiterations":10000,
    "preconditioning.method": "jacobi",
    "newton.verbose": True,
    "newton.linear.verbose": False}

#####################################################
print("Solve a problem on the bulk grid")
#####################################################
space = lagrange(gridView, order=3)
x = ufl.SpatialCoordinate(space)
u = ufl.TrialFunction(space)
v = ufl.TestFunction(space)
exact = ufl.sin(x[0]*x[1]*4*ufl.pi)
a = ufl.inner(ufl.grad(u),ufl.grad(v)) * ufl.dx
b = -ufl.div(ufl.grad(exact))*v * ufl.dx
uh = space.interpolate(0,name="solution")
scheme = galerkin([a==b,DirichletBC(space,exact)], solver="cg",
                  parameters=solverParameters)
scheme.solve(target=uh)
print("  error", integrate(gridView,ufl.dot(uh-exact,uh-exact),order=5))
gridView.writeVTK("laplace-bulk",
    pointdata={"solution":uh,"exact":exact})

###########################################################
print("Solve a problem on the interface grid (no coupling)")
###########################################################
ispace = lagrange(igridView, order=3)
ix = ufl.SpatialCoordinate(ispace)
iu = ufl.TrialFunction(ispace)
iv = ufl.TestFunction(ispace)
if dim==3:
    iexact = ufl.sin(0.5*ix[1]*4*ufl.pi)
else:
    iexact = ufl.sin(ix[0]*0.5*4*ufl.pi)
ia = ufl.inner(ufl.grad(iu),ufl.grad(iv)) * ufl.dx
ib = -ufl.div(ufl.grad(iexact))*iv * ufl.dx
iuh = ispace.interpolate(0,name="solution")
ischeme = galerkin([ia==ib,DirichletBC(ispace,iexact)])
ischeme.solve(target=iuh)
print("  error", integrate(igridView, ufl.dot(iuh-iexact,iuh-iexact),order=5))
igridView.writeVTK("laplace-interface", pointdata=[iuh])

###############################
print("couple bulk to surface")
###############################
nBulk = ufl.FacetNormal(space)
nInterface = ufl.FacetNormal(ispace)
assert not nBulk == nInterface

# replace 'iexact' with the traces of 'uh' in the bilinear form and bc
ib = ufl.inner(ufl.grad(ufl.avg(trace(uh))), ufl.grad(iv)) * ufl.dx
ischeme = galerkin([ia==ib,DirichletBC(ispace,ufl.avg(trace(uh)))])
iuh.interpolate(0)
ischeme.solve(target=iuh)
print("  error", integrate(igridView, ufl.dot(iuh-iexact,iuh-iexact),order=5))
igridView.writeVTK("laplace-bulk2interface",
                   pointdata=[iuh, trace(uh) ]) # TODO: ufl.grad(iuh) fails for some reason...


###############################
print("couple surface to bulk")
###############################
skeleton = skeleton(iuh)
a = ufl.inner(ufl.grad(u),ufl.grad(v)) * ufl.dx
b = ufl.avg(skeleton)*ufl.avg(v) * ufl.dS
uh = space.interpolate(0,name="solution")
scheme = galerkin([a==b,DirichletBC(space,0)], solver="cg",
                  parameters=solverParameters)
scheme.solve(target=uh)
gridView.writeVTK("laplace-interface2bulk",
    pointdata={"skeleton":uh})
