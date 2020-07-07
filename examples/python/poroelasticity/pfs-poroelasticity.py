## @example pfs-poroelasticity.py
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from dune.grid import reader
from dune.mmesh import mmesh, skeleton, domainMarker

from dune.fem.space import *
from dune.fem.scheme import galerkin

from ufl import *
import dune.ufl

from dune.fem import parameter
parameter.append({"fem.verboserank": 0})

dim = 2
file = "middlefull" + str(dim) + "d.msh"
gridView = mmesh((reader.gmsh, file), dim)

space = lagrange(gridView, dimRange=dim+1, order=1)
x = SpatialCoordinate(space)
n = FacetNormal(space)
bc = dune.ufl.DirichletBC(space, as_vector([0]*(dim+1)))

dm = domainMarker(gridView)
pfspace = lagrange(gridView, order=1)
pf = pfspace.interpolate( dm, name="pf" )

mu = 1 # 2.58e7
lamb = 0.3
alpha = 1.0
K = 1 # 9.8e-12
Kf = 1/12.

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma = lambda u: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)

trial = TrialFunction(space)
test = TestFunction(space)
u = as_vector([trial[0], trial[1]])
uu = as_vector([test[0], test[1]])
p = trial[2]
pp = test[2]

beta = 1e6
k = 1e-14
a = inner(sigma(u) - alpha*p*Identity(dim), epsilon(uu)) * ((1-pf)**2 + k) * dx
a += inner( (K*(1-pf)+Kf*pf) * grad(p), grad(pp) ) * dx

# DirichletBC
a += beta * inner(as_vector([0]*dim) - u, uu) * ds
a += beta * (0 - p) * pp * ds
a += beta * (0.1 - p) * pp * pf*dx

scheme = galerkin([a == 0],
    solver=('suitesparse', 'umfpack'),
    parameters = {"newton.verbose": True} )

solution = space.interpolate([0]*(dim+1), name="solution")
scheme.solve(target=solution)

gridView.writeVTK('pfs-poroelasticity', pointdata={"displacement": [solution[0], solution[1], 0], "pressure": solution[2], "pf": pf})
