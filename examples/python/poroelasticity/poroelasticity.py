## @example poroelasticity.py
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from dune.grid import reader
from dune.mmesh import mmesh, skeleton

from dune.fem.space import *
from dune.fem.scheme import galerkin

from ufl import *
import dune.ufl

from dune.fem import parameter
parameter.append({"fem.verboserank": 0})

dim = 2
file = "middle" + str(dim) + "d.msh"
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

space = dglagrange(gridView, dimRange=dim+1, order=1)
x = SpatialCoordinate(space)
n = FacetNormal(space)

lamb = 1.5
mu = 1
alpha = 1.0
K = 1e-3

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma_eff = lambda u: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)
sigma = lambda u, p: sigma_eff(u) - alpha*p*Identity(dim)

trial = TrialFunction(space)
test = TestFunction(space)
u = as_vector([trial[0], trial[1]])
uu = as_vector([test[0], test[1]])
p = trial[2]
pp = test[2]

space_gamma = dglagrange(igridView, order=0)
one = space_gamma.interpolate(1, name="one")
I = avg(skeleton(one))
beta = 1e3
normal = FacetNormal(space.cell())

a = inner(sigma(u, p), epsilon(uu)) * dx
a += inner( K * grad(p), grad(pp) ) * dx
# enforce solution to be continuous except for the interface
a += beta * inner(jump(u), jump(uu)) * (1-I)*dS
a -= dot(dot(avg(sigma(u,p)), n('+')), jump(uu)) * (1-I)*dS
a -= dot(dot(avg(epsilon(uu)), n('+')), jump(u)) * (1-I)*dS
a += beta * inner(jump(p), jump(pp)) * (1-I)*dS
a -= inner(avg(K*grad(p)), n('+')) * jump(pp) * (1-I)*dS
a -= inner(avg(K*grad(pp)), n('+')) * jump(p) * (1-I)*dS

# DirichletBC
a += beta * inner(as_vector([0]*dim) - u, uu) * ds
a -= dot(dot(sigma(u,p), n), uu) * ds
a -= dot(dot(epsilon(uu), n), u) * ds
a += beta * (0 - p) * pp * ds
a -= inner(K*grad(p), n) * pp * ds
a -= inner(K*grad(pp), n) * p * ds

p0 = 0.1
ispace = lagrange(igridView, dimRange=1, order=1)
pgamma = ispace.interpolate(p0, name="pgamma")
pg = skeleton( pgamma )('+')
b = 0
# normal stress is -p
b += -pg * inner(uu('+'), normal('+')) * I*dS
b += -pg * inner(uu('-'), normal('-')) * I*dS
# p_bulk = p_gamma
b += beta * (pg - p('+')) * pp('+') * I*dS
b += beta * (pg - p('-')) * pp('-') * I*dS

scheme = galerkin([a == b],
    solver=('suitesparse', 'umfpack'),
    parameters = {"newton.verbose": True} )

solution = space.interpolate([0]*(dim+1), name="solution")
scheme.solve(target=solution)

c = 5
E = mu * (3*lamb + 2*mu) / (lamb + mu)
w = conditional(abs(x[0]-50)<5, 2 * p0 / E * sqrt(c**2 - (x[0]-50)**2), 0)

gridView.writeVTK('poroelasticity', pointdata={"displacement": [solution[0], solution[1], 0], "pressure": solution[2], "w": w},
    nonconforming=True, subsampling=space.order-1)
igridView.writeVTK('poroelasticity-interface', pointdata={"pressure": pgamma})
