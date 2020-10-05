## @example dfc-poroelasticity.py
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from dune.grid import reader
from dune.mmesh import mmesh, skeleton, trace, normal, interfaceIndicator

from dune.fem.space import *
from dune.fem.scheme import galerkin

from ufl import *
import dune.ufl
import numpy as np

dim = 2
file = "middle" + str(dim) + "d.msh"
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

space = dglagrange(gridView, dimRange=dim+1, order=1)
x = SpatialCoordinate(space)
n = FacetNormal(space)

E = 67e6
nu = 0.3
lamb = E*nu/((1+nu)*(1-2*nu))
mu = E/(2*(1+nu))

eta = 1e-3
alpha = 1
K = 9.8e-12 / eta
aperture = 1e-3
q = 1e-3
d = dune.ufl.Constant(aperture, name="d")
K_gamma = d**2 / (12 * eta)

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma_eff = lambda u: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)
sigma = lambda u, p: sigma_eff(u) - alpha*p*Identity(dim)

trial = TrialFunction(space)
test = TestFunction(space)
u = as_vector([trial[0], trial[1]])
uu = as_vector([test[0], test[1]])
p = trial[2]
pp = test[2]

I = interfaceIndicator(igridView)
beta = dune.ufl.Constant(1e-6, name="beta")

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

ispace = lagrange(igridView, order=1)
pgamma = ispace.interpolate(0, name="pgamma")
zeta = aperture / K

pg = skeleton( pgamma )('+')
b = 0

# normal stress is -p
b += -pg * inner(uu('+'), n('+')) * I*dS
b += -pg * inner(uu('-'), n('-')) * I*dS

# p_bulk = p_gamma
b += beta * (pg - p('+')) * pp('+') * I*dS
# b -= pg * inner(K*grad(pp('+')), n('+')) * I*dS
b += inner(K*grad(p('+')), n('+')) * pp('+') * I*dS
b += inner(K*grad(pp('+')), n('+')) * p('+') * I*dS
b += beta * (pg - p('-')) * pp('-') * I*dS
# b -= pg * inner(K*grad(pp('-')), n('-')) * I*dS
b += inner(K*grad(p('-')), n('-')) * pp('-') * I*dS
b += inner(K*grad(pp('-')), n('-')) * p('-') * I*dS

scheme = galerkin([a == b], solver=('suitesparse', 'umfpack'))
solution = space.interpolate([0]*(dim+1), name="solution")

p_gamma = TrialFunction(ispace)
pp_gamma = TestFunction(ispace)

inormal = normal(igridView)
bp_p = trace(solution)[2]('+')
bp_m = trace(solution)[2]('-')

n_gamma = FacetNormal(ispace)

a_gamma = inner( d * K_gamma * grad(p_gamma), grad(pp_gamma) ) * dx

# coupling
a_gamma -= -0.5 * zeta * ( - K_gamma * inner( grad(bp_p), inormal ) * pp_gamma ) * dx
a_gamma -= beta * (bp_p - p_gamma) * pp_gamma * dx
a_gamma -= -0.5 * zeta * ( - K_gamma * inner( grad(bp_m), inormal ) * pp_gamma ) * dx
a_gamma -= (bp_m - p_gamma) * pp_gamma * dx

# source
f = dune.ufl.Constant(q, name="f")
b_gamma = f * pp_gamma * d * dx

scheme_gamma = galerkin([a_gamma == b_gamma], solver=('suitesparse', 'umfpack'))

from dune.mmesh import monolithicNewton
monolithicNewton(schemes=(scheme, scheme_gamma), targets=(solution, pgamma), verbose=True)

gridView.writeVTK('dfc-poroelasticity', pointdata={"displacement": [solution[0], solution[1], 0], "pressure": solution[2]},
    nonconforming=True, subsampling=space.order-1)
igridView.writeVTK('dfc-poroelasticity-interface', pointdata={"pressure": pgamma})
