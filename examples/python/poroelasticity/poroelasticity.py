## @example poroelasticity.py
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from dune.grid import reader
from dune.mmesh import mmesh, skeleton, trace, normal

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

E = 1 #15.96e9
nu = 0.2

lamb = E*nu/((1+nu)*(1-2*nu))
mu = E/(2*(1+nu))
alpha = 0.79
K = 1 #2e-11
aperture = 1
K_gamma = 83.3
# f = 2e-3
p0 = 0.1

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

ispace = lagrange(igridView, order=1)
pgamma = ispace.interpolate(0, name="pgamma")
eta = aperture / K_gamma    # actually: normal component of K_gamma

pg = skeleton( pgamma )('+')
b = 0
# normal stress is -p
b += -pg * inner(uu('+'), normal('+')) * I*dS
b += -pg * inner(uu('-'), normal('-')) * I*dS
# p_bulk = p_gamma
b -= -0.5 * eta * ( - K_gamma * inner( grad(p)('+'), normal('+') ) * pp('+') ) * I*dS
b -= (p('+') - pg) * pp('+') * I*dS
b -= -0.5 * eta * ( - K_gamma * inner( grad(p)('-'), normal('+') ) * pp('-') ) * I*dS
b -= (p('-') - pg) * pp('-') * I*dS

# b += beta * (pg - p('+')) * pp('+') * I*dS
# b -= pg * inner(K*grad(pp('+')), n('+')) * I*dS
# b += inner(K*grad(p('+')), n('+')) * pp('+') * I*dS
# b += inner(K*grad(pp('+')), n('+')) * p('+') * I*dS
# b += beta * (pg - p('-')) * pp('-') * I*dS
# b -= pg * inner(K*grad(pp('-')), n('-')) * I*dS
# b += inner(K*grad(p('-')), n('-')) * pp('-') * I*dS
# b += inner(K*grad(pp('-')), n('-')) * p('-') * I*dS

scheme = galerkin([a == b], solver=('suitesparse', 'umfpack'))
solution = space.interpolate([0]*(dim+1), name="solution")

p_gamma = TrialFunction(ispace)
pp_gamma = TestFunction(ispace)

inormal = igridView.normal
bp_p = trace(solution)[2]('+')
bp_m = trace(solution)[2]('-')

a_gamma = inner( aperture * K_gamma * grad(p_gamma), grad(pp_gamma) ) * dx
a_gamma -= -0.5 * eta * ( - K_gamma * inner( grad(bp_p), inormal ) * pp_gamma ) * dx
a_gamma -= (bp_p - p_gamma) * pp_gamma * dx
a_gamma -= -0.5 * eta * ( - K_gamma * inner( grad(bp_m), inormal ) * pp_gamma ) * dx
a_gamma -= (bp_m - p_gamma) * pp_gamma * dx

# a_gamma += aperture * f * pp_gamma * dx

bc = dune.ufl.DirichletBC( ispace, p0 )
scheme_gamma = galerkin([a_gamma == 0, bc], solver=('suitesparse', 'umfpack'))


solution_old = solution.copy()
solution_gamma_old = pgamma.copy()
for i in range(100):
    solution_old.assign(solution)
    solution_gamma_old.assign(pgamma)

    scheme.solve(target=solution)
    scheme_gamma.solve(target=pgamma)

    sol = solution_old.as_numpy[:]
    sol -= solution.as_numpy
    error = np.dot(sol, sol)

    sol_gamma = solution_gamma_old.as_numpy[:]
    sol_gamma -= pgamma.as_numpy
    error_gamma = np.dot(sol_gamma, sol_gamma)

    print("["+str(i)+"]: errors=", [error, error_gamma], flush=True)

    if max(error, error_gamma) < 1e-12:
        break

c = 5
w = conditional(abs(x[0]-50)<5, 2 * p0 / E * sqrt(c**2 - (x[0]-50)**2), 0)

gridView.writeVTK('poroelasticity', pointdata={"displacement": [solution[0], solution[1], 0], "pressure": solution[2], "w": w},
    nonconforming=True, subsampling=space.order-1)
igridView.writeVTK('poroelasticity-interface', pointdata={"pressure": pgamma})
