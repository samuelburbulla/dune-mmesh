## @example pfs-poroelasticity.py
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from dune.grid import reader
from dune.mmesh import mmesh, skeleton, interfaceIndicator

from dune.fem.space import *
from dune.fem.scheme import galerkin

from ufl import *
import dune.ufl

import numpy as np

from dune.fem import parameter
parameter.append({"fem.verboserank": 0})
solverParameters =\
   {"newton.tolerance": 1e-6,
    "newton.linear.tolerance": 1e-12,
    "newton.linear.preconditioning.method": "jacobi",
    "newton.linear.preconditioning.iterations": 1,
    "newton.linear.preconditioning.relaxation": 1.2,
    "newton.linear.maxiterations": 10000,
    "newton.verbose": True,
    "newton.linear.verbose": False}

dim = 2
file = "middle" + str(dim) + "d.msh"
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid
I = interfaceIndicator(igridView)

pfspace = lagrange(gridView, order=1)
pt = TrialFunction(pfspace)
ppt = TestFunction(pfspace)
l = 1
pfa = pt * ppt * dx
pfa += l**2 * inner( grad(pt), grad(ppt) ) * dx
pfa += 200*avg( (pt - 1) * ppt ) * I*dS
pfscheme = galerkin([pfa == 0], solver=('suitesparse', 'umfpack'))
pf = pfspace.interpolate(0, name="pf")
pfscheme.solve(target=pf)

E = 15.96e9
nu = 0.2
lamb = E*nu/((1+nu)*(1-2*nu))
mu = E/(2*(1+nu))
M = 12.5e9

alpha = 0.79
d = 1
eta = 1e-3
K = 2e-14 / eta
Kf = d**2/12 / eta
exp = 50

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma = lambda u: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)

space = lagrange(gridView, dimRange=dim+1, order=1)
normal = FacetNormal(space.cell())

trial = TrialFunction(space)
test = TestFunction(space)
u = as_vector([trial[0], trial[1]])
uu = as_vector([test[0], test[1]])
p = trial[2]
pp = test[2]

uh = space.interpolate([0]*(dim+1), name="u")
uh_old = uh.copy()

n = grad(pf) / sqrt(dot(grad(pf), grad(pf)))
Kd = (K * Identity(dim)) + pf**exp * Kf*(Identity(dim) - outer(n,n))

dt = 0.5
tau = dune.ufl.Constant(dt, name="dt")
k = 1e-5

storage = lambda s: s[2] / M + alpha * div(as_vector([s[0], s[1]]))

a = inner(sigma(u), epsilon(uu)) * ((1-pf)**2 + k) * dx - inner(alpha*p*Identity(dim), epsilon(uu)) * dx
a += (storage(trial) - storage(uh_old)) * pp * dx
a += tau*dot( dot(Kd, grad(p)), grad(pp) ) * dx

# boundary conditions
bc = dune.ufl.DirichletBC(space, [0,0,0])

# interior NeumannBC
f = dune.ufl.Constant(0.002, name="f")
b = 2*avg( f * pp ) * I*dS

scheme = galerkin([a == b, bc], solver=('suitesparse', 'umfpack'), parameters=solverParameters)

for step in range(100):
    print("step", step)

    uh_old.assign(uh)
    scheme.solve(target=uh)

    gridView.writeVTK(
        'pfs-poroelasticity-'+str(step),
        pointdata={
            "displacement": [uh[0], uh[1], 0],
            "pressure": uh[2],
            "phasefield": pf,
            "n": [n[0], n[1], 0]
        }
    )
