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
Kf = d**2 / (12 * eta)

# pfs params
exp = 50
k = 1e-5

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma_eff = lambda u: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)

space = lagrange(gridView, dimRange=dim+1, order=1)
normal = FacetNormal(space.cell())

trial = TrialFunction(space)
test = TestFunction(space)
u = as_vector([trial[0], trial[1]])
uu = as_vector([test[0], test[1]])
p = trial[2]
pp = test[2]

uh = space.interpolate([0]*(dim+1), name="u")

n = grad(pf) / sqrt(dot(grad(pf), grad(pf)))
Kd = (K * Identity(dim)) + pf**exp * Kf*(Identity(dim) - outer(n,n))

a = inner(sigma_eff(u), epsilon(uu)) * ((1-pf)**2 + k) * dx
a -= inner(alpha*p*Identity(dim), epsilon(uu)) * dx
a += dot( dot(Kd, grad(p)), grad(pp) ) * dx

# boundary conditions
bc = dune.ufl.DirichletBC(space, [0,0,0])

# source
f = dune.ufl.Constant(q, name="f")
b = f * pp * d * pf * dx

scheme = galerkin([a == b, bc], solver=('suitesparse', 'umfpack'))

scheme.solve(target=uh)

gridView.writeVTK(
    'pfs-poroelasticity',
    pointdata={
        "displacement": [uh[0], uh[1], 0],
        "pressure": uh[2],
        "phasefield": pf,
        "n": [n[0], n[1], 0]
    }
)
