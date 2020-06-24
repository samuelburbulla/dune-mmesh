## @example poroelasticity.py
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from dune.grid import reader
from dune.mmesh import mmesh

from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.function import integrate

from ufl import *
import dune.ufl

dim = 2
file = "../grids/horizontal" + str(dim) + "d.msh"
gridView = mmesh((reader.gmsh, file), dim)


# Spaces
uspace = lagrange(gridView, dimRange=dim, order=2)
uh = uspace.interpolate([0]*dim, name="uh")
uh_old = uspace.interpolate([0]*dim, name="uh_old")

u = TrialFunction(uspace)
v = TestFunction(uspace)

pspace = lagrange(gridView, dimRange=1, order=1, scalar=True)
ph = pspace.interpolate([0], name="ph")
ph_old = pspace.interpolate([0], name="ph_old")

p = TrialFunction(pspace)
q = TestFunction(pspace)

dt = 0.1
tau = dune.ufl.Constant(dt, name="dt")


# Model
mu = 1.0
lamb = 0.1
alpha = 1.0
K = 1.0
M = 1.0
L = 1e-2

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma = lambda u, p: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u) - alpha*p*Identity(dim)

a = inner(sigma(u, ph), epsilon(v)) * dx

b = inner(p - ph_old, q) / M * dx
b += alpha * inner( div(uh - uh_old), q ) * dx
b += L * inner( p - ph, q ) * dx
b += tau * inner(K * grad(p), grad(q)) * dx


# Scheme (fixed-stress splitting)
x = SpatialCoordinate(uspace)
ubc0 = dune.ufl.DirichletBC(uspace, as_vector([0]*dim), x[0] > -1.0)
uscheme = galerkin([a == 0, ubc0])

pbc0 = dune.ufl.DirichletBC(pspace, 0.0, x[0] > 1e-10)
pbc1 = dune.ufl.DirichletBC(pspace, sin(3.14*x[1]), x[0] < 1e-10)
pscheme = galerkin([b == 0, pbc0, pbc1])

intermediateu = uh.copy()
intermediatep = ph.copy()


# Timeloop
for t in range(11):
    print("t =", t*dt)
    uh_old.assign(uh)
    ph_old.assign(ph)

    for i in range(100):
        print(" i =", i)
        intermediateu.assign(uh)
        intermediatep.assign(ph)
        uscheme.solve(target=uh)
        pscheme.solve(target=ph)

        if integrate(gridView, inner(uh - intermediateu, uh - intermediateu), uspace.order) < 1e-12**2 \
          and integrate(gridView, (ph - intermediatep)**2, pspace.order) < 1e-12**2:
            break

    space3 = lagrange(gridView, dimRange=3, order=2)
    uh3 = space3.interpolate([uh[0], uh[1], 0], name="displacement")
    gridView.writeVTK('poroelasticity-'+str(t), pointdata=[uh3, ph])
