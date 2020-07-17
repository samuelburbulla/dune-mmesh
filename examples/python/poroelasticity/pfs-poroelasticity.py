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
fvspace = finiteVolume(gridView)
dmfv = fvspace.interpolate(dm, name="dm")
pfspace = lagrange(gridView, order=1)
pt = TrialFunction(pfspace)
ppt = TestFunction(pfspace)
l = 0.1
pfa = pt * ppt * dx + l**2 * inner( grad(pt), grad(ppt) ) * dx + 1e2*(pt-1) * ppt * dmfv*dx
pfscheme = galerkin([pfa == 0], solver=('suitesparse', 'umfpack'))
pf = pfspace.interpolate(dm, name="pf")
pfscheme.solve(target=pf)

lamb = 1.5
mu = 1
alpha = 1.0
K = 1e-3
d = 1
Kf = d**2/12

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma = lambda u: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)

trial = TrialFunction(space)
test = TestFunction(space)
u = as_vector([trial[0], trial[1]])
uu = as_vector([test[0], test[1]])
p = trial[2]
pp = test[2]

n = grad(pf) / sqrt(dot(grad(pf), grad(pf)))
Kd = (K * Identity(dim)) + pf * (Kf-K)*(Identity(dim) - outer(n,n))

a = inner(sigma(u), epsilon(uu)) * ((1-pf)**2 + 1e-5) * dx - inner(alpha*p*Identity(dim), epsilon(uu)) * dx
a += dot( dot(Kd, grad(p)), grad(pp) ) * dx

# boundary conditions
bc = dune.ufl.DirichletBC(space, [0,0,0])
# source
p0 = 0.1
b = (p0 - p) * pp * dmfv*dx

scheme = galerkin([a == b, bc],
    solver=('suitesparse', 'umfpack'),
    parameters = {"newton.verbose": True} )

solution = space.interpolate([0]*(dim+1), name="solution")
scheme.solve(target=solution)

c = 5
E = mu * (3*lamb + 2*mu) / (lamb + mu)
sigma = 0 # TODO
w = conditional(abs(x[0]-50)<5, 2 * (1 - sigma**2) * p0 / E * sqrt(c**2 - (x[0]-50)**2), 0)
gridView.writeVTK('pfs-poroelasticity', pointdata={"displacement": [solution[0], solution[1], 0], "pressure": solution[2], "phasefield": pf, "w": w},
    celldata={"domainMarker": dmfv})
