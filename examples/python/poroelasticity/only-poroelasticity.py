## @example only-poroelasticity.py
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from dune.grid import reader
from dune.mmesh import mmesh, skeleton

from dune.fem.space import *
from dune.fem.scheme import galerkin

from ufl import *
import dune.ufl

from dune.fem.function import integrate

import pygmsh
import meshio

h = 0.01
dim = 2
file = "box" + str(dim) + "d.tmp.msh"
geom = pygmsh.built_in.Geometry()
geom.add_rectangle(0, 1, 0, 1, 0, lcar=h)
mesh = pygmsh.generate_mesh(geom, verbose=False)
meshio.gmsh.write(filename=file, mesh=mesh, fmt_version="2.2", binary=False)
gridView = mmesh((reader.gmsh, file), dim)

space = lagrange(gridView, dimRange=dim+1, order=1)
x = SpatialCoordinate(space)

lamb = 1
mu = 1
alpha = 1
K = 1
d = 1
p0 = 0.1

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma = lambda u: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)

# manufactured solution
hat = 16*x[0]*(1-x[0])*x[1]*(1-x[1])
sgn = x - as_vector([0.5]*dim)
uexact = as_vector([hat*(1-hat)*sgn[0], hat*(1-hat)*sgn[1]])
pexact = p0 * hat
f = div( sigma(uexact) - alpha*pexact*Identity(dim) )
q = div( -K * grad(pexact) )

trial = TrialFunction(space)
test = TestFunction(space)
u = as_vector([trial[0], trial[1]])
uu = as_vector([test[0], test[1]])
p = trial[2]
pp = test[2]

a = -inner(sigma(u) - alpha*p*Identity(dim), epsilon(uu)) * dx
a += inner(K * grad(p), grad(pp)) * dx
b = inner(f, uu) * dx
b += inner(q, pp) * dx

bc = dune.ufl.DirichletBC(space, [uexact[0], uexact[1], pexact])
scheme = galerkin([a == b, bc], solver=('suitesparse', 'umfpack'))

solution = space.interpolate([0]*(dim+1), name="solution")
scheme.solve(target=solution)

error = sqrt( integrate(gridView, (solution[0]-uexact[0])**2 + (solution[1]-uexact[1])**2 + (solution[2]-pexact)**2, order=5) )
print("error =", error)

gridView.writeVTK('only-poroelasticity',
    pointdata={
        "u":      [solution[0], solution[1], 0],
        "p":      solution[2],
        "uexact": [uexact[0], uexact[1], 0],
        "pexact": pexact,
        "udiff":  [uexact[0]-solution[0], uexact[1]-solution[1], 0],
        "pdiff":  pexact-solution[2]
    }
)
