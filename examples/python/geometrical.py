from dune.grid import reader
from dune.mmesh import mmesh, domainMarker, skeleton, trace, iterativeSolve, normal
from dune.fem.space import finiteVolume, lagrange, dglagrange
from dune.fem.scheme import galerkin
from dune.ufl import DirichletBC
from dune.fem.view import filteredGridView, geometryGridView
from ufl import *

dim = 2
file = "vertical.tmp.msh"
from onep import vertical
vertical.create(file, h=0.1, hf=0.01)

gridView  = mmesh((reader.gmsh, file), dim)
gamma = gridView.hierarchicalGrid.interfaceGrid

fvspace = finiteVolume(gridView)
dm = fvspace.interpolate(domainMarker(gridView), name="dm")

# Move domain
d0 = 0.1
def d(x):
    return 2 * d0 + d0 * cos(8 * pi * x[1])

positionspace = dglagrange(gridView, dimRange=dim, order=1)
x = SpatialCoordinate(positionspace)
position = positionspace.interpolate(x + (2*dm-1) * 0.5 * d(x) * (dm - (2*dm-1) * x[0]) / 0.5 * as_vector([1, 0]), name="position")
omega = geometryGridView(position)

# Solve coupled problem
fvspaceG = finiteVolume(gamma)
one = fvspaceG.interpolate(1, name="one")
I = avg(skeleton(one, grid=omega))

space = dglagrange(omega, order=1)
spaceG = dglagrange(gamma, order=1)

n = FacetNormal(space.cell())
nG = FacetNormal(spaceG.cell())
xG = SpatialCoordinate(spaceG)

uh = space.interpolate(0, name="solution")
uhG = spaceG.interpolate(0, name="solution")

mu = 1e3

u = TrialFunction(space)
v = TestFunction(space)
a = inner(grad(u), grad(v)) * dx

a += mu * jump(u) * jump(v) * (1-I)*dS
a -= inner(avg(grad(u)), n('+')) * jump(v) * (1-I)*dS
a -= inner(avg(grad(v)), n('+')) * jump(u) * (1-I)*dS

a += mu * (u('+') - skeleton(uhG, grid=omega)('+')) * v('+') * I*dS
a += mu * (u('-') - skeleton(uhG, grid=omega)('-')) * v('-') * I*dS

a += mu * (u - 0) * v * ds
a -= u * inner(grad(v), n) * ds
a -= v * inner(grad(u), n) * ds
scheme = galerkin([a == 0], solver=('suitesparse', 'umfpack'))

uG = TrialFunction(spaceG)
vG = TestFunction(spaceG)
aG = inner(grad(uG), grad(vG)) * dx
aG += mu * jump(uG) * jump(vG) * dS
aG -= inner(avg(grad(uG)), nG('+')) * jump(vG) * dS
aG -= inner(avg(grad(vG)), nG('+')) * jump(uG) * dS
aG -= (1+avg(trace(uh))) * vG * dx
aG += mu * (uG - 0) * vG * ds
schemeG = galerkin([aG == 0], solver=('suitesparse', 'umfpack'))

iterativeSolve((scheme, schemeG), (uh, uhG), verbose=True)

omega.writeVTK("geometrical-omega", pointdata=[uh], nonconforming=True)
gamma.writeVTK("geometrical-gamma", pointdata=[uhG], nonconforming=True)
