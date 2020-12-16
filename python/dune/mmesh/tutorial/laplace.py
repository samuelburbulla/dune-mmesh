from dune.grid import reader
from dune.mmesh import mmesh, trace, skeleton, normal

# world dimension
dim = 2

# create grid
import grids.horizontal as gridfile
gridfile.create(h=0.1, dim=dim)

gridView  = mmesh((reader.gmsh, gridfile.filename), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

from ufl import *
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.ufl import DirichletBC
from dune.fem.function import integrate

#########################################
print("Solve a problem on the bulk grid")
#########################################
space = lagrange(gridView, order=3)
x = SpatialCoordinate(space)
u = TrialFunction(space)
v = TestFunction(space)
exact = sin(x[0]*x[1]*4*pi)
a = inner(grad(u), grad(v)) * dx
b = -div(grad(exact)) * v * dx
uh = space.interpolate(0, name="solution")
scheme = galerkin([a==b, DirichletBC(space, exact)], solver="cg")
scheme.solve(target=uh)
print("  error", integrate(gridView, dot(uh-exact, uh-exact), order=5))
gridView.writeVTK("laplace-bulk",
    pointdata={"solution": uh, "exact": exact})

############################################################
print("Solve a problem on the interface grid (no coupling)")
############################################################
ispace = lagrange(igridView, order=3)
ix = SpatialCoordinate(ispace)
iu = TrialFunction(ispace)
iv = TestFunction(ispace)
if dim == 3:
    iexact = sin(0.5*ix[1]*4*pi)
else:
    iexact = sin(ix[0]*0.5*4*pi)
ia = inner(grad(iu), grad(iv)) * dx
ib = -div(grad(iexact)) * iv * dx
iuh = ispace.interpolate(0, name="solution")
ischeme = galerkin([ia==ib, DirichletBC(ispace, iexact)])
ischeme.solve(target=iuh)
print("  error", integrate(igridView, dot(iuh-iexact, iuh-iexact), order=5))
igridView.writeVTK("laplace-interface",
    pointdata={"solution": iuh, "exact": iexact})

###############################
print("Couple bulk to surface")
###############################
nBulk = FacetNormal(space)
nInterface = FacetNormal(ispace)
assert not nBulk == nInterface

# replace 'iexact' with the traces of 'uh' in the bilinear form and bc
ib = inner(grad(avg(trace(uh))), grad(iv)) * dx
ischeme = galerkin([ia==ib, DirichletBC(ispace, avg(trace(uh)))])
iuh.interpolate(0)
ischeme.solve(target=iuh)
print("  error", integrate(igridView, dot(iuh-iexact, iuh-iexact), order=5))
igridView.writeVTK("laplace-bulk2interface",
    pointdata={"iuh": iuh, "grad(iuh)": grad(iuh)})


###############################
print("Couple surface to bulk")
###############################
skeleton = skeleton(iuh)
a = inner(grad(u), grad(v)) * dx
b = avg(skeleton) * avg(v) * dS
uh = space.interpolate(0, name="solution")
scheme = galerkin([a==b, DirichletBC(space,0)], solver="cg")
scheme.solve(target=uh)
gridView.writeVTK("laplace-interface2bulk",
    pointdata={"skeleton": uh})


################################################
print("Compute jump of grad trace times normal")
################################################
from dune.fem.space import dglagrange
space = dglagrange(gridView, order=1)
x = SpatialCoordinate(space)
u = TrialFunction(space)
v = TestFunction(space)
exact = conditional(x[1]<0.5, x[1], -2*x[1])
uh = space.interpolate(exact, name="uh")
gridView.writeVTK("laplace-tracenormalexact",
    pointdata={"uh": uh}, nonconforming=True)

n = normal(igridView)
jumpgrad = jump(grad(trace(uh)), n)

iexact = -3
print("  error", integrate(igridView, dot(jumpgrad-iexact, jumpgrad-iexact), order=5))
igridView.writeVTK("laplace-tracenormaljump",
    pointdata={"jumpgrad": jumpgrad, "normals": n}, nonconforming=True)
