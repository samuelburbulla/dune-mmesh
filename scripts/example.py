# Generate mesh
import gmsh
gmsh.initialize()
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
name = "grid.msh"
h = 0.01
gmsh.model.add(name)
kernel = gmsh.model.occ
box = kernel.addRectangle(0, 0, 0, 1, 1)
p0 = kernel.addPoint(0.5, 0.5, 0, h)
p1 = kernel.addPoint(0.25, 0.5, 0, h)
p2 = kernel.addPoint(0.75, 0.75, 0, h)
p3 = kernel.addPoint(0.75, 0.25, 0, h)
lf1 = kernel.addLine(p0, p1)
lf2 = kernel.addLine(p0, p2)
lf3 = kernel.addLine(p0, p3)
kernel.synchronize()
gmsh.model.mesh.embed(1, [lf1, lf2, lf3], 2, box)
gmsh.model.mesh.generate(dim=2)
gmsh.write(name)
gmsh.finalize()


# Grid creation
from dune.grid import reader
from dune.mmesh import mmesh
gridView  = mmesh((reader.gmsh, name), 2)
igridView = gridView.hierarchicalGrid.interfaceGrid


# Coupled problem
from dune.mmesh import trace, skeleton, interfaceIndicator, monolithicSolve

from dune.fem.space import dglagrange, lagrange
space = dglagrange(gridView, order=3)
ispace = lagrange(igridView, order=3)

from ufl import *
u = TrialFunction(space)
v = TestFunction(space)
n = FacetNormal(space)
h = MinFacetEdgeLength(space)
uh = space.interpolate(0, name="uh")

iu = TrialFunction(ispace)
iv = TestFunction(ispace)
iuh = ispace.interpolate(0, name="iuh")

from dune.mmesh import interfaceIndicator
I = interfaceIndicator(igridView)

from dune.ufl import Constant
q = Constant(1, name="q")
beta = Constant(1e2, name="beta")
omega = Constant(1e-6, name="omega")

a  = inner(grad(u), grad(v)) * dx
a += beta / h * inner(jump(u), jump(v)) * (1-I)*dS
a -= dot(dot(avg(grad(u)), n('+')), jump(v)) * (1-I)*dS

a += beta / h * inner(u - 0, v) * ds
a -= dot(dot(grad(u), n), v) * ds

ia  = inner(grad(iu), grad(iv)) * dx
ib  = q * iv * dx

from dune.mmesh import skeleton, trace
a -= (skeleton(iuh)('+') - u('+')) / omega * v('+') * I*dS
a -= (skeleton(iuh)('-') - u('-')) / omega * v('-') * I*dS

ia += (iu - trace(uh)('+')) / omega * iv * dx
ia += (iu - trace(uh)('-')) / omega * iv * dx

from dune.fem.scheme import galerkin
scheme = galerkin([a == 0])
ischeme = galerkin([ia == ib])

from dune.mmesh import monolithicSolve
monolithicSolve(schemes=(scheme, ischeme), targets=(uh, iuh), verbose=True)

# Plot
gridView.writeVTK("example", pointdata=[uh], nonconforming=True)
igridView.writeVTK("interface", pointdata=[iuh])

import matplotlib.pyplot as plt
from dune.fem.plotting import plotPointData as plot
fig = plt.figure()
plot(uh, linewidth=0, clim=[0, 0.17], figure=fig)
plot(iuh, linewidth=0.01, colorbar=None, clim=[0, 0.17], figure=fig)
plt.savefig("plot.png")
