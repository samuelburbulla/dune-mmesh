## @example elasticity.py
#  This is an example of how to use finite elements with MMesh and dune-python

from dune.grid import reader
from dune.mmesh import mmesh

import matplotlib
matplotlib.rc( 'image', cmap='jet' )
from matplotlib import pyplot
from dune.fem.plotting import plotPointData as plot
from dune.fem.space import lagrange as solutionSpace
from dune.fem.scheme import galerkin as solutionScheme

from ufl import *
import dune.ufl

dim = 2
file = "../grids/beam" + str(dim) + "d.msh"
gridView = mmesh((reader.gmsh, file), dim)

space = solutionSpace(gridView, dimRange=dim, order=2, storage="istl")
displacement = space.interpolate([0]*dim, name="displacement")

x = SpatialCoordinate(space)
dbc = dune.ufl.DirichletBC(space, as_vector([0]*dim), x[0] < 1e-10)

mu = 1.0
lamb = 0.1
rho = 1e-3
g = 9.81

epsilon = lambda u: 0.5*(nabla_grad(u) + nabla_grad(u).T)
sigma = lambda u: lamb*nabla_div(u)*Identity(dim) + 2*mu*epsilon(u)

u = TrialFunction(space)
v = TestFunction(space)
equation = inner(sigma(u), epsilon(v))*dx == dot(as_vector([0]*(dim-1) + [-rho*g]), v)*dx

scheme = solutionScheme([equation, dbc], solver='cg',
            parameters = {"newton.linear.preconditioning.method": "ilu"} )
info = scheme.solve(target=displacement)

gridView.writeVTK('elasticity', pointdata=[displacement])

if dim == 2:
    fig = pyplot.figure(figsize=(20,10))
    displacement.plot(gridLines=None, figure=(fig, 121), colorbar="horizontal")
    s = sigma(displacement) - (1./3)*tr(sigma(displacement))*Identity(dim)
    von_Mises = sqrt(3./2*inner(s, s))
    plot(von_Mises, grid=gridView, gridLines=None, figure=(fig, 122), colorbar="horizontal")

    from dune.fem.view import geometryGridView
    position = space.interpolate( x+displacement, name="position" )
    beam = geometryGridView( position )
    beam.plot()
