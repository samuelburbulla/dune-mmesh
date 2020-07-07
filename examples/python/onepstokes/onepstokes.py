import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

import horizontal
from dune.grid import reader
from dune.mmesh import mmesh, trace, skeleton, domainMarker, interfaceIndicator, normal

import numpy as np
from scipy.sparse import linalg
from math import log

from ufl import *
import dune.ufl
from dune.fem.function import integrate
from dune.fem.space import dglagrange, lagrange, finiteVolume
from dune.fem.scheme import galerkin
from dune.fem.operator import linear as linearOperator
from dune.fem.view import adaptiveLeafGridView as adaptive

########
# Grid #
########
file = "horizontal.msh"
horizontal.create(file, h=0.5, hf=0.05)

dim = 2
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

##########
# Spaces #
##########
order = 1
space = lagrange(gridView, order=order)
p   = TrialFunction(space)
phi = TestFunction(space)
ph = space.interpolate(0, name="pressure")

ispace = lagrange(igridView, dimRange=3, order=order)
trial  = TrialFunction(ispace)
test   = TestFunction(ispace)
P   = trial[0]
Phi = test[0]
V_n = trial[1]
W_n = test[1]
V_t = trial[2]
W_t = test[2]
ih   = ispace.interpolate([0,0,0], name="isolution")
# Ph   = ih[0]
# Vh_n = ih[1]
# Vh_t = ih[2]

x = SpatialCoordinate(space)
ix = SpatialCoordinate(ispace)
I = interfaceIndicator(igridView)

fvspace = finiteVolume(igridView, dimRange=dim)
inormal = fvspace.interpolate(normal(igridView), name="normal")
tangential = as_vector([-inormal[1], inormal[0]])


######################
# Problem parameters #
######################

mu = 1e-3   # [Pa s]
K = 1e-8    # [m^2]
aBJ = 1
d = 1e-2    # [m]

right = conditional(x[0] > 1-1e-6, 1, 0)
neumann = 1-right
vbar = 0
dirichlet = right
pbar = 0

ileft = conditional(ix[0] < 1e-6, 1, 0)
idirichlet = ileft
vbarf_n = 0
vbarf_t = 0.1 + 1e-6 - d/30
ineumann = 1-ileft
H = as_vector([0,0])

K1 = K
K2 = K
largeBJterm = ( 1 + sqrt(K1) / (d*aBJ) + sqrt(K2) / (d*aBJ) )
largeBJterm /= d/12. + sqrt(K1) / (3*aBJ) + sqrt(K2) / (3*aBJ) + sqrt(K1 * K2) / (d*aBJ**2)


################
# Bulk problem #
################

skel = skeleton(ih)
skelPh = skel[0]('+')
skelVh_n = skel[1]('+')
# (4.4) - (4.7)
a = K / mu * inner(grad(p), grad(phi)) * dx
a += d / (2 * mu) * avg(p * phi) * I*dS
f_pm = skelVh_n * jump(phi) * I*dS
e_pm = d / (2 * mu) * avg(skelPh) * avg(phi) * I*dS
L_pm = 0 if vbar == 0 else vbar * phi * neumann * ds

BC = dune.ufl.DirichletBC(space, pbar, dirichlet)

# consistency terms
# beta0 = 10
# h = MaxFacetEdgeLength(space.cell())
# beta = beta0 * (order + 1) * (order + dim) / h
# a += beta * jump(p) * jump(phi) * (1-I)*dS
# a -= jump(phi) * K * inner(avg(grad( p )), n('+')) * (1-I)*dS
# a -= jump( p ) * K * inner(avg(grad(phi)), n('+')) * (1-I)*dS
# a += beta * (p - pbar) * phi * dirichlet * ds


#####################
# Interface problem #
#####################

# (4.17) - (4.22)
c = d * mu * inner(grad(V_n), tangential) * inner(grad(W_n), tangential) * dx
c += d * mu * inner(grad(V_t), tangential) * inner(grad(W_t), tangential) * dx
c += mu * largeBJterm * V_t * W_t * dx

bW = - d * inner(grad(W_t), tangential) * P * dx
bV = - d * inner(grad(V_t), tangential) * Phi * dx
g = d / (2 * mu) * P * Phi * dx
L_f = 0 if H == as_vector([0,0]) else d * inner(H, as_vector([W_n, W_t])) * ineumann * ds

f_f = 0.5 * jump(trace(ph)) * W_n * dx
e_f = d / (2 * mu) * avg(trace(ph)) * Phi * dx

iBC = dune.ufl.DirichletBC(ispace, [None, vbarf_n, vbarf_t], idirichlet)

#########
# Forms #
#########
# (4.8)
scheme = galerkin([a == L_pm + f_pm + e_pm, BC])
# (4.23)
ischeme = galerkin([(c + bW) + (bV - g) == (L_f - f_f) + (-e_f), iBC])



A      = linearOperator(scheme).as_numpy
rhs    = ph.copy()
zero   = ph.copy()
ph_old = ph.copy()

iA     = linearOperator(ischeme).as_numpy
irhs   = ih.copy()
izero  = ih.copy()
ih_old = ih.copy()

def solve():
    for i in range(100):

        # bulk
        scheme(zero, rhs)
        ph_old.assign(ph)
        ph.as_numpy[:] = linalg.spsolve(A, -rhs.as_numpy)

        phnp = ph_old.as_numpy[:]
        phnp -= ph.as_numpy
        error = np.dot(phnp, phnp)

        # interface
        ischeme(izero, irhs)
        ih_old.assign(ih)
        ih.as_numpy[:] = linalg.spsolve(iA, -irhs.as_numpy)

        ihnp = ih_old.as_numpy[:]
        ihnp -= ih.as_numpy
        ierror = np.dot(ihnp, ihnp)

        print("["+str(i)+"]: errors=", [error, ierror], flush=True)

        if max(error, ierror) < 1e-12:
            break

solve()

gridView.writeVTK("onepstokes-bulk", pointdata={"p": ph, "v": -K/mu*grad(ph)})
    # nonconforming=True, subsampling=order-1)
igridView.writeVTK("onepstokes-interface", pointdata={"P": ih[0], "V_n": ih[1], "V_t": ih[2]})
    # nonconforming=True, subsampling=order-1)
