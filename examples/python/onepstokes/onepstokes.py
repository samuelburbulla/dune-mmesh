import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

import horizontal
from dune.grid import reader
from dune.mmesh import mmesh, trace, skeleton, domainMarker, interfaceIndicator

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
horizontal.create(file, h=0.1)

dim = 2
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

##########
# Spaces #
##########
order = 1
space = dglagrange(gridView, order=order)
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
ih = ispace.interpolate([0,0,0], name="isolution")
Ph   = ih[0]
Vh_n = ih[1]
Vh_t = ih[2]

x = SpatialCoordinate(space)
n = FacetNormal(space)
ix = SpatialCoordinate(ispace)
in = FacetNormal(ispace)
I = interfaceIndicator(igridView)

################
# Bulk problem #
################

mu = 1
K = 1
q = 0
g = 0
vbar = 0

# (4.4) - (4.7)
a = K / mu * grad(p) * grad(phi) * dx
a += d / (2 * mu) * avg(p * phi) * I*dS
f_pm = skeleton(Vh_n) * jump(phi) * I*dS
e_pm = d / (2 * mu) * skeleton(Ph) * avg(phi) * I*dS
L_pm = q * phi * dx - vbar * phi * ds

# consistency terms (TODO)
beta0 = 10
h = MaxFacetEdgeLength(space.cell())
beta = beta0 * K * (order + 1) * (order + dim) / h
a += beta * jump(p) * jump(phi) * (1-I)*dS
# a -= jump(phi) * dot(avg(dot(grad( p ), K)), n('+')) * (1-I)*dS
# a -= jump( p )  * dot(avg(dot(grad(phi), K)), n('+')) * (1-I)*dS


#####################
# Interface problem #
#####################
d = 0.01
alpha_BJ = 1
q_gamma = 0
g_gamma = 0
largeBJterm = 1 # TODO
F = as_vector([0,0])
W = as_vector([0,0])

# (4.17) - (4.22)
c = d * mu * grad(V_n) * grad(W_n) * dx
c += d * mu * grad(V_t) * grad(W_t) * dx
c += mu * largeBJterm * V_t * W_t * dx

b = - d * div(V_t) * Phi * dx
g = d / (2 * mu) * P * Phi * dx
L_f = d * inner(F, as_vector([W_n, W_t])) * dx + d * inner(H, as_vector([W_n, W_t])) * ds

f_f = 0.5 * jump(trace(ph)) * W_n * dx
e_f = d / (2 * mu) * avg(trace(ph)) * Phi * dx



#########
# Forms #
#########
# (4.8)
scheme = galerkin([a == L_pm + f_pm + e_pm])
# (4.23)
ischeme = galerkin([(c + b) + (b - g) == (L_f - f_f) + (-e_f)])



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

gridView.writeVTK("onepstokes-bulk-"+str(i), pointdata={"p": ph},
    nonconforming=True, subsampling=order-1)
igridView.writeVTK("onepstokes-interface-"+str(i), pointdata={"P": Ph, "V_n": Vh_n, "V_t": Vh_t},
    nonconforming=True, subsampling=order-1)
