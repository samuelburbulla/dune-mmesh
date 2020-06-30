## @example onep.py
#  Single-phase flow in fractured porous media with IPDG

import io
from dune.grid import reader
from dune.mmesh import mmesh, trace, skeleton
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

dim = 2
file = "vertical2d.msh"

# MMesh
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

from ufl import *
import dune.ufl
from dune.fem import parameter, adapt
from dune.fem.function import integrate, uflFunction
from dune.fem.space import dglagrange
from dune.fem.scheme import galerkin
from dune.ufl import DirichletBC
from math import log
parameter.append({"fem.verboserank": 0, "fem.adaptation.method": "callback"})
solverParameters =\
   {"newton.tolerance": 1e-8,
    "newton.linear.tolerance": 1e-12,
    "newton.linear.preconditioning.method": "jacobi",
    "newton.linear.maxiterations": 10000,
    "preconditioning.method": "jacobi",
    "newton.verbose": True,
    "newton.linear.verbose": False}

errors = []
eocs = []
for i in range(3):
    # Bulk problem
    order = 1
    space = dglagrange(gridView, order=order)
    p = TrialFunction(space)
    phi = TestFunction(space)

    x = SpatialCoordinate(space)
    n = FacetNormal(space)

    K = as_matrix([[1, 0], [0, 1]])
    pexact = conditional( x[0]<0.5, sin(4*x[0])*cos(pi*x[1]), cos(4*x[0])*cos(pi*x[1]) )
    q = (16.+pi**2) * pexact
    g = pexact

    mu = 1000 # TODO: penalty parameter

    space_gamma = dglagrange(igridView, order=order)
    one = space_gamma.interpolate(1, name="one")
    interface = skeleton(one)('+')

    L = q * phi * dx
    L += mu * g * phi * ds
    L -= g * dot(dot(grad(phi), K), n) * ds

    B = dot(dot(grad(p), K), grad(phi)) * dx
    B += mu * jump(p) * jump(phi) * (1-interface) * dS
    B -= dot(avg(dot(grad(p), K)), n('+')) * jump(phi) * (1-interface) * dS
    B -= jump(p) * dot(avg(dot(grad(phi), K)), n('+')) * (1-interface) * dS
    B += mu * p * phi * ds
    B -= p * dot(dot(grad(phi), K), n) * ds
    B -= phi * dot(dot(grad(p), K), n) * ds


    # Interface problem
    space_gamma = dglagrange(igridView, order=order)
    p_gamma = TrialFunction(space_gamma)
    phi_gamma = TestFunction(space_gamma)

    x_gamma = SpatialCoordinate(space_gamma)
    n_gamma = FacetNormal(space_gamma)

    d = 1 # TODO d != 1 is not correct
    K_gamma = as_matrix([[2, 0], [0, 1]])
    xi = 3./4.
    p_gammaexact = 3./4.*(cos(2.)+sin(2.))*cos(pi*x_gamma[1])
    q_gamma = 4./3.*p_gammaexact*(4. + 3./4.* d**2 * pi**2)
    g_gamma = p_gammaexact

    mu_gamma = mu

    L_gamma = d * q_gamma * phi_gamma * dx
    L_gamma += mu_gamma * d**2 * g_gamma * phi_gamma * ds
    L_gamma -= d**2 * g_gamma * dot(dot(grad(d * phi_gamma), K_gamma), n_gamma) * ds

    B_gamma = d * dot(dot(grad(d*p_gamma), K_gamma), grad(d*phi_gamma)) * dx
    B_gamma += d * mu_gamma * d * jump(p_gamma) * jump(phi_gamma) * dS
    B_gamma -= d * jump(phi_gamma) * dot(avg(dot(grad(d *  p_gamma ), K_gamma)), n_gamma('+')) * dS
    B_gamma -= d *  jump(p_gamma)  * dot(avg(dot(grad(d * phi_gamma), K_gamma)), n_gamma('+')) * dS
    B_gamma += mu_gamma * d**2 * p_gamma * phi_gamma * ds
    B_gamma -= d**2 * p_gamma * dot(dot(grad(d * phi_gamma), K_gamma), n_gamma) * ds
    B_gamma -= d**2 * phi_gamma * dot(dot(grad(d * p_gamma), K_gamma), n_gamma) * ds


    ph = space.interpolate(0, name="pressure")
    ph_gamma = space_gamma.interpolate(0, name="pressure")

    # Coupling
    igvNormal = igridView.normal
    nK_gamman = dot( dot(K_gamma, igvNormal), igvNormal)
    skeletonnK_gamman = skeleton(ph_gamma)
    beta = 4. * nK_gamman / (2. * xi - 1)
    skeletonbeta = 4. * skeletonnK_gamman / (2. * xi - 1)
    skeletonp_gamma = skeleton(ph_gamma)
    C = skeletonnK_gamman('+') * jump(p) * jump(phi) * dS
    C += skeletonbeta * ( avg(p) - skeletonp_gamma ) * dS

    tracep = trace(ph)
    C_gamma = nK_gamman * jump(tracep) * dx
    C_gamma += beta * ( avg(tracep) - p_gamma ) * dx


    # Scheme
    scheme = galerkin([B == L], solver="cg", parameters=solverParameters)
    scheme_gamma = galerkin([B_gamma == L_gamma], solver="cg", parameters=solverParameters)

    def solve():
        for i in range(100):
          print("iteration", i)

          ph_old = ph.copy()
          ph_gamma_old = ph_gamma.copy()

          print("  solve bulk")
          scheme.solve(target=ph)

          print("  solve interface")
          scheme_gamma.solve(target=ph_gamma)

          if integrate(gridView, dot(ph-ph_old, ph-ph_old), order=5) + \
             integrate(igridView, dot(ph_gamma-ph_gamma_old, ph_gamma-ph_gamma_old), order=5) < 1e-14:
                break


    print("\nEOC", i, "\n")
    solve()

    errorbulk = integrate(gridView, dot(ph-pexact, ph-pexact), order=5)
    print("  error bulk", errorbulk)
    errorinterface = integrate(igridView, dot(ph_gamma-p_gammaexact, ph_gamma-p_gammaexact), order=5)
    print("  error interface", errorinterface)
    error = errorbulk + errorinterface
    print("  error bulk + interface", error)
    errors += [error]
    if i > 0:
        eocs += [ log( errors[i-1] / errors[i] ) / log(2) ]


    gridView.writeVTK("onep-bulk", pointdata={"p": ph}, nonconforming=True, subsampling=order-1)
    igridView.writeVTK("onep-interface", pointdata={"p": ph_gamma}, nonconforming=True, subsampling=order-1)

    gridView.hierarchicalGrid.globalRefine(dim)
    igridView.hierarchicalGrid.globalRefine(dim-1)

print(errors)
print(eocs)
