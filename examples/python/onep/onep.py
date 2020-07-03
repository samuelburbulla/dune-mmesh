# <markdowncell>
# ## Single-phase flow in fractured porous media with IPDG
# Credit: Maximilian Hoerl
#
# \begin{align}
# \newcommand{\llbracket}{[\mskip-5mu[}
# \newcommand{\rrbracket}{]\mskip-5mu]}
# \newcommand{\jump}[1]{\llbracket #1 \rrbracket}
# \newcommand{\ldblbrace}{\{\mskip-5mu\{}
# \newcommand{\rdblbrace}{\}\mskip-5mu\}}
# \newcommand{\avg}[1]{\ldblbrace #1 \rdblbrace}
# \end{align}
#
#
# ## Model equations
# \begin{align}
# -\nabla\cdot\left(\mathbf{K}_i\nabla p_i\right) &= q_i &\qquad &\text{in } \Omega_i, \quad &i&=1,2 ,\\
# -\nabla_\parallel \cdot \left( d\mathbf{K}^\Gamma_\parallel \nabla_\parallel \left( d p^\Gamma\right) \right) &= q^\Gamma - \jump{\mathbf{K}\nabla p} &\qquad &\text{in }\Gamma ,\\
# -\avg{\mathbf{K} \nabla p} \cdot {n} &= K^\Gamma_\perp \left( p_1 - p_2 \right)  &\qquad &\text{on }\Gamma ,\\
# -\jump{\mathbf{K} \nabla p} &= \beta^\Gamma \left( \avg{p} - p^\Gamma  \right) &\qquad &\text{on } \Gamma ,\\
# p_i &= g_i &\qquad &\text{on } \rho_i , &i&=1,2 , \\
# p^\Gamma &= g^\Gamma &\qquad &\text{on } \partial\Gamma, \\
# \end{align}
# where $\beta = \tfrac{4 K_\perp^\Gamma}{2 \xi -1}$.
#
#
# ## IPDG formulation
#
# Find $\left( p_h , p_h^\Gamma \right) \in  \Phi_h^{PM} \times \Phi_h^\Gamma $, s.t.
# \begin{align}
# \mathcal{B}_{DG} \left( \left( p_h  , p_h^\Gamma \right) , \left( \phi_h , \phi_h^\Gamma \right)\right) = \mathcal{L}_{DG} \left( \phi_h , \phi_h^\Gamma \right)
# \end{align}
# for all $\left( \phi_h , \phi_h^\Gamma \right) \in  \Phi_h^{PM} \times \Phi_h^\Gamma$, where

# \begin{align}
# \begin{split}
# \mathcal{B}_{DG} \left( \left( p_h  , p_h^\Gamma \right) , \left( \phi_h , \phi_h^\Gamma \right)\right) &:= \mathcal{B}_{DG}^{PM} \left( p_h , \phi_h \right) + \mathcal{B}^\Gamma_{DG} \left( p_h^\Gamma , \phi_h^\Gamma \right) + \mathcal{K}_{DG}  \left( \left( p_h  , p_h^\Gamma \right) , \left( \phi_h , \phi_h^\Gamma \right)\right), \\
# \mathcal{L}_{DG} \left( \phi_h , \phi_h^\Gamma \right) &:= \mathcal{L}^{PM}_{DG} \left( \phi_h \right) + \mathcal{L}^\Gamma_{DG} \left( \phi_h^\Gamma \right) ,
# \end{split}
# \end{align}
# and
# \begin{align}
# \mathcal{B}^{PM}_{DG} \left( p_h , \phi_h \right) &:= \int_\Omega \mathbf{K} \nabla_h p_h \cdot \nabla_h \phi_h \,{d} V \\
# &\quad + \int_{\mathcal{F}^{\circ}\setminus\Gamma} \left( \mu \jump{p_h}\cdot\jump{\phi_h} - \avg{\mathbf{K} \nabla_h p_h} \cdot \jump{\phi_h} - \jump{p_h} \cdot \avg{\mathbf{K} \nabla_h \phi_h }\right) \,{d}\sigma \\
# &\quad + \int_{\partial\Omega} \mu  p_h \phi_h \,{d}\sigma -  \int_{\partial\Omega} \left( p_h  \mathbf{K} \nabla_h \phi_h + \phi_h \mathbf{K} \nabla_h p_h \right) \cdot {d}{\sigma}, \\
# \mathcal{L}^{PM}_{DG} \left( \phi_h\right) &:= \int_\Omega q \phi_h \,{d} V +\int_{\partial\Omega} \mu  g \phi_h \,{d}\sigma - \int_{\partial\Omega} g \mathbf{K} \nabla_h \phi_h \cdot {d}{\sigma}, \\
# \mathcal{B}^\Gamma_{DG} \left( p_h^\Gamma , \phi_h^\Gamma \right) &:= \int_\Gamma \mathbf{K}_\parallel^\Gamma \nabla_{\parallel h}\left( d p_h^\Gamma \right)  \cdot \nabla_{\parallel h} \left( d \phi_h^\Gamma \right) \,{d}\sigma
# + \int_{\mathcal{E}_\Gamma^\circ} d \Big[ \mu^\Gamma  d \jump{ p_h^\Gamma} \cdot \jump{\phi_h^\Gamma} \\
# &\quad - \avg{\mathbf{K}_\parallel^\Gamma \nabla_{\parallel h} \left( d p_h^\Gamma \right) } \cdot \jump{\phi_h^\Gamma} - \jump{p_h^\Gamma} \cdot \avg{\mathbf{K}_\parallel^\Gamma \nabla_{\parallel h}\left( d \phi_h^\Gamma\right) } \Big] \,{d}r \\
# &\quad + \int_{\partial\Gamma} \mu^\Gamma d^2 p_h^\Gamma \phi_h^\Gamma \,{d} r - \int_{\partial \Gamma} d \left[ p_h^\Gamma \mathbf{K}_\parallel^\Gamma \nabla_{\parallel h} \left( d \phi_h^\Gamma \right) + \phi_h^\Gamma \mathbf{K}_\parallel^\Gamma \nabla_{\parallel h} \left( d p_h^\Gamma\right)\right] \cdot {d}{r},\\
# \mathcal{L}_{DG}^\Gamma \left( \phi_h^\Gamma \right) &:= \int_\Gamma q^\Gamma \phi_h^\Gamma \,{d}\sigma + \int_{\partial\Gamma} \mu^\Gamma d^2 g^\Gamma \phi_h^\Gamma \,{d} r - \int_{\partial\Gamma} d g^\Gamma \mathbf{K}_\parallel^\Gamma \nabla_{\parallel h} \left( d \phi_h^\Gamma \right) \cdot {d} {r}, \\
# \mathcal{K}_{DG} \left( \left( p_h , p_h^\Gamma \right) , \left( \phi_h , \phi_h^\Gamma \right) \right) & := \int_\Gamma K_\perp^\Gamma \jump{p_h} \cdot \jump{\phi_h} \,{d} \sigma +  \int_\Gamma \beta^\Gamma \left( \avg{p_h} - p_h^\Gamma \right) \left( \avg{\phi_h} - \phi_h^\Gamma \right) \,{d}\sigma .
# \end{align}
#
# Penalty
# \begin{align}
# \mu := \mu_0 \bar{\lambda}_{T} \frac{\left( k_T + 1 \right) \left( k_T + n\right)}{h_T},  \mu_0 > 0
# \end{align}


# ## Test problem
#
# $\Omega_1 = [0,0.5] \times [0,1], \Omega_2 = [0.5,1] \times [0,1], \Gamma = \{ 0.5 \} \times [0,1]$

# Exact solutions
# \begin{align*}
# p_1 \left( {x} \right) &=  \sin\left( 4x_1 \right) \cos\left( \pi x_2 \right) ,\\
# p_2 \left( {x} \right) &=  \cos\left( 4x_1 \right) \cos\left( \pi x_2 \right) ,\\
# p^\Gamma \left( {s}\right) &= \tfrac{3}{4} \left( \cos\left( 2\right) + \sin\left( 2\right) \right) \cos\left( \pi x_2 \right)
# \end{align*}
#
# $\mathbf{K} = \mathbf{I}$, $\mathbf{K}_\parallel^\Gamma = \mathbf{I}, K_\perp^\Gamma = 2$, $d>0$, $\xi = \tfrac{3}{4}$
#
# Source terms
# \begin{align*}
# q_1 \left( {x} \right) &= \left( 16 + \pi^2 \right) \sin\left( 4x_1 \right) \cos\left( \pi x_2 \right) ,\\
# q_2 \left( {x} \right) &=  \left( 16 + \pi^2 \right) \cos\left( 4x_1 \right) \cos\left( \pi x_2 \right),\\
# q^\Gamma \left( {s} \right) &= \left( \cos\left( 2\right) + \sin\left( 2 \right) \right) \left( 4 + \tfrac{3}{4} d^2 \pi^2 \right) \cos\left( \pi s_2 \right).
# \end{align*}

# <codecell>

import vertical
from dune.grid import reader
from dune.mmesh import mmesh, trace, skeleton, domainMarker
import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

from ufl import *
import dune.ufl
from dune.fem import parameter, globalRefine, adapt
from dune.fem.function import integrate, uflFunction
from dune.fem.space import dglagrange, finiteVolume
from dune.fem.scheme import galerkin
from dune.fem.view import adaptiveLeafGridView as adaptive
from dune.ufl import DirichletBC
from math import log
parameter.append({"fem.verboserank": 0,
                  "fem.adaptation.method": "callback"})
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
for i in range(4):

    file = "vertical.msh"
    vertical.create(file, h=0.5*2**-i)

    dim = 2
    gridView = mmesh((reader.gmsh, file), dim)
    igridView = gridView.hierarchicalGrid.interfaceGrid

    fvspace = finiteVolume(gridView)
    dm = fvspace.interpolate(domainMarker(gridView), name="dm")

    # Bulk problem
    order = 2
    space = dglagrange(gridView, order=order)
    p = TrialFunction(space)
    phi = TestFunction(space)

    x = SpatialCoordinate(space)
    n = FacetNormal(space)

    K = as_matrix([[1, 0], [0, 1]])

    pexact = (1-dm) * sin(4*x[0])*cos(pi*x[1]) + (dm) * cos(4*x[0])*cos(pi*x[1])
    q = (16.+pi**2) * pexact
    g = pexact

    mu0 = 10
    lambdaMax = 1
    h = MaxFacetEdgeLength(space.cell())
    mu = mu0 * lambdaMax * (order + 1) * (order + dim) / h

    space_gamma = dglagrange(igridView, order=order)
    one = space_gamma.interpolate(1, name="one")
    I = avg(skeleton(one))

    L = q * phi * dx
    L += mu * g * phi * ds
    L -= g * dot(dot(grad(phi), K), n) * ds

    B = dot(dot(grad(p), K), grad(phi)) * dx
    B += mu * jump(p) * jump(phi) * (1-I)*dS
    B -= jump(phi) * dot(avg(dot(grad( p ), K)), n('+')) * (1-I)*dS
    B -= jump( p )  * dot(avg(dot(grad(phi), K)), n('+')) * (1-I)*dS
    B += mu * p * phi * ds
    B -= p * dot(dot(grad(phi), K), n) * ds
    B -= phi * dot(dot(grad(p), K), n) * ds


    # Interface problem
    space_gamma = dglagrange(igridView, order=order)
    p_gamma = TrialFunction(space_gamma)
    phi_gamma = TestFunction(space_gamma)

    x_gamma = SpatialCoordinate(space_gamma)
    n_gamma = FacetNormal(space_gamma)

    d = 0.01
    K_gamma = as_matrix([[2, 0], [0, 1]])
    xi = 3./4.
    p_gammaexact = 3./4.*(cos(2.)+sin(2.))*cos(pi*x_gamma[1])
    q_gamma = 4./3.*p_gammaexact*(4. + 3./4.* d**2 * pi**2)
    g_gamma = p_gammaexact

    lambdaMax_gamma = 1
    h_gamma = MaxFacetEdgeLength(space_gamma.cell())
    mu_gamma = mu0 * lambdaMax_gamma * (order + 1) * (order + dim) / h_gamma

    L_gamma = q_gamma * phi_gamma * dx
    L_gamma += mu_gamma * d**2 * g_gamma * phi_gamma * ds
    L_gamma -= d * g_gamma * dot(dot(grad(d * phi_gamma), K_gamma), n_gamma) * ds

    B_gamma = dot(dot(grad(d*p_gamma), K_gamma), grad(d*phi_gamma)) * dx
    B_gamma += d * mu_gamma * d * jump(p_gamma) * jump(phi_gamma) * dS
    B_gamma -= d * jump(phi_gamma) * dot(avg(dot(grad(d *  p_gamma ), K_gamma)), n_gamma('+')) * dS
    B_gamma -= d * jump( p_gamma ) * dot(avg(dot(grad(d * phi_gamma), K_gamma)), n_gamma('+')) * dS
    B_gamma += mu_gamma * d**2 * p_gamma * phi_gamma * ds
    B_gamma -= d * p_gamma * dot(dot(grad(d * phi_gamma), K_gamma), n_gamma) * ds
    B_gamma -= d * phi_gamma * dot(dot(grad(d * p_gamma), K_gamma), n_gamma) * ds


    ph = space.interpolate(0, name="pressure")
    ph_gamma = space_gamma.interpolate(0, name="pressure")

    # Coupling
    skelK_perp = dot(dot(K_gamma, n('+')), n('+'))
    skelBeta = 4. * skelK_perp / (2. * xi - 1)

    C = skelK_perp * jump(p) * jump(phi) * I*dS
    C += skelBeta * ( avg(p) * avg(phi) - skeleton(ph_gamma)('+') * avg(phi) ) * I*dS

    inormal = igridView.normal
    K_perp = dot(dot(K_gamma, inormal), inormal)
    beta = 4. * K_perp / (2. * xi - 1)

    C_gamma = beta * ( -avg(trace(ph)) * phi_gamma + p_gamma * phi_gamma ) * dx


    # Scheme
    scheme = galerkin([B + C == L], solver='cg', parameters=solverParameters)
    scheme_gamma = galerkin([B_gamma + C_gamma == L_gamma], solver='cg', parameters=solverParameters)

    def solve():
        for i in range(100):
          print("# Iteration", i)

          ph_old = ph.copy()
          ph_gamma_old = ph_gamma.copy()

          print("Solve bulk")
          scheme.solve(target=ph)

          print("Solve interface")
          scheme_gamma.solve(target=ph_gamma)

          if integrate(gridView, dot(ph-ph_old, ph-ph_old), order=order) + \
             integrate(igridView, dot(ph_gamma-ph_gamma_old, ph_gamma-ph_gamma_old), order=order) < 1e-14:
                break


    print("\nEOC", i, "\n")
    solve()

    errorbulk = sqrt(integrate(gridView, dot(ph-pexact, ph-pexact), order=order))
    print("  error bulk", errorbulk)
    errorinterface = sqrt(integrate(igridView, dot(ph_gamma-p_gammaexact, ph_gamma-p_gammaexact), order=order))
    print("  error interface", errorinterface)
    error = errorbulk + errorinterface
    print("  error bulk + interface", error)
    errors += [error]
    if i > 0:
        eocs += [ log( errors[i-1] / errors[i] ) / log(2) ]


    gridView.writeVTK("onep-bulk-"+str(i), pointdata={"p": ph, "exact": pexact},
        nonconforming=True, subsampling=order-1)
    igridView.writeVTK("onep-interface-"+str(i), pointdata={"p": ph_gamma, "exact": p_gammaexact},
        nonconforming=True, subsampling=order-1)

print(errors)
print("EOC:", eocs)
