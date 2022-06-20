#!/usr/bin/env python
# coding: utf-8

# # Moving and adapting
#
# This is an example of how to move the interface and adapt the mesh.

# We implement the finite volume moving mesh method presented in [CMR+18]_. 
# 
# .. [CMR+18] C. Chalons, J. Magiera, C. Rohde, M. Wiebe. A Finite-Volume Tracking Scheme for Two-Phase Compressible Flow. Theory, Numerics and Applications of Hyperbolic Problems I, pp. 309--322, 2018.

# ## Grid creation
# 
# We use the [vertical](grids/vertical.rst) grid file that contains an interface $\Gamma = {0.5} \times [0, 1]$ embedded in a domain $\Omega = [0,1]^2$.
# For this example, we have to construct an _adaptive_ leaf grid view and we will need to obtain the hierarchical grid object.

# In[1]:


from dune.grid import reader
from dune.mmesh import mmesh
from dune.fem.view import adaptiveLeafGridView as adaptive
dim = 2
file = "grids/vertical.msh"
gridView = adaptive( mmesh((reader.gmsh, file), dim) )
hgrid = gridView.hierarchicalGrid
igridView = hgrid.interfaceGrid


# ## Problem
# 
# Let us consider the following transport problem.
# \begin{align}
# u_t + \operatorname{div} f(u) = 0, & \qquad \text{in } \Omega \times [0,T], \\
# u(\cdot, 0) = u_0, & \qquad \text{in } \Omega
# \end{align}
# where
# \begin{align}
# f(u) &= [1,0]^T u, \\
# u_0(x,y) &= (0.5+x) \chi_{x<0.5}.
# \end{align}
# 
# Further, the interface is supposed to move with the transport speed in $f$, i.e. $m = [1,0]^T$.
# \begin{align*}
# \renewcommand{\jump}[1]{[\mskip-5mu[ #1 ]\mskip-5mu]}
# \end{align*}

# In[2]:


from ufl import *
from dune.ufl import Constant

tEnd = 0.4

def speed():
    return as_vector([1.0, 0.0])

def movement(x):
    return as_vector([1.0, 0.0])

def f(u):
    return speed() * u

def u0(x):
    return conditional(x[0] < 0.5, 0.5+x[0], 0.0)

def uexact(x, t):
    return u0( x - t * speed() )


# ## Finite Volume Moving Mesh Method
# 
# We use a Finite Volume Moving Mesh method to keep the discontinuity sharp. It can be formulated by
# \begin{align}
# \int_\Omega (u^{n+1} |det(\Psi)| - u^n) v\ dx + \Delta t \int_\mathcal{F} \big( g(u^{n}, n) - h(u^{n}, n) \big) \jump{v}\ dS = 0
# \end{align}
# where $\Psi := x + \Delta t s$ and s is a linear interpolation of the interface’s vertex movement m on the bulk triangulation.
# 
# The numerical fluxes $g(u, n)$ and $h(u, n)$ are assumed to be consistent with the flux functions $f(u) \cdot n$ and $u s \cdot n$, respectively.

# In[3]:


from dune.fem.space import finiteVolume

space = finiteVolume(gridView)
u = TrialFunction(space)
v = TestFunction(space)

x = SpatialCoordinate(space)
n = FacetNormal(space)

uh = space.interpolate(u0(x), name="uh")
uh_old = uh.copy()


# In[4]:


import numpy as np
from dune.geometry import vertex
from dune.mmesh import edgeMovement

def getShifts():
    mapper = igridView.mapper({vertex: 1})
    shifts = np.zeros((mapper.size, dim))
    for v in igridView.vertices:
        shifts[ mapper.index(v) ] = as_vector(movement( v.geometry.center ))
    return shifts

em = edgeMovement(gridView, getShifts())
t = Constant(0, name="time")

def g(u, n):
    sgn = inner(speed(), n('+'))
    return inner( conditional( sgn > 0, f( u('+') ), f( u('-') ) ), n('+') )

def gBnd(u, n):
    sgn = inner(speed(), n)
    return inner( conditional( sgn > 0, f(u), f(uexact(x, t)) ), n )

def h(u, n):
    sgn = inner(em('+'), n('+'))
    return conditional( sgn > 0, sgn * u('+'), sgn * u('-') )


# In[5]:


from dune.fem.scheme import galerkin

tau = Constant(0, name="tau")
detPsi = abs(det(nabla_grad(x + tau * em)))

a = (u * detPsi - uh_old) * v * dx
a += tau * (g(uh_old, n) - h(uh_old, n)) * jump(v) * dS
a += tau * gBnd(uh_old, n) * v * ds

scheme = galerkin([a == 0], solver=("suitesparse","umfpack"))


# Compute CFL time step using cell edge lengths
l, ll = TrialFunction(space), TestFunction(space)
b = (l - MinCellEdgeLength(space)) * ll * dx
lscheme = galerkin([b == 0], solver=("suitesparse","umfpack"))
lh = space.interpolate(0, name="lh")

def dtCFL():
    lscheme.solve(lh)
    lMin = np.min(lh.as_numpy)
    return 0.5 * lMin / speed()[0]


# ## Timeloop without adaptation
# 
# For comparison, we run the finite volume scheme once without adaptation of the mesh.

# In[6]:


from time import time
from dune.fem.function import integrate
from dune.fem.plotting import plotPointData as plot
import matplotlib.pyplot as plt
plots = 4
fig, axs = plt.subplots(1, plots, figsize=(3*plots,3))
def shouldPlot(a):
    t_a = a * tEnd / (plots-1)
    return t.value > t_a and t.value - t_a <= tau.value

runtime = 0

uh.interpolate(u0(x))

def L2Error(uh):
    u = uexact(x, t)
    return sqrt(integrate(gridView, dot(uh-u, uh-u), order=5))

l2s = []

# Disable edge movement
em.interpolate([0,0])

tau.value = dtCFL()
print(f"dt = {tau.value:.4f}")

i = 0
t.assign(0)
while t.value < tEnd:
    runtime -= time()

    t.value += tau.value

    uh_old.assign(uh)
    scheme.solve(target=uh)
    l2s += [(t.value, L2Error(uh))]

    runtime += time()

    i += 1
    for a in range(plots):
        if shouldPlot(a):
            plot(uh, figure=(fig, axs[a]), clim=[0,1], colorbar=None)

print(f"Runtime: {runtime:.3f}s")

L2 = L2Error(uh)
print(f"L2-Error: {L2:.3f}")


# ## Timeloop with adaptation
# 
# We need to set the `"fem.adaptation.method"` parameter to `"callback"` in order to use the non-hierarchical adaptation strategy of Dune-MMesh. Then, within the time loop, we can adapt the mesh according to the following strategy.

# In[7]:


from dune.fem import parameter, adapt
parameter.append( { "fem.adaptation.method": "callback" } )
fig.clear()
fig, axs = plt.subplots(1, plots, figsize=(3*plots,3))
runtimeAdapted = 0
dtmin = 1e3

uh.interpolate(u0(x))

l2sAdapted = []

i = 0
t.assign(0)
while t.value < tEnd:
    runtimeAdapted -= time()

    tau.value = dtCFL()
    dtmin = min(dtmin, tau.value)

    t.value += tau.value
    
    mark = hgrid.markElements()
    ensure = hgrid.ensureInterfaceMovement(getShifts()*tau.value)
    if mark or ensure:
        adapt([uh])
    
    shifts = getShifts()
    em.assign(edgeMovement(gridView, shifts))

    uh_old.assign(uh)
    scheme.solve(target=uh)

    hgrid.moveInterface(shifts*tau.value)
    l2sAdapted += [(t.value, L2Error(uh))]

    runtimeAdapted += time()

    i += 1
    for a in range(plots):
        if shouldPlot(a):
            plot(uh, figure=(fig, axs[a]), clim=[0,1], colorbar=None)


print(f"dt_min = {dtmin:.4f}")
print(f"Runtime: {runtimeAdapted:.3f}s")

L2Adapted = L2Error(uh)
print(f"L2-Error: {L2Adapted:.3f}")

print(f"\nRuntime factor: {runtimeAdapted / runtime:.2f}x")
print(f"Error improvement: {L2 / L2Adapted:.2f}x")

import matplotlib.pyplot as plt
fig = plt.figure()
plt.plot([l2[0] for l2 in l2s], [l2[1] for l2 in l2s], label="Without adaptation")
plt.plot([l2[0] for l2 in l2sAdapted], [l2[1] for l2 in l2sAdapted], label="Adapted")
plt.legend(loc='best')
plt.title('L2-Error')
plt.show()
