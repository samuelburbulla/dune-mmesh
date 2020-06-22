import io
import math
from dune.grid import reader
from dune.mmesh import mmesh, domainMarker
# from dune.femdg import createLimiter

dim = 2
file = "../grids/cube2d.msh"

gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

import ufl
from dune.fem.function import integrate
from dune.fem.space import dglagrange, finiteVolume
from dune.fem.scheme import galerkin
from dune.ufl import DirichletBC, Constant, Space

dt = 0.001
T = 1

indicator = domainMarker(gridView)

Kb = 1
Kf = 1
K = indicator*Kf + (1-indicator)*Kb

# solve laplace
space = dglagrange(gridView, dimRange=2, order=1)
trial = ufl.TrialFunction(space)
test = ufl.TestFunction(space)
p = trial[0]
phi = test[0]
s = trial[1]
psi = test[1]

solution_old = space.interpolate([0,0], name="solution_old")

tau = Constant(dt, name="timeStep")
x = ufl.SpatialCoordinate(space)

dirichlet = ufl.conditional(x[0] < 1e-6, 1.0, 0.0) + ufl.conditional(x[0] > 1.0 - 1e-6, 1.0, 0.0)
p_D = ufl.conditional(x[0] < 0.5, 1.0, 0.0)
s_D = ufl.conditional(x[0] < 0.5, 1.0, 0.0)

dBulk_p = K * ufl.grad(p)
dBulk_s = s * dBulk_p
a = ufl.inner(dBulk_p, ufl.grad(phi)) * ufl.dx
a += (s - solution_old[1]) * psi * ufl.dx + tau * ufl.inner(dBulk_s, ufl.grad(psi)) * ufl.dx

# DG terms
penalty   = 5 * (space.order * ( space.order + 1 ))
beta     = Constant(penalty, name="penalty")
def sMax(a): return ufl.max_value(a('+'), a('-'))
cell     = space.cell()
n         = ufl.FacetNormal(cell)
hT        = ufl.MaxCellEdgeLength(cell)
he        = ufl.avg( ufl.CellVolume(cell) ) / ufl.FacetArea(cell)
heBnd     = ufl.CellVolume(cell) / ufl.FacetArea(cell)
k         = ufl.dot(K*n,n)
lambdaMax = k('+')*k('-')/ufl.avg(k)
def wavg(z): return (k('-')*z('+')+k('+')*z('-'))/(k('+')+k('-'))
# Penalty terms (including dirichlet boundary treatment)
l_n0 = Constant(0.5, name="l_n0")
l_w0 = Constant(0.5, name="l_w0")
penalty_p = [beta*lambdaMax*sMax(l_n0+l_w0), beta*k*(l_n0+l_w0)]
penalty_s = [beta*lambdaMax*sMax(l_n0), beta*k*(l_n0)]
a += penalty_p[0]/he * ufl.jump(p)*ufl.jump(phi) * ufl.dS
a += penalty_s[0]/he * ufl.jump(s)*ufl.jump(psi) * ufl.dS
a += penalty_p[1]/heBnd * (p-p_D) * phi * dirichlet * ufl.ds
a += penalty_s[1]/heBnd * (s-s_D) * psi * dirichlet * ufl.ds
# Consistency terms
a -= ufl.inner(wavg(dBulk_p),n('+')) * ufl.jump(phi) * ufl.dS
a -= ufl.inner(wavg(dBulk_s),n('+')) * ufl.jump(psi) * ufl.dS
a -= ufl.inner(dBulk_p,n) * phi * dirichlet * ufl.ds
a -= ufl.inner(dBulk_s,n) * psi * dirichlet * ufl.ds

newtonParameters = {"tolerance": 1e-8,
                    "verbose": "true", "linear.verbose": "false",
                    "linabstol": 1e-8, "reduction": 1e-8}
scheme = galerkin([a==0], solver=("suitesparse","umfpack"), parameters={"newton." + k: v for k, v in newtonParameters.items()})

solution = space.interpolate(ufl.as_vector([0,0]), name="solution")
# limiter = createLimiter( space, limiter="scaling" )
# tmp = solution.copy()
# def limit(target):
#     tmp.assign(target)
#     limiter(tmp,target)

intermediate = solution.copy()
def solveStep():
    solution_old.assign(solution)
    iter = 0
    while True:
        intermediate.assign(solution)
        scheme.solve(target=solution)
        # limit(solution)
        iter += 1
        print(" iter =", iter)
        if integrate(gridView, solution[1] - intermediate[1], 5 ) < 1e-8**2 or iter > 20:
            break

print("Start time loop...")
t = 0
step = 0

vspace = dglagrange(gridView, dimRange=3, order=space.order-1)
velocity = vspace.project([-K*ufl.grad(solution[0])[0], -K*ufl.grad(solution[0])[1], 0], name="velocity")
gridView.writeVTK("solution-"+str(step), pointdata={"pressure": solution[0], "saturation": solution[1], "velocity": velocity}, celldata=[indicator])

while t < T:
    solveStep()
    step += 1

    velocity = vspace.project([-K*ufl.grad(solution[0])[0], -K*ufl.grad(solution[0])[1], 0], name="velocity")
    gridView.writeVTK("solution-"+str(step), pointdata={"pressure": solution[0], "saturation": solution[1], "velocity": velocity}, celldata=[indicator])
    igridView.writeVTK("interface-"+str(step))

    print("step =", step)
    t += dt
