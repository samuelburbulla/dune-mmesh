from dune.grid import yaspGrid, reader, cartesianDomain
from dune.alugrid import aluConformGrid as alu
from dune.mmesh import mmesh
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin as scheme
from dune.fem.operator import galerkin as op
from dune.ufl import DirichletBC
from ufl import *

def warn(grid):
    try:
        return grid.hierarchicalGrid._cartesianConstructionWithIds
    except AttributeError:
        return False

def test(lgv, domain):
    grid = lgv(domain, dimgrid=2)
    print("if only id is '1' need warning?\t\t", warn(grid))
    spc = lagrange(grid)
    u = TrialFunction(spc)
    v = TestFunction(spc)
    eqn = u * v * dx == 0
    dbc = DirichletBC(spc, 0)
    dbc1 = DirichletBC(spc, 0, 1)
    dbc2 = DirichletBC(spc, 0, 2)
    print("only id=1:")
    scheme([eqn, dbc1])
    print("no id:")
    scheme([eqn, dbc])
    print("id=1 and 2:")
    scheme([eqn, dbc1, dbc2])

print("alu, old boundary treatment", end=" -> ")
test( alu, cartesianDomain([0,0], [1,1], [10,10], boundary=False) )
print("alu, new boundary treatment", end=" -> ")
test( alu, cartesianDomain([0,0], [1,1], [10,10], boundary=True) )
print("alu, default treatment", end=" -> ")
test( alu, cartesianDomain([0,0], [1,1], [10,10]) )

print("yasp, old boundary treatment", end=" -> ")
test( yaspGrid, cartesianDomain([0,0], [1,1], [10,10], boundary=False) )
print("yasp, new boundary treatment", end=" -> ")
test( yaspGrid, cartesianDomain([0,0], [1,1], [10,10], boundary=True) )
print("yasp, default treatment", end=" -> ")
test( yaspGrid, cartesianDomain([0,0], [1,1], [10,10]) )

print("mmesh, old boundary treatment", end=" -> ")
test( mmesh, cartesianDomain([0,0], [1,1], [10,10], boundary=False) )
print("mmesh, new boundary treatment", end=" -> ")
test( mmesh, cartesianDomain([0,0], [1,1], [10,10], boundary=True) )
print("mmesh, default treatment", end=" -> ")
test( mmesh, cartesianDomain([0,0], [1,1], [10,10]) )
