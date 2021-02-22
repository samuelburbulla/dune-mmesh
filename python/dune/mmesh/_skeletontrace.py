import logging, traceback
logger = logging.getLogger(__name__)

import hashlib

from dune.generator.generator import SimpleGenerator


################################################################################
# Skeleton function
################################################################################
def skeleton(interfaceFunction, grid=None):
    if grid == None:
        grid = interfaceFunction.space.grid.hierarchicalGrid.bulkGrid

    includes = ["dune/mmesh/misc/pyskeletontrace.hh"]
    includes += interfaceFunction._includes + grid._includes
    generator = SimpleGenerator("SkeletonGF", "Dune::Fem")

    typeName = "Dune::Fem::SkeletonGF< " + grid._typeName + ", " + interfaceFunction._typeName + " >"
    moduleName = "skeleton_" + hashlib.md5(typeName.encode('utf8')).hexdigest()
    cls = generator.load(includes, typeName, moduleName)
    skeleton = cls.SkeletonGF(grid, interfaceFunction)

    import dune.ufl
    skeleton = dune.ufl.GridFunction(skeleton)

    interfaceFunction.skeleton = skeleton
    return skeleton
################################################################################


################################################################################
# Trace function
################################################################################
def trace(bulkFunction, igrid=None):
    if igrid == None:
      igrid = bulkFunction.space.grid.hierarchicalGrid.interfaceGrid

    traces = {}
    includes = ["dune/mmesh/misc/pyskeletontrace.hh"]
    includes += bulkFunction._includes
    generator = SimpleGenerator(["TraceGF","TraceGF"], "Dune::Fem", pythonname=["TraceGFP","TraceGFM"])
    typeName = []
    for side in [ "in", "out" ]:
        sideStr = "Dune::Fem::IntersectionSide::" + side
        typeName += ["Dune::Fem::TraceGF< " + igrid._typeName + ", " + bulkFunction._typeName + ", " + sideStr + " >"]
    moduleName = "skeleton_"+hashlib.md5(''.join(typeName).encode('utf8')).hexdigest()
    module = generator.load(includes, typeName, moduleName)
    traces["in"]  = module.TraceGFP(igrid, bulkFunction)
    traces["out"] = module.TraceGFM(igrid, bulkFunction)

    import ufl
    import dune.ufl
    trace_p = dune.ufl.GridFunction(traces["in"])
    trace_m = dune.ufl.GridFunction(traces["out"])
    if bulkFunction.scalar:
        trace_p = trace_p.toVectorCoefficient()
        trace_m = trace_m.toVectorCoefficient()
    trace = trace_p
    predefined = {}
    predefined[trace('+')]                     = trace_p
    predefined[trace('-')]                     = trace_m
    predefined[ufl.grad(trace)('+')]           = ufl.grad(trace_p)
    predefined[ufl.grad(trace)('-')]           = ufl.grad(trace_m)
    predefined[ufl.grad(ufl.grad(trace))('+')] = ufl.grad(ufl.grad(trace_p))
    predefined[ufl.grad(ufl.grad(trace))('-')] = ufl.grad(ufl.grad(trace_m))
    trace.predefined = predefined
    if bulkFunction.scalar:
        trace = trace[0]

    bulkFunction.trace = trace
    return trace
################################################################################
