import logging, traceback
logger = logging.getLogger(__name__)

import hashlib

from dune.generator.generator import SimpleGenerator

def skeleton(interfaceFunction, grid=None):
    """Return the skeleton representation of a discrete function on the interface grid.

    Args:
        interfaceFunction: The discrete function on the interface grid.
        grid (Grid, optional): The bulk grid. Necessary, if wrapped.

    Returns:
        Skeleton representation of given interface function.

    Note:
        This function has to be restricted when evaluated on facets, e.g. using avg(skeleton).
    """
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


def trace(bulkFunction, igrid=None):
    """Return the trace representation of a discrete function on the bulk grid.

    Args:
        bulkFunction: The discrete function on the bulk grid.
        igrid (InterfaceGrid, optional): The interface grid.

    Returns:
        Trace representation of a given interface function.
        This function has to be restricted to positive ('+') or negative side ('-').
    """
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
