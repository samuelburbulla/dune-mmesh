## @example laplace.py
#  This is an example of how to use MMesh with dune-python

import io
from dune.grid import reader
from dune.mmesh import mmesh # , trace, skeletonFunction

import logging
logger = logging.getLogger('dune')
logger.setLevel(logging.INFO)

dim = 2

if dim == 3:
    file = "../grids/horizontal3d.msh"
else:
    file = "../grids/horizontal2d.msh"

# MMesh
gridView = mmesh((reader.gmsh, file), dim)
igridView = gridView.hierarchicalGrid.interfaceGrid

print( igridView.size(0) )

igridView.writeVTK("test-python-mmesh-interface")

import ufl
import dune.ufl
from dune.fem import parameter
from dune.fem.function import integrate, uflFunction
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.ufl import DirichletBC
parameter.append({"fem.verboserank": 0})
solverParameters =\
   {"newton.tolerance": 1e-8,
    "newton.linear.tolerance": 1e-12,
    "newton.linear.preconditioning.method": "jacobi",
    "newton.linear.maxiterations":10000,
    "preconditioning.method": "jacobi",
    "newton.verbose": True,
    "newton.linear.verbose": False}

#####################################################
print("Solve a problem on the bulk grid")
#####################################################
space = lagrange(gridView, order=3)
x = ufl.SpatialCoordinate(space)
u = ufl.TrialFunction(space)
v = ufl.TestFunction(space)
exact = ufl.sin(x[0]*x[1]*4*ufl.pi)
a = ufl.inner(ufl.grad(u),ufl.grad(v)) * ufl.dx
b = -ufl.div(ufl.grad(exact))*v * ufl.dx
uh = space.interpolate(0,name="solution")
scheme = galerkin([a==b,DirichletBC(space,exact)], solver="cg",
                  parameters=solverParameters)
scheme.solve(target=uh)
print("  error", integrate(gridView,ufl.dot(uh-exact,uh-exact),order=5))
gridView.writeVTK("test-python-mmesh-bulk",
    pointdata={"solution":uh,"exact":exact})

###########################################################
print("Solve a problem on the interface grid (no coupling)")
###########################################################
ispace = lagrange(igridView, order=3)
ix = ufl.SpatialCoordinate(ispace)
iu = ufl.TrialFunction(ispace)
iv = ufl.TestFunction(ispace)
if dim==3:
    iexact = ufl.sin(0.5*ix[1]*4*ufl.pi)
else:
    iexact = ufl.sin(ix[0]*0.5*4*ufl.pi)
ia = ufl.inner(ufl.grad(iu),ufl.grad(iv)) * ufl.dx
ib = -ufl.div(ufl.grad(iexact))*iv * ufl.dx
iuh = ispace.interpolate(0,name="solution")
ischeme = galerkin([ia==ib,DirichletBC(ispace,iexact)])
ischeme.solve(target=iuh)
print("  error", integrate(igridView, ufl.dot(iuh-iexact,iuh-iexact),order=5))
igridView.writeVTK("test-python-mmesh-interface", pointdata=[iuh])

###############################
print("couple bulk to surface")
###############################
# hide the following in dune.mmesh
def _trace(bulkGF):
    if False: # try:
        return bulkGF.trace
    else: # except:
        interfaceGV = igridView # bulkGF.space.grid.hierarchicalGrid.interfaceGrid
        traceGF="""
        #include <dune/python/common/typeregistry.hh>
        #include <dune/mmesh/misc/pyskeletonfunction.hh>

        template <Dune::Fem::IntersectionSide side, class InterfaceGV, class BulkGridFunction>
        auto traceGF(const InterfaceGV &iGV, const BulkGridFunction &bgf) {
          using GFType = Dune::Fem::TraceGF<InterfaceGV, BulkGridFunction, side>;
          std::string sideStr = (side==Dune::Fem::IntersectionSide::in)?
                                "Dune::Fem::IntersectionSide::in":"Dune::Fem::IntersectionSide::out";
          pybind11::object pygf = pybind11::cast( &bgf );
          auto cls = Dune::Python::insertClass<GFType>(pygf,"SkeletonPlus",
                     Dune::Python::GenerateTypeName("TraceGF",
                          Dune::MetaType<InterfaceGV>(),Dune::MetaType<BulkGridFunction>(), sideStr),
                     Dune::Python::IncludeFiles({"dune/mmesh/misc/pyskeletonfunction.hh"})).first;
          Dune::FemPy::registerGridFunction( pygf, cls );
          bool scalar = pygf.attr("scalar").template cast<bool>();
          cls.def_property_readonly( "scalar", [scalar] ( GFType &self) { return scalar; } );
          return GFType(iGV, bgf);
        }
        template <class InterfaceGV, class BulkGridFunction>
        auto traceGF(const InterfaceGV &iGV, const BulkGridFunction &bgf) {
          return std::make_pair(
                 traceGF<Dune::Fem::IntersectionSide::in>(iGV, bgf ),
                 traceGF<Dune::Fem::IntersectionSide::out>(iGV, bgf ) );
        }
        """
        from dune.generator.algorithm import run
        from dune.generator import path
        print("loading",flush=True)
        trace   = run("traceGF", io.StringIO(traceGF), interfaceGV, bulkGF)
        # all of this needs hiding in some way - also need to think about
        # how to add this to integrants. Possibly the 'predefined' needs ot
        # be attached to the grid function `trace` in this case and then
        # extracted by the codegen. Then it does not have to be passed in from
        # the outside. So
        #    trace = dune.mmesh(igridView, uh)
        # returns a `ufl Coefficient` that acts like `trace` from above
        # and includes the predefines given below. Then we could write
        #    ib = ufl.inner(grad(avg(trace)), ufl.grad(iv)) * ufl.dx
        #    ischeme = galerkin([ia==ib,DirichletBC(ispace,avg(trace))])
        trace_p = dune.ufl.GridFunction(trace[0])
        trace_m = dune.ufl.GridFunction(trace[1])
        if uh.scalar:
            trace_p = trace_p.toVectorCoefficient()
            trace_m = trace_m.toVectorCoefficient()
        trace   = trace_p # just for now
        predefined = {}
        # predefined[trace]                          = ufl.dx # ufl.FacetNormal(space)
        predefined[trace('+')]                     = trace_p
        predefined[trace('-')]                     = trace_m
        predefined[ufl.grad(trace)('+')]           = ufl.grad(trace_p)
        predefined[ufl.grad(trace)('-')]           = ufl.grad(trace_m)
        predefined[ufl.grad(ufl.grad(trace))('+')] = ufl.grad(ufl.grad(trace_p))
        predefined[ufl.grad(ufl.grad(trace))('-')] = ufl.grad(ufl.grad(trace_m))
        trace.predefined = predefined
        if uh.scalar:
            trace = trace[0]
        uh.trace = trace
        return trace

nBulk = ufl.FacetNormal(space)
nInterface = ufl.FacetNormal(ispace)
assert not nBulk == nInterface

trace = _trace(uh)
# replace 'iexact' with the traces of 'uh' in the bilinear form and bc
ib = ufl.inner(ufl.grad(ufl.avg(trace)), ufl.grad(iv)) * ufl.dx
ischeme = galerkin([ia==ib,DirichletBC(ispace,ufl.avg(trace))])
iuh.interpolate(0)
ischeme.solve(target=iuh)
print("  error", integrate(igridView, ufl.dot(iuh-iexact,iuh-iexact),order=5))
igridView.writeVTK("test-python-mmesh-bulk2interface",
                   pointdata=[iuh, trace ]) # TODO: ufl.grad(iuh) fails for some reason...


###############################
print("couple surface to bulk")
###############################
def _skeleton(bulkGV, interfaceGF):
    skeletonGF="""
    #include <dune/python/common/typeregistry.hh>
    #include <dune/mmesh/misc/pyskeletonfunction.hh>

    template <class BulkGV, class InterfaceGridFunction>
    auto skeletonGF(const BulkGV &bulkGV, const InterfaceGridFunction &igf) {
      using SkeletonGFType = Dune::Fem::SkeletonGF<BulkGV,InterfaceGridFunction>;
      pybind11::object pygf = pybind11::cast( &igf );
      auto cls = Dune::Python::insertClass<SkeletonGFType>(pygf,"SkeletonFunction",
                 Dune::Python::GenerateTypeName("SkeletonGF",
                      Dune::MetaType<BulkGV>(),Dune::MetaType<InterfaceGridFunction>()),
                 Dune::Python::IncludeFiles({"dune/mmesh/misc/pyskeletonfunction.hh"})).first;
      Dune::FemPy::registerGridFunction( pygf, cls );
      bool scalar = pygf.attr("scalar").template cast<bool>();
      cls.def_property_readonly( "scalar", [scalar] ( SkeletonGFType &self) { return scalar; } );
      return SkeletonGFType(bulkGV, igf );
    }
    """
    from dune.generator.algorithm import run
    from dune.generator import path
    print("loading",flush=True)
    skeleton = run("skeletonGF", io.StringIO(skeletonGF), bulkGV, interfaceGF)
    skeleton = dune.ufl.GridFunction(skeleton)
    print("done",flush=True)
    return skeleton

skeleton = _skeleton(gridView, iuh)
a = ufl.inner(ufl.grad(u),ufl.grad(v)) * ufl.dx
b = ufl.avg(skeleton)*ufl.avg(v) * ufl.dS
uh = space.interpolate(0,name="solution")
scheme = galerkin([a==b,DirichletBC(space,0)], solver="cg",
                  parameters=solverParameters)
scheme.solve(target=uh)
gridView.writeVTK("test-python-mmesh-interface2bulk",
    pointdata={"skeleton":uh})
