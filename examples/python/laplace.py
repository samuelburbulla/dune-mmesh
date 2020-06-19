## @example mmesh.py
#  This is an example of how to use MMesh with dune-python

import io
from dune.grid import reader
from dune.mmesh import mmesh

dim = 2

if dim == 3:
    file = "../grids/horizontal3d.msh"
else:
    file = "../grids/horizontal2d.msh"

# MMesh
gridView = mmesh((reader.gmsh, file), dim)

print( gridView.size(0) )

gridView.writeVTK("test-python-mmesh")

# InterfaceGrid
igridView = gridView.hierarchicalGrid.interfaceGrid

print( igridView.size(0) )

igridView.writeVTK("test-python-mmesh-interface")

# try:
if True:
    import ufl
    from dune.fem import parameter
    from dune.fem.function import integrate
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
    print("Solve a problem on the interface grid (no coupling")
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
    from dune.fem.function import cppFunction
    code="""
    #include <functional>
    template <class GV, class IGV, class BulkGF>
    auto trace(const IGV &igv, const BulkGF &bulkgf, bool inside) {
      using IElement = typename IGV::template Codim<0>::Entity;
      using LocalCoordinate = typename IElement::Geometry::LocalCoordinate;
      if (inside)
      {
        std::function<typename BulkGF::RangeType(const IElement&,LocalCoordinate)>
          ret = [&bulkgf, lbulk=localFunction(bulkgf)]
          (const auto& interfaceEntity,const auto& xLocal) mutable -> auto {
          const auto intersection =
                bulkgf.space().gridPart().grid().asIntersection( interfaceEntity );
          lbulk.bind(intersection.inside()); // lambda must be mutable so that the non const function can be called
          return lbulk(intersection.geometryInInside().global(xLocal));
        };
        return ret;
      }
      else
      {
        std::function<typename BulkGF::RangeType(const IElement&,LocalCoordinate)>
          ret = [&bulkgf, lbulk=localFunction(bulkgf)]
          (const auto& interfaceEntity,const auto& xLocal) mutable -> auto {
          typename BulkGF::RangeType ret(0);
          const auto intersection =
                bulkgf.space().gridPart().grid().asIntersection( interfaceEntity );
          if (intersection.neighbor())
          {
            lbulk.bind(intersection.outside()); // lambda must be mutable so that the non const function can be called
            ret = lbulk(intersection.geometryInOutside().global(xLocal));
          }
          return ret;
        };
        return ret;
      }
    }
    """

    from dune.fem.function import uflFunction
    # TODO: this uses knowledge about the surface normal direction
    if dim == 3:
        surfaceGrad = lambda w: ufl.as_vector([ufl.grad(w)[1],ufl.grad(w)[2]])
    else:
        surfaceGrad = lambda w: ufl.as_vector([ufl.grad(w)[0]])
    graduh = uflFunction(gridView, name="gradUh", order=space.order-1,
                         ufl=ufl.grad(uh))
    gradInside  = cppFunction(igridView, name="gradInside", order=ispace.order-1,
                              fctName="trace",includes=io.StringIO(code),
                              args=[igridView, graduh,True])
    gradOutside = cppFunction(igridView, name="gradOutside", order=ispace.order-1,
                             fctName="trace",includes=io.StringIO(code),
                             args=[igridView, graduh,False])
    trace = cppFunction(igridView, name="trace", order=ispace.order,
                        fctName="trace",includes=io.StringIO(code),
                        args=[igridView, uh,True])
    # replace 'iexact' with the traces of 'uh' in the bilinear form and bc
    ib = ufl.inner(gradInside+gradOutside, ufl.grad(iv))/2 * ufl.dx
    ischeme = galerkin([ia==ib,DirichletBC(ispace,trace)])
    iuh.interpolate(0)
    ischeme.solve(target=iuh)
    print("  error", integrate(igridView, ufl.dot(iuh-iexact,iuh-iexact),order=5))
    igridView.writeVTK("test-python-mmesh-bulk2interface",
    pointdata=[iuh, gradInside, trace ])


    # showcase coupling - surface to bulk (TODO assumes scalar interface gf)
    ###############################
    print("couple surface to bulk")
    ###############################
    skeletonGF="""
    #include <dune/python/common/typeregistry.hh>
    #include <dune/mmesh/misc/pyskeletonfunction.hh>

    template <class BulkGV, class InterfaceGridFunction>
    auto skeletonGF(const BulkGV &bulkGV, const InterfaceGridFunction &igf) {
      pybind11::object pygf = pybind11::cast( &igf );
      auto cls = Dune::Python::insertClass<SkeletonGF<BulkGV,InterfaceGridFunction>>(pygf,"SkeletonFunction",
                 Dune::Python::GenerateTypeName("SkeletonGF",
                      Dune::MetaType<BulkGV>(),Dune::MetaType<InterfaceGridFunction>()),
                 Dune::Python::IncludeFiles({"dune/mmesh/misc/pyskeletonfunction.hh"})).first;
      Dune::FemPy::registerGridFunction( pygf, cls );
      return SkeletonGF(bulkGV, igf );
    }
    """
    from dune.generator.algorithm import run
    from dune.generator import path
    import dune.ufl
    print("loading",flush=True)
    skeleton = run("skeletonGF", io.StringIO(skeletonGF),gridView,iuh)
    skeleton = dune.ufl.GridFunction(skeleton, scalar=True)
    print("done",flush=True)

    a = ufl.inner(ufl.grad(u),ufl.grad(v)) * ufl.dx
    b = ufl.avg(skeleton)*ufl.avg(v) * ufl.dS
    uh = space.interpolate(0,name="solution")
    scheme = galerkin([a==b,DirichletBC(space,0)], solver="cg",
                      parameters=solverParameters)
    scheme.solve(target=uh)
    gridView.writeVTK("test-python-mmesh-interface2bulk",
        pointdata={"skeleton":uh})

    # Goal 1:
    # dI = gridView.hierarchicalGrid.interfaceMeasure
    # use this for skeleton/boundary intergrals - needs a bit of work for the codegen
    # Goal 2:
    # constrain solution on interface, i.e., extend DirichletConstraints to
    # include general interfaces (see dirichletconstraints feature branch in dune-fem?)
# except ImportError:
#     pass
