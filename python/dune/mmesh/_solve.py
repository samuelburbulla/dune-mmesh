import logging, traceback
logger = logging.getLogger(__name__)

from time import time
import numpy as np

# Return dof vector
def as_vector(df):
  try:
    return df.as_numpy
  except:
    try:
      return df.as_istl
    except:
      raise


def iterativeSolve(schemes, targets, callback=None, iter=100, tol=1e-8, f_tol=None, verbose=False, accelerate=False):
    """Helper function to solve bulk and interface scheme coupled iteratively.

    Args:
        schemes:  pair of schemes
        targets:  pair of discrete functions that should be solved for AND that are used in the coupling forms
        callback: update function that is called every time before solving a scheme
        iter:     maximum number of iterations
        tol:      objective tolerance between two iterates in two norm
        verbose:  print residuum for each iteration
        accelerate: use a vector formulation of Aitken's fix point acceleration proposed by Irons and Tuck.

    Returns:
        if converged and number of iterations

    Note:
        The targets also must be used in the coupling forms.
    """
    if f_tol is not None:
        print("f_tol is deprecated, use tol instead!")
        if tol != 1e-8:
            tol = f_tol

    assert len(schemes) == 2
    assert len(targets) == 2

    (A, B) = schemes
    (u, v) = targets

    a = u.copy()
    b = v.copy()

    if accelerate:
        Fa = a.copy()
        Fb = b.copy()

        FFa = a.copy()
        FFb = b.copy()

    def residuum(a, b):
        return np.sqrt(np.dot(a, a)+np.dot(b, b))

    if callback is not None:
        callback()

    converged = False
    for i in range(iter):
        a.assign(u)
        b.assign(v)

        # Evaluation: a, b -> Fa, Fb
        A.solve(u)
        B.solve(v)

        if accelerate:
            Fa.assign(u)
            Fb.assign(v)

            # Evaluation: Fa, Fb -> FFa, FFb
            if callback is not None:
                callback()
            A.solve(u)
            B.solve(v)
            FFa.assign(u)
            FFb.assign(v)

            # Irons-Tuck update
            DA  = as_vector(FFa) - as_vector(Fa)
            D2a = DA - as_vector(Fa) + as_vector(a)
            as_vector(u)[:] = as_vector(FFa)
            if np.dot(D2a, D2a) != 0:
                as_vector(u)[:] -= np.dot(DA, D2a) / np.dot(D2a, D2a) * DA

            DB  = as_vector(FFb) - as_vector(Fb)
            D2b  = DB - as_vector(Fb) + as_vector(b)
            as_vector(v)[:] = as_vector(FFb)
            if np.dot(D2b, D2b) != 0:
                as_vector(v)[:] -= np.dot(DB, D2b) / np.dot(D2b, D2b) * DB

        if callback is not None:
            callback()

        res = residuum(as_vector(u) - as_vector(a), as_vector(v) - as_vector(b))

        if verbose:
            print("{:3d}:".format(i), "[", "{:1.2e}".format(res), "]", flush=True)

        if res < tol:
            converged = True
            break

    return {'converged': converged, 'iterations': i}



def monolithicSolve(schemes, targets, callback=None, iter=30, tol=1e10, f_tol=1e-7, eps=1.49012e-8, verbose=0, iterative=False):
    """Helper function to solve bulk and interface scheme coupled monolithically.
       A newton method assembling the underlying jacobian matrix.
       The coupling jacobian blocks are evaluated by finite differences.
       We provide an implementation with a fast C++ backend.

    Args:
        schemes:  pair of schemes
        targets:  pair of discrete functions that should be solved for AND that are used in the coupling forms
        callback: update function that is called every time before solving a scheme
        iter:     maximum number of iterations
        tol:      objective residual of iteration step in two norm
        f_tol:    objective residual of function value in two norm
        eps:      step size for finite difference
        verbose:  1: print residuum for each newton iteration, 2: print details
        iterative: Use the experimental iterative solver backend instead of UMFPack. Remark that the solver arguments of the bulk scheme are taken and preconditioning is not supported yet.

    Returns:
        if converged
    """
    assert len(schemes) == 2
    assert len(targets) == 2

    (scheme, ischeme) = schemes
    (uh, th) = targets

    n = len(as_vector(uh))
    m = len(as_vector(th))

    def call():
        if callback is not None:
            callback()

    from numpy.linalg import norm

    # Evaluate
    f = uh.copy()
    g = th.copy()

    call()
    scheme(uh, f)
    ischeme(th, g)

    # Load C++ jacobian implementation
    import hashlib
    from dune.generator import Constructor
    from dune.generator.generator import SimpleGenerator
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if iterative:
      header = "jacobian_iterative.hh"
    else:
      header = "jacobian.hh"
    typeName = "Dune::Python::MMesh::Jacobian< " + scheme.cppTypeName + ", " \
        + ischeme.cppTypeName + ", " + uh.cppTypeName + ", " + th.cppTypeName + " >"
    includes = scheme.cppIncludes + ischeme.cppIncludes + ["dune/python/mmesh/"+header]
    moduleName = "jacobian_" + hashlib.md5(typeName.encode('utf8')).hexdigest()
    constructor = Constructor(['const '+scheme.cppTypeName+'& scheme','const '+ischeme.cppTypeName+' &ischeme', 'const '+uh.cppTypeName+' &uh', 'const '+th.cppTypeName+' &th', 'const double eps', 'const std::function<void()> &callback'],
                              ['return new ' + typeName + '( scheme, ischeme, uh, th, eps, callback );'],
                              ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()', 'pybind11::keep_alive< 1, 6 >()'])

    generator = SimpleGenerator("Jacobian", "Dune::Python::MMesh")
    module = generator.load(includes, typeName, moduleName, constructor)
    jacobian = module.Jacobian(scheme, ischeme, uh, th, eps, call)
    jacobian.init()

    ux = uh.copy()
    tx = th.copy()

    def norm(u, t):
      return np.sqrt(u.scalarProductDofs(u) + t.scalarProductDofs(t))

    for i in range(1, iter+1):

        jacobian.update(uh, th)


        solveTime = -time()
        jacobian.solve(f, g, ux, tx)
        solveTime += time()

        uh -= ux
        th -= tx
        xres = norm(ux, tx)

        call()
        scheme(uh, f)
        ischeme(th, g)

        fres = norm(f, g)

        if verbose > 0 and rank == 0:
            print(" i:", i, " |Δx| =", "{:1.8e}".format(xres), "",  "|f| =", "{:1.8e}".format(fres))
            if verbose > 1:
                print(f"Solve took {solveTime:.6f} seconds.\n")
                file = open('runtime.txt', 'w')
                file.write(str(solveTime))
                file.close()

        if xres < tol and fres < f_tol:
            return True

    return False
