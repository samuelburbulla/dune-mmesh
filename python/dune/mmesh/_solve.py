import logging, traceback
logger = logging.getLogger(__name__)

import numpy as np

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
            DA  = FFa.as_numpy - Fa.as_numpy
            D2a = DA - Fa.as_numpy + a.as_numpy
            u.as_numpy[:] = FFa.as_numpy
            if np.dot(D2a, D2a) != 0:
                u.as_numpy[:] -= np.dot(DA, D2a) / np.dot(D2a, D2a) * DA

            DB  = FFb.as_numpy - Fb.as_numpy
            D2b  = DB - Fb.as_numpy + b.as_numpy
            v.as_numpy[:] = FFb.as_numpy
            if np.dot(D2b, D2b) != 0:
                v.as_numpy[:] -= np.dot(DB, D2b) / np.dot(D2b, D2b) * DB

        if callback is not None:
            callback()

        res = residuum(u.as_numpy - a.as_numpy, v.as_numpy - b.as_numpy)

        if verbose:
            print("{:3d}:".format(i), "[", "{:1.2e}".format(res), "]", flush=True)

        if res < tol:
            converged = True
            break

    return {'converged': converged, 'iterations': i}



def monolithicSolve(schemes, targets, callback=None, iter=30, tol=1e-8, f_tol=1e-5, eps=1.49012e-8, verbose=0):
    """Helper function to solve bulk and interface scheme coupled monolithically.
       A newton method assembling the underlying jacobian matrix.
       The coupling jacobian blocks are evaluated by finite differences.
       We provide a fast version with a C++ backend using UMFPACK.

    Args:
        schemes:  pair of schemes
        targets:  pair of discrete functions that should be solved for AND that are used in the coupling forms
        callback: update function that is called every time before solving a scheme
        iter:     maximum number of iterations
        tol:      objective residual of iteration step in two norm
        f_tol:    objective residual of function value in two norm
        eps:      step size for finite difference
        verbose:  1: print residuum for each newton iteration, 2: print details

    Returns:
        if converged
    """
    assert len(schemes) == 2
    assert len(targets) == 2

    (scheme, ischeme) = schemes
    (uh, th) = targets

    n = len(uh.as_numpy)
    m = len(th.as_numpy)

    if m == 0:
        if verbose:
            print("second scheme is empty, forward to scheme.solve()", flush=True)
        scheme.solve(target=uh)
        return

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

    typeName = "Dune::Python::MMesh::Jacobian< " + scheme._typeName + ", " \
        + ischeme._typeName + ", " + uh._typeName + ", " + th._typeName + " >"
    includes = scheme._includes + ischeme._includes + ["dune/python/mmesh/jacobian.hh"]
    moduleName = "jacobian_" + hashlib.md5(typeName.encode('utf8')).hexdigest()
    constructor = Constructor(['const '+scheme._typeName+'& scheme','const '+ischeme._typeName+' &ischeme', 'const '+uh._typeName+' &uh', 'const '+th._typeName+' &th', 'const double eps', 'const std::function<void()> &callback'],
                              ['return new ' + typeName + '( scheme, ischeme, uh, th, eps, callback );'],
                              ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()', 'pybind11::keep_alive< 1, 6 >()'])

    generator = SimpleGenerator("Jacobian", "Dune::Python::MMesh")
    module = generator.load(includes, typeName, moduleName, constructor)
    jacobian = module.Jacobian(scheme, ischeme, uh, th, eps, call)
    jacobian.init();

    ux = uh.copy()
    tx = th.copy()

    for i in range(1, iter+1):

        jacobian.update(uh, th)
        jacobian.solve(f, g, ux, tx)

        uh -= ux
        th -= tx
        xres = np.sqrt(norm(ux.as_numpy)**2 + norm(tx.as_numpy)**2)

        call()
        scheme(uh, f)
        ischeme(th, g)

        fres = np.sqrt(norm(f.as_numpy)**2 + norm(g.as_numpy)**2)

        if verbose > 0:
            print(" i:", i, " |Δx| =", "{:1.8e}".format(xres), "",  "|f| =", "{:1.8e}".format(fres))

        if xres < tol and fres < f_tol:
            return True

    return False
