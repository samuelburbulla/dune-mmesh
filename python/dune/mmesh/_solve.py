import logging, traceback
logger = logging.getLogger(__name__)

import numpy as np

def iterativeSolve(schemes, targets, callback=None, iter=100, tol=1e-8, f_tol=None, factor=1.0, verbose=False):
    """Helper function to solve bulk and interface scheme coupled iteratively.

    Args:
        schemes:  pair of schemes
        targets:  pair of discrete functions
        callback: update function that is called every time before solving a scheme
        iter:     maximum number of iterations
        tol:    objective tolerance between two iterates in two norm
        verbose:  print residuum for each iteration
    Note:
        The targets also must be used in the coupling forms.
        We use a vector formulation of Aitken's fix point acceleration proposed by Irons and Tuck.
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

    Fa = a.copy()
    Fb = b.copy()

    FFa = a.copy()
    FFb = b.copy()

    def residuum(a, b):
        return np.sqrt(np.dot(a, a)+np.dot(b, b))

    converged = False
    for i in range(iter):
        a.assign(u)
        b.assign(v)

        # Evaluation: a, b -> Fa, Fb
        if callback is not None:
            callback()
        A.solve(u)
        B.solve(v)
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

        res = residuum(u.as_numpy - a.as_numpy, v.as_numpy - b.as_numpy)

        if verbose:
            print("{:3d}:".format(i), "[", "{:1.2e}".format(res), "]", flush=True)

        if res < tol:
            converged = True
            break

    return {'converged': converged, 'iterations': i}



def monolithicSolve(schemes, targets, callback=None, iter=30, tol=1e-8, f_tol=1e-5, eps=1e-8, verbose=0, python=False):
    """Helper function to solve bulk and interface scheme coupled monolithically.
       A newton method based on scipy.
       The coupling jacobian blocks are evalutaed by finite difference on demand.

    Args:
        schemes:  pair of schemes
        targets:  pair of discrete functions that should be solved for AND that are used in the coupling forms
        callback: update function that is called every time before solving a scheme
        iter:     maximum number of iterations
        tol:      objective residual of iteration step in infinity norm
        f_tol:    objective residual of function value in infinity norm
        eps:      step size for finite difference
        verbose:  1: print residuum for each newton iteration, 2: for each gmres iteration

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

    if python:
        from dune.fem.operator import linear as linearOperator
        A = linearOperator(scheme)
        D = linearOperator(ischeme)

        def updateJacobians():
            scheme.jacobian(uh, A)
            ischeme.jacobian(th, D)

        # Evaluate the coupling blocks by finite difference on-demand
        Bz = uh.copy()
        def B(z):
            if norm(z) == 0:
                return np.zeros(n)
            zh = z * eps / norm(z)
            th.as_numpy[:] += zh
            call()
            scheme(uh, Bz)
            th.as_numpy[:] -= zh
            Bz.as_numpy[:] -= f.as_numpy
            Bz.as_numpy[:] /= eps
            return Bz.as_numpy * norm(z)

        Cz = th.copy()
        def C(z):
            if norm(z) == 0:
                return np.zeros(m)
            zh = z * eps / norm(z)
            uh.as_numpy[:] += zh
            call()
            ischeme(th, Cz)
            uh.as_numpy[:] -= zh
            Cz.as_numpy[:] -= g.as_numpy
            Cz.as_numpy[:] /= eps
            return Cz.as_numpy * norm(z)

    # Evaluate
    f = uh.copy()
    g = th.copy()

    call()
    scheme(uh, f)
    ischeme(th, g)

    if not python:
        # Load C++ jacobian implementation
        import hashlib
        from dune.generator import Constructor
        from dune.generator.generator import SimpleGenerator

        typeName = "Dune::Python::MMesh::Jacobian< " + scheme._typeName + ", " \
            + ischeme._typeName + ", " + uh._typeName + ", " + th._typeName + " >"
        includes = scheme._includes + ischeme._includes + ["dune/python/mmesh/jacobian.hh"]
        moduleName = "jacobian_" + hashlib.md5(typeName.encode('utf8')).hexdigest()
        constructor = Constructor(['const '+scheme._typeName+'& scheme','const '+ischeme._typeName+' &ischeme', 'const '+uh._typeName+' &uh', 'const '+th._typeName+' &th', 'const std::function<void()> &callback'],
                                  ['return new ' + typeName + '( scheme, ischeme, uh, th, callback );'],
                                  ['pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()', 'pybind11::keep_alive< 1, 5 >()'])

        generator = SimpleGenerator("Jacobian", "Dune::Python::MMesh")
        module = generator.load(includes, typeName, moduleName, constructor)
        jacobian = module.Jacobian(scheme, ischeme, uh, th, call)
        jacobian.init();

        ux = uh.copy()
        tx = th.copy()

    for i in range(1, iter+1):

        if python:
            updateJacobians()

            def J(x):
                return (A.as_numpy.dot(x[:n]) + B(x[n:])).tolist() + (C(x[:n]) + D.as_numpy.dot(x[n:])).tolist()

            data = []
            row = []
            col = []
            v = np.zeros(n+m)
            for k in range(n+m):
                v[k] = 1
                dv = J(v)
                v[k] = 0
                for j in range(n+m):
                    if abs(dv[j]) > 1e-14:
                        data += [dv[j]]
                        row += [j]
                        col += [k]

            from scipy.sparse import csr_matrix
            jac = csr_matrix((data, (row, col)))

            r = np.concatenate((f.as_numpy, g.as_numpy))

            from scipy.sparse.linalg import spsolve
            x = spsolve(jac, r)

            uh.as_numpy[:] -= x[:n]
            th.as_numpy[:] -= x[n:]
            xres = norm(x)

        if not python:
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
