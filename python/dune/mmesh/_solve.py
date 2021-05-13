import logging, traceback
logger = logging.getLogger(__name__)

import numpy as np

def iterativeSolve(schemes, targets, callback=None, iter=100, f_tol=1e-8, factor=1.0, verbose=False):
    """Helper function to solve bulk and interface scheme coupled iteratively.

    Args:
        schemes:  pair of schemes
        targets:  pair of discrete functions
        callback: update function that is called every time before solving a scheme
        iter:     maximum number of iterations
        f_tol:    objective tolerance between two iterates in two norm
        verbose:  print residuum for each iteration
    Note:
        The targets also must be used in the coupling forms.
        We use a vector formulation of Aitken's fix point acceleration proposed by Irons and Tuck.
    """
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

        if res < f_tol:
            converged = True
            break

    return {'converged': converged, 'iterations': i}



def monolithicSolve(schemes, targets, callback=None, iter=10, f_tol=1e-8, eps=1e-6, verbose=0):
    """Helper function to solve bulk and interface scheme coupled monolithically.
       A newton method based on scipy.
       The coupling jacobian blocks are evalutaed by finite difference on demand.
       The Schur complement is used for the inversion of the block-wise jacobian.

    Args:
        schemes:  pair of schemes
        targets:  pair of discrete functions that should be solved for AND that are used in the coupling forms
        callback: update function that is called every time before solving a scheme
        iter:     maximum number of iterations
        f_tol:    objective residual two norm
        eps:      step size for finite difference
        verbose:  1: print residuum for each iteration, 2: and for each Schur iteration

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

    # We use the block structure to solve with the Schur complement
    # jac  = [[ A, B ],  = [[ D_u scheme(u,t),  D_t scheme(u,t)  ],
    #         [ C, D ]]     [ D_u ischeme(u,t), D_t ischeme(u,t) ]]

    from dune.fem.operator import linear as linearOperator
    A = linearOperator(scheme)
    D = linearOperator(ischeme)

    def updateJacobians():
        scheme.jacobian(uh, A)
        ischeme.jacobian(th, D)

    def call():
        if callback is not None:
            callback()

    # Evaluate the coupling blocks by finite difference on-demand
    Bz = uh.copy()
    def B(z):
        norm = np.sqrt(z.dot(z))
        if norm == 0:
            return np.zeros(n)
        th.as_numpy[:] += z * eps / norm
        call()
        scheme(uh, Bz)
        th.as_numpy[:] -= z * eps / norm
        Bz.as_numpy[:] -= f.as_numpy
        Bz.as_numpy[:] /= eps
        return Bz.as_numpy * norm

    Cz = th.copy()
    def C(z):
        norm = np.sqrt(z.dot(z))
        if norm == 0:
            return np.zeros(m)
        uh.as_numpy[:] += z * eps / norm
        call()
        ischeme(th, Cz)
        uh.as_numpy[:] -= z * eps / norm
        Cz.as_numpy[:] -= g.as_numpy
        Cz.as_numpy[:] /= eps
        return Cz.as_numpy * norm

    # Evaluate
    f = uh.copy()
    g = th.copy()

    call()
    scheme(uh, f)
    ischeme(th, g)

    i = 0
    def checkResiduum():
        a = f.as_numpy
        b = g.as_numpy
        res = np.sqrt(np.dot(a, a) + np.dot(b, b))

        if verbose > 0:
            print(" i:", i, " | f | =", res)

        if res < f_tol:
            return True

    if checkResiduum():
        return True

    x = uh.copy().as_numpy
    y = th.copy().as_numpy

    from scipy.sparse.linalg import spsolve, LinearOperator
    from scipy.sparse.linalg import lgmres as iterativeSolver

    def iterCallback(xk):
        if verbose == 2:
            t = evalS(xk)-r
            print("  |Sx-r| =", np.sqrt(t.dot(t)))

    for i in range(1, iter):

        updateJacobians()

        # Compute Newton update
        def evalS(y):
            AinvBy = spsolve(A.as_numpy, B(y))
            return D.as_numpy.dot(y) - C(AinvBy)

        S = LinearOperator((m,m), evalS)
        Ainvf = spsolve(A.as_numpy, f.as_numpy)
        r = g.as_numpy - C(Ainvf)
        y, _ = iterativeSolver(S, r, x0=th.as_numpy, tol=f_tol, callback=iterCallback)
        x = spsolve(A.as_numpy, f.as_numpy - B(y))

        uh.as_numpy[:] -= x
        th.as_numpy[:] -= y

        call()
        scheme(uh, f)
        ischeme(th, g)

        if checkResiduum():
            return True

    return False
