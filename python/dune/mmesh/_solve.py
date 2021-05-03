import logging, traceback
logger = logging.getLogger(__name__)

import numpy as np

def iterativeSolve(schemes, targets, callback=None, iter=100, f_tol=1e-8, factor=1.0, verbose=False):
    """Helper function to solve bulk and interface scheme coupled iteratively.

    Args:
        schemes:  pair of schemes
        targets:  pair of discrete functions
        iter:     maximum number of iterations
        f_tol:    objective tolerance between two iterates in two norm
        callback: update function that is called every time before solving a scheme
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
