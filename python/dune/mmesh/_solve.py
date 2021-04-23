import logging, traceback
logger = logging.getLogger(__name__)

import numpy as np

def iterativeSolve(schemes, targets, callback=None, iter=100, f_tol=1e-8, factor=1.0, verbose=False):
    """Helper function to solve bulk and interface scheme coupled iteratively.

    Args:
        schemes:  pair of schemes
        targets:  pair of target discrete functions
        callback: update function that is called every time before solving a scheme
        iter:     maximum number of iterations
        f_tol:    objective tolerance between two iterates
        factor:   iteration update damping or acceleration (1: new, 0: old)
        verbose:  print residuum for each iteration

    Note:
        The targets also must be used in the coupling forms.

    """
    assert len(schemes) == 2
    assert len(targets) == 2

    (scheme, scheme_gamma) = schemes
    (ph, ph_gamma) = targets

    ph_old = ph.copy()
    ph_gamma_old = ph_gamma.copy()

    converged = False
    for i in range(iter):
        if callback is not None:
            callback()

        ph_old.assign(ph)
        scheme.solve(target=ph)
        ph.as_numpy[:] = factor*ph.as_numpy + (1-factor)*ph_old.as_numpy

        phnp = ph_old.as_numpy[:]
        phnp -= ph.as_numpy
        error = np.sqrt(np.dot(phnp, phnp))

        if callback is not None:
            callback()

        ph_gamma_old.assign(ph_gamma)
        scheme_gamma.solve(target=ph_gamma)
        ph_gamma.as_numpy[:] = factor*ph_gamma.as_numpy + (1-factor)*ph_gamma_old.as_numpy

        ph_gammanp = ph_gamma_old.as_numpy[:]
        ph_gammanp -= ph_gamma.as_numpy
        error_gamma = np.sqrt(np.dot(ph_gammanp, ph_gammanp))

        if verbose:
            print("{:3d}:".format(i), "[", "{:1.2e}".format(error), " {:1.2e}".format(error_gamma), "]", flush=True)

        if max(error, error_gamma) < f_tol:
            converged = True
            break

    return {'converged': converged, 'iterations': i}
