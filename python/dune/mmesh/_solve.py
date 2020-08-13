import logging, traceback
logger = logging.getLogger(__name__)

import numpy as np
from scipy.sparse import linalg, csc_matrix

import dune.ufl
from dune.fem.operator import linear as linearOperator

def monolithicNewton(schemes, targets, iter=100, f_tol=1e-8, verbose=False, **kwargs):
    """solve bulk and interface scheme coupled monolithically using a newton method

    Args:
        schemes: pair of schemes
        targets: pair of discrete functions that should be solved for AND that are used in the coupling forms
        iter:    maximum number of iterations
        f_tol:   absolute tolerance in maximum norm
        rdiff:   step size for numerical differentiation
        verbose: print residuum for each iteration
        kwargs:  additional arguments passed to the newton_krylov call
    """
    assert len(schemes) == 2
    assert len(targets) == 2

    (scheme, ischeme) = schemes
    (ph, ih) = targets

    n = len(ph.as_numpy)
    m = len(ih.as_numpy)

    if m == 0:
        if verbose:
            print("second scheme is empty, forward to scheme.solve()", flush=True)
        scheme.solve(target=ph)
        return

    Ax = ph.copy()
    iAx = ih.copy()
    def F(x):
        ph.as_numpy[:] = x[:n]
        ih.as_numpy[:] = x[n:]
        scheme(ph, Ax)
        ischeme(ih, iAx)
        return np.concatenate((Ax.as_numpy, iAx.as_numpy))

    x = np.concatenate((ph.as_numpy, ih.as_numpy))

    from scipy.sparse import coo_matrix


    f = F(x)
    for i in range(iter):
        row = []
        col = []
        data = []
        nonzero = 0

        for k in range(n+m):
            xk = x.copy()
            xk[k] += 1e-6
            dFk = F(xk) - f
            dFk /= 1e-6

            for l in range(n+m):
                if abs(dFk[l]) > 1e-12:
                    nonzero += 1
                    row += [l]
                    col += [k]
                    data += [dFk[l]]

        jac = coo_matrix((data, (row, col)), shape=(m+n,m+n))
        print("nonzero", nonzero, "of", (n+m)**2, "entries")

        jac = linalg.spsolve(jac.tocsc(), f)
        jacp = jac[:n]
        ijacp = jac[n:]

        x[:n] -= jacp
        x[n:] -= ijacp

        f = F(x)
        if verbose:
            print("i:", i, "  res =", np.sqrt(np.dot(f, f)))

        if np.max(np.abs(f)) < f_tol:
            break

    ph.as_numpy[:] = x[:n]
    ih.as_numpy[:] = x[n:]


def monolithicSolve(schemes, targets, solver=linalg.gmres, preconditioner=linalg.spilu, **kwargs):
    """solve bulk and interface scheme coupled monolithically

    Args:
        schemes: pair of schemes
        targets: pair of discrete functions that should be solved for AND that are used in the coupling forms
        solver: iterative solver to solve coupled system (default: linalg.gmres)
        preconditioner: preconditioner applied to diagonal blocks (default: linalg.spilu)
        kwargs:  additional arguments passed to the solver
    Returns:
        int: return value of solver
    """
    assert len(schemes) == 2
    assert len(targets) == 2

    (scheme, ischeme) = schemes
    (ph, ih) = targets

    rhs = ph.copy()
    irhs = ih.copy()

    zero = ph.copy()
    zero.clear()
    izero = ih.copy()
    izero.clear()

    scheme(zero, rhs)
    ischeme(izero, irhs)
    rhscon = np.concatenate((-rhs.as_numpy, -irhs.as_numpy))

    n = len(ph.as_numpy)
    m = len(ih.as_numpy)

    Ax = ph.copy()
    iAx = ih.copy()

    def matrixEvaluation(x_coeff):
        ph.as_numpy[:] = x_coeff[:n]
        ih.as_numpy[:] = x_coeff[n:]

        scheme(ph, Ax)
        Ax.as_numpy[:] -= rhs.as_numpy

        ischeme(ih, iAx)
        iAx.as_numpy[:] -= irhs.as_numpy

        return np.concatenate((Ax.as_numpy, iAx.as_numpy))

    Mat = linalg.LinearOperator((n+m,n+m), matrixEvaluation)

    # Preconditioner
    precA = preconditioner( csc_matrix(linearOperator(scheme).as_numpy) )
    preciA = preconditioner( csc_matrix(linearOperator(ischeme).as_numpy) )
    a = np.zeros(n)
    ia = np.zeros(m)
    def precEvaluation(x_coeff):
        a = precA.solve( x_coeff[:n] )
        ia = preciA.solve( x_coeff[n:] )
        return np.concatenate((a, ia))

    prec = linalg.LinearOperator((n+m,n+m), precEvaluation)

    r, info = solver(Mat, rhscon, M=prec, **kwargs)

    ph.as_numpy[:] = r[:n]
    ih.as_numpy[:] = r[n:]

    return info
