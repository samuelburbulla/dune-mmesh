import logging, traceback
logger = logging.getLogger(__name__)

import numpy as np
from scipy.sparse import linalg, lil_matrix

import dune.ufl
from dune.fem.operator import galerkin as galerkinOperator

def monolithicSolve(schemes, targets, **kwargs):
    """solve bulk and interface scheme coupled monolithically

    Args:
        schemes: pair of schemes
        targets: pair of discrete functions that should be solved for and that are used in the coupling forms
        kwargs:  additional arguments passed to the solver
        # TODO add possibility to choose solver
    Returns:
        int: return value of solver
    """
    assert len(schemes) == 2
    assert len(targets) == 2

    (scheme, ischeme) = schemes
    (ph, ih) = targets

    A = scheme
    iA = ischeme

    rhs = ph.copy(); rhs.clear()
    irhs = ih.copy(); irhs.clear()

    zero = ph.copy()
    zero.clear()
    izero = ih.copy()
    izero.clear()

    n = len(ph.as_numpy)
    m = len(ih.as_numpy)

    Ax = ph.space.interpolate(0, name="Ax")
    iAx = ih.space.interpolate(0, name="iAx")

    scheme(zero, rhs)
    ischeme(izero, irhs)
    rhscon = np.concatenate((-rhs.as_numpy, -irhs.as_numpy))

    def M(x_coeff):
        ph.as_numpy[:] = x_coeff[:n]
        ih.as_numpy[:] = x_coeff[n:]
        A(ph, Ax)
        Ax.as_numpy[:] -= rhs.as_numpy
        iA(ih, iAx)
        iAx.as_numpy[:] -= irhs.as_numpy
        return np.concatenate((Ax.as_numpy, iAx.as_numpy))

    Mat = lil_matrix((n+m, n+m))
    for i in range(n+m):
        xv = np.zeros(n+m)
        xv[i] = 1
        Mat[i] = M(xv)
    Mat.transpose()
    Matc = Mat.tocsc()
    Mat = Mat.tocsr()

    ILU = linalg.spilu(Matc)
    M_x = lambda x: ILU.solve(x)
    M = linalg.LinearOperator((n+m,n+m), M_x)

    r = linalg.gmres(Mat, rhscon, M=M, **kwargs)
    print(r[1])

    if r[1] == 0:
      r = r[0]

      ph.as_numpy[:] = r[:n]
      ih.as_numpy[:] = r[n:]

    return r[1]
