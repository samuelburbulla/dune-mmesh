"""The solver functionality.

This module include the methods iterativeSolve and monolithicSolve.
"""

import logging
from time import time
import numpy as np
import hashlib
from dune.generator import Constructor
from dune.generator.generator import SimpleGenerator

try:
  from mpi4py import MPI
except ImportError:
  pass

logger = logging.getLogger(__name__)

# Return dof vector
def as_vector(df):
  try:
    return df.as_numpy
  except AttributeError:
    return df.as_istl

#pylint: disable=redefined-builtin
def iterativeSolve(schemes, targets, callback=None, iter=100, tol=1e-8, f_tol=None, verbose=False, accelerate=False):
  """Helper function to solve bulk and interface scheme coupled iteratively.

  Args:
    schemes:  pair of schemes
    targets:  pair of discrete functions that should be solved for AND that are used in the coupling forms
    callback: update function that is called every time before solving a scheme
    iter:   maximum number of iterations
    tol:    objective tolerance between two iterates in two norm
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

  (scheme, ischeme) = schemes
  (u, v) = targets

  a = u.copy()
  b = v.copy()

  if accelerate:
    fa = a.copy()
    fb = b.copy()

    ffa = a.copy()
    ffb = b.copy()

  def residuum(a, b):
    return np.sqrt(np.sum(np.square(a)) + np.sum(np.square(b)))

  if callback is not None:
    callback()

  converged = False
  for i in range(iter):
    a.assign(u)
    b.assign(v)

    # Evaluation: a, b -> fa, fb
    scheme.solve(u)
    ischeme.solve(v)

    if accelerate:
      fa.assign(u)
      fb.assign(v)

      # Evaluation: fa, fb -> ffa, ffb
      if callback is not None:
        callback()
      scheme.solve(u)
      ischeme.solve(v)
      ffa.assign(u)
      ffb.assign(v)

      # Irons-Tuck update
      da  = as_vector(ffa) - as_vector(fa)
      d2a = da - as_vector(fa) + as_vector(a)
      as_vector(u)[:] = as_vector(ffa)
      if np.dot(d2a, d2a) != 0:
        as_vector(u)[:] -= np.dot(da, d2a) / np.dot(d2a, d2a) * da

      db  = as_vector(ffb) - as_vector(fb)
      d2b  = db - as_vector(fb) + as_vector(b)
      as_vector(v)[:] = as_vector(ffb)
      if np.dot(d2b, d2b) != 0:
        as_vector(v)[:] -= np.dot(db, d2b) / np.dot(d2b, d2b) * db

    if callback is not None:
      callback()

    res = residuum(as_vector(u) - as_vector(a), as_vector(v) - as_vector(b))

    if verbose:
      print(f"{i:3d}: [ {res:1.2e} ]", flush=True)

    if res < tol:
      converged = True
      break

  return {"converged": converged, "iterations": i}



def monolithicSolve(schemes, targets, callback=None, iter=30, tol=1e10, f_tol=1e-7, eps=1.49012e-8, verbose=0, iterative=False):
  """Helper function to solve bulk and interface scheme coupled monolithically.
     A newton method assembling the underlying jacobian matrix.
     The coupling jacobian blocks are evaluated by finite differences.
     We provide an implementation with a fast C++ backend.

  Args:
    schemes:  pair of schemes
    targets:  pair of discrete functions that should be solved for AND that are used in the coupling forms
    callback: update function that is called every time before solving a scheme
    iter:   maximum number of iterations
    tol:    objective residual of iteration step in two norm
    f_tol:  objective residual of function value in two norm
    eps:    step size for finite difference
    verbose:  1: print residuum for each newton iteration, 2: print details
    iterative: Use the experimental iterative solver backend instead of UMFPack. Remark that the solver arguments of the bulk scheme are taken
               and preconditioning is not supported yet.

  Returns:
    if converged
  """
  assert len(schemes) == 2
  assert len(targets) == 2

  (scheme, ischeme) = schemes
  (uh, th) = targets

  def call():
    if callback is not None:
      callback()

  # Evaluate
  f = uh.copy()
  g = th.copy()

  call()
  scheme(uh, f)
  ischeme(th, g)

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  # Load C++ jacobian implementation
  if iterative:
    header = "jacobian_iterative.hh"
  else:
    header = "jacobian.hh"
  typeName = "Dune::Python::MMesh::Jacobian< " + scheme.cppTypeName + ", " \
    + ischeme.cppTypeName + ", " + uh.cppTypeName + ", " + th.cppTypeName + " >"
  includes = scheme.cppIncludes + ischeme.cppIncludes + ["dune/python/mmesh/"+header]
  moduleName = "jacobian_" + hashlib.md5(typeName.encode("utf8")).hexdigest()
  constructor = Constructor(["const "+scheme.cppTypeName+"& scheme","const "+ischeme.cppTypeName+" &ischeme", "const "+uh.cppTypeName+" &uh",
                 "const "+th.cppTypeName+" &th", "const double eps", "const std::function<void()> &callback"],
                ["return new " + typeName + "( scheme, ischeme, uh, th, eps, callback );"],
                ["pybind11::keep_alive< 1, 2 >()", "pybind11::keep_alive< 1, 3 >()", "pybind11::keep_alive< 1, 4 >()", "pybind11::keep_alive< 1, 6 >()"])

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
      print(" i:", i, f" |Δx| = {xres:1.8e}  |f| = {fres:1.8e}")
      if verbose > 1:
        print(f"Solve took {solveTime:.6f} seconds.\n")
        with open("runtime.txt", "w", encoding="utf-8") as file:
          file.write(str(solveTime))

    if xres < tol and fres < f_tol:
      return True

  return False
#pylint: enable=redefined-builtin
