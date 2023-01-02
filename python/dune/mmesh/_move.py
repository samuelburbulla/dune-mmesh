"""The moving mesh functionality.
"""

import logging
logger = logging.getLogger(__name__)

def moveMesh(grid, gridFunctions, getShifts=None, igrid=None, igridFunctions=None):
  """Move the mesh and adapt given discrete functions.

  Args:
    grid: The adaptive bulk grid view.
    gridFunctions: List of bulk discrete functions to be adapted.
    getShifts: Callback function to obtain a list with shift per interface vertex index.
    igrid (optional): The adaptive interface grid view.
    igridFunctions (optional): List of interface discrete functions to be adapted.

  Returns:
    Number of ensure and mark steps.
  """
  #pylint: disable=import-outside-toplevel
  from dune.fem import adapt
  #pylint: enable=import-outside-toplevel

  hgrid = grid.hierarchicalGrid

  ensure = 0
  if getShifts is not None:
    shifts = getShifts()
    assert len(shifts) == igrid.size(igrid.dimension)

    while hgrid.ensureInterfaceMovement(2 * shifts) and ensure < 10:
      ensure += 1
      if len(gridFunctions) > 0:
        adapt(gridFunctions)
      if igridFunctions is not None:
        adapt(igridFunctions)
      shifts = getShifts()

    hgrid.moveInterface(shifts)

  mark = 0
  while hgrid.markElements() and mark < 10:
    mark += 1
    if len(gridFunctions) > 0:
      adapt(gridFunctions)
    if igridFunctions is not None:
      adapt(igridFunctions)

  return {"ensure": ensure, "mark": mark}
