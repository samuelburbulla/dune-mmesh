import logging, traceback
logger = logging.getLogger(__name__)

from dune.fem import adapt

def moveMesh(grid, gridFunctions, getShifts, igrid=None, igridFunctions=[]):
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
    igridFunctions += [igrid.hierarchicalGrid.one]
    hgrid = grid.hierarchicalGrid
    shifts = getShifts()

    ensure = 0
    while hgrid.ensureInterfaceMovement(2 * shifts) and ensure < 10:
        ensure += 1
        adapt(gridFunctions)
        adapt(igridFunctions)
        shifts = getShifts()

    hgrid.moveInterface(shifts)

    mark = 0
    while hgrid.markElements() and mark < 10:
        mark += 1
        adapt(gridFunctions)
        adapt(igridFunctions)

    return {"ensure": ensure, "mark": mark}
