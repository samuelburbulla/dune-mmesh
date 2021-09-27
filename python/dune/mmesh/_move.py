import logging, traceback
logger = logging.getLogger(__name__)

from ufl import *
from ._utility import edgeMovement, interfaceEdgeMovement
from dune.fem.scheme import galerkin
from dune.fem import parameter, adapt

def moveMesh(grid, gridFunctions, getShifts, igrid=None, igridFunctions=None):
    """Move the mesh while keeping a discrete function mass conservative.

    Args:
        grid: The adaptive bulk grid view.
        igrid (optional): The adaptive interface grid view.
        gridFunctions: List [uh, uhOld, ...] of bulk discrete functions to be adapted. uh is the `moved` function while uhOld is used as storage for the old time step.
        igridFunctions: List [iuh, iuhOld, ...] of interface discrete functions to be adapted. iuh is the `moved` function while iuhOld is used as storage for the old time step.
        getShift: Callback function to obtain a list with shift per interface vertex index.

    Returns:
        Number of ensure and mark steps.
    """
    assert len(gridFunctions) >= 2
    uh = gridFunctions[0]
    uhOld = gridFunctions[1]
    uhOld.assign(uh)

    if igrid is not None:
        interface = True
        assert igridFunctions is not None
        assert len(igridFunctions) >= 2

        iuh = igridFunctions[0]
        iuhOld = igridFunctions[1]
        iuhOld.assign(iuh)

    hgrid = grid.hierarchicalGrid

    space = uh.space
    x = SpatialCoordinate(space)
    n = FacetNormal(space)
    u = TrialFunction(space)
    phi = TestFunction(space)
    shifts = getShifts()
    em = edgeMovement(grid, shifts)

    if interface:
        ispace = iuh.space
        ix = SpatialCoordinate(ispace)
        ni = FacetNormal(ispace)
        iu = TrialFunction(ispace)
        iphi = TestFunction(ispace)
        iem = interfaceEdgeMovement(igrid, shifts)

    def h(u, m, n):
        sgn = inner(m('+'), n('+'))
        return conditional( sgn > 0, sgn * u('+'), sgn * u('-') )

    def adv(u, phi, em, n):
        return inner( outer(u, em), grad(phi) ) * dx - inner( h(u, em, n), jump(phi) ) * dS

    scheme = galerkin([
        inner(u * abs(det(grad(x + em))), phi) * dx
        - inner(uhOld, phi) * dx
        + 0.5*(
              adv(uhOld, phi, em, n)
            + adv(u, phi, em, n)
        ) == 0
    ])

    if interface:
        ischeme = galerkin([
            inner(iu * abs(det(grad(ix + iem))), iphi) * dx
            - inner(iuhOld, iphi) * dx
            + 0.5*(
                  adv(iuhOld, iphi, iem, ni)
                + adv(iu, iphi, iem, ni)
            ) == 0
        ])

    parameter["fem.adaptation.method"] = "callback"

    ensure = 0
    while hgrid.ensureInterfaceMovement(shifts) and ensure < 10:
        ensure += 1
        adapt(gridFunctions)
        if interface:
            adapt(igridFunctions + [igrid.hierarchicalGrid.one])
        shifts = getShifts()
        em = edgeMovement(grid, shifts)
        if interface:
            iem = interfaceEdgeMovement(igrid, shifts)

    scheme.solve(uh)
    ischeme.solve(iuh)

    hgrid.moveInterface(shifts)

    mark = 0
    while hgrid.markElements() and mark < 10:
        mark += 1
        adapt(gridFunctions)
        if interface:
            adapt(igridFunctions + [igrid.hierarchicalGrid.one])

    return {"ensure": ensure, "mark": mark}
