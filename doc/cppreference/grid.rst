****
Grid
****


.. doxygenclass:: Dune::MMesh
  :members: MMesh(HostGrid),
    getHostGrid,
    interfaceGrid,
    isInterface(Vertex),
    isInterface(Intersection),
    isOnInterface,
    asInterfaceEntity(Intersection),
    asIntersection,
    addInterface,
    addToInterface,
    indicator,
    markElements,
    adapt(),
    ensureVertexMovement,
    moveVertices,
    ensureInterfaceMovement
    moveInterface,
    locate,
