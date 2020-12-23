# A cube grid with a horizontal interface
filename = "horizontal.msh"

def create(h, dim=2):
    assert dim == 2 or dim == 3

    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", h)
    gmsh.model.add(filename)

    kernel = gmsh.model.occ

    if dim == 2:
        box = kernel.addRectangle(0, 0, 0, 1, 1)
        interface = kernel.addLine(
            kernel.addPoint(0.25, 0.5, 0),
            kernel.addPoint(0.75, 0.5, 0)
        )

    if dim == 3:
        box = kernel.addBox(0, 0, 0, 1, 1, 1)
        interface = kernel.addDisk(0.5, 0.5, 0.5, 0.25, 0.25)

    kernel.synchronize()
    gmsh.model.mesh.embed(dim-1, [interface], dim, box)

    gmsh.model.mesh.generate(dim=dim)
    gmsh.write(filename)
    gmsh.finalize()
