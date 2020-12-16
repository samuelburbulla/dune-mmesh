# A cube grid with a horizontal interface
filename = "spheres.msh"

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
        sphere1 = kernel.addCircle(0.2, 0.5, 0, 0.15)
        sphere2 = kernel.addCircle(0.7, 0.5, 0, 0.25)

    if dim == 3:
        box = kernel.addBox(0, 0, 0, 1, 1, 1)
        sp1 = kernel.addSphere(0.2, 0.5, 0.5, 0.15)
        sp2 = kernel.addSphere(0.7, 0.5, 0.5, 0.25)
        kernel.synchronize()
        kernel.cut([(3,box)], [(3,sp1),(3,sp2)], removeTool=False)

        kernel.synchronize()
        sphere1 = gmsh.model.getEntities(dim=2)[-2][1]
        sphere2 = gmsh.model.getEntities(dim=2)[-1][1]

    kernel.synchronize()
    gmsh.model.mesh.embed(dim-1, [sphere1, sphere2], dim, box)

    gmsh.model.mesh.generate(dim=dim)
    gmsh.write(filename)
    gmsh.finalize()
