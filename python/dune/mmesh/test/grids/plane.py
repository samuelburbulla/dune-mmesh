filename = "plane.msh"
h = 0.2
hf = 0.05

import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.add(filename)

occ = gmsh.model.occ
box = occ.addBox(0, 0, 0, 1, 1, 1)
boxPoints = gmsh.model.occ.getEntities(0)

plane = occ.addRectangle(0, 0, 0.5, 1, 1)
planePoints = gmsh.model.occ.getEntities(0)

occ.mesh.setSize(planePoints, hf)
occ.mesh.setSize(boxPoints, h)

occ.fragment([(3, box)], [(2, plane)])
occ.synchronize()

gmsh.model.mesh.generate(dim=3)
gmsh.write(filename)
gmsh.finalize()
