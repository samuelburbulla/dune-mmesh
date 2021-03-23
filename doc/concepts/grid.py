name = "grid.msh"
h = 0.5
hf = 0.2

import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.add(name)

p1 = gmsh.model.geo.addPoint(0, 0, 0, h, 1)
p2 = gmsh.model.geo.addPoint(3, 0, 0, h, 2)
p3 = gmsh.model.geo.addPoint(3, 3, 0, h, 3)
p4 = gmsh.model.geo.addPoint(0, 3, 0, h, 4)

l1 = gmsh.model.geo.addLine(p1, p2, 1)
l2 = gmsh.model.geo.addLine(p2, p3, 2)
l3 = gmsh.model.geo.addLine(p3, p4, 3)
l4 = gmsh.model.geo.addLine(p4, p1, 4)

p5 = gmsh.model.geo.addPoint(1, 1, 0, hf, 5)
p6 = gmsh.model.geo.addPoint(1.5, 1.5, 0, hf, 6)
p7 = gmsh.model.geo.addPoint(2, 2, 0, hf, 7)
p8 = gmsh.model.geo.addPoint(2, 1, 0, hf, 8)
lf1 = gmsh.model.geo.addLine(p5, p6, 10)
lf2 = gmsh.model.geo.addLine(p7, p6, 11)
lf3 = gmsh.model.geo.addLine(p8, p6, 12)

gmsh.model.geo.addCurveLoop([l1, l2, l3, l4], 1)
gmsh.model.geo.addPlaneSurface([1], 0)

gmsh.model.geo.synchronize()
gmsh.model.mesh.embed(1, [lf1, lf2, lf3], 2, 0)

gmsh.model.mesh.generate(dim=2)
gmsh.write(name)
gmsh.finalize()
