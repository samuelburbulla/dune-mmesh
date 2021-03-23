# A rectangle grid with a horizontal centered interface

name = "horizontal.msh"
h = 0.02
hf = h

import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.add(name)

p1 = gmsh.model.geo.addPoint(0, 0, 0, h, 1)
p2 = gmsh.model.geo.addPoint(1, 0, 0, h, 2)
p3 = gmsh.model.geo.addPoint(1, 1, 0, h, 3)
p4 = gmsh.model.geo.addPoint(0, 1, 0, h, 4)

l1 = gmsh.model.geo.addLine(p1, p2, 1)
l2 = gmsh.model.geo.addLine(p2, p3, 2)
l3 = gmsh.model.geo.addLine(p3, p4, 3)
l4 = gmsh.model.geo.addLine(p4, p1, 4)

p5 = gmsh.model.geo.addPoint(0.25, 0.5, 0, hf, 5)
p6 = gmsh.model.geo.addPoint(0.75, 0.5, 0, hf, 6)
lf = gmsh.model.geo.addLine(p5, p6, 10)

gmsh.model.geo.addCurveLoop([l1, l2, l3, l4], 1)
gmsh.model.geo.addPlaneSurface([1], 0)

gmsh.model.geo.synchronize()
gmsh.model.mesh.embed(1, [lf], 2, 0)

gmsh.model.mesh.generate(dim=2)
gmsh.write(name)
gmsh.finalize()
