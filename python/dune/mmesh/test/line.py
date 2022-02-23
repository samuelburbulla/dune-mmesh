filename = "line.msh"
h = 0.1
hf = 0.01

import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.add(filename)

geo = gmsh.model.geo
p1 = geo.addPoint(0,  0, 0, h)
p2 = geo.addPoint(1,  0, 0, h)
p3 = geo.addPoint(1, .5, 0, hf)
p4 = geo.addPoint(1,  1, 0, h)
p5 = geo.addPoint(0,  1, 0, h)
p6 = geo.addPoint(0, .5, 0, hf)

l1 = geo.addLine(p1, p2)
l2 = geo.addLine(p2, p3)
l3 = geo.addLine(p3, p4)
l4 = geo.addLine(p4, p5)
l5 = geo.addLine(p5, p6)
l6 = geo.addLine(p6, p1)

lf = geo.addLine(p3, p6)

cl = geo.addCurveLoop([l1, l2, l3, l4, l5, l6])
s1 = geo.addPlaneSurface([cl])

geo.synchronize()
gmsh.model.mesh.embed(1, [lf], 2, s1)

gmsh.model.mesh.generate(dim=2)
gmsh.write(filename)
gmsh.finalize()
