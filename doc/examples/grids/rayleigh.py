# A rectangle grid with a horizontal centered interface
name = "rayleigh.msh"
h = 0.05
hf = 0.025

import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Verbosity", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

gmsh.model.add(name)

p1 = gmsh.model.geo.addPoint(0, 0, 0,  h, 1)
p2 = gmsh.model.geo.addPoint(1, 0, 0,  h, 2)
p3 = gmsh.model.geo.addPoint(1, 1, 0, hf, 3)
p4 = gmsh.model.geo.addPoint(1, 2, 0,  h, 4)
p5 = gmsh.model.geo.addPoint(0, 2, 0,  h, 5)
p6 = gmsh.model.geo.addPoint(0, 1, 0, hf, 6)

l1 = gmsh.model.geo.addLine(p1, p2, 1)
l2 = gmsh.model.geo.addLine(p2, p3, 2)
l3 = gmsh.model.geo.addLine(p3, p4, 3)
l4 = gmsh.model.geo.addLine(p4, p5, 4)
l5 = gmsh.model.geo.addLine(p5, p6, 5)
l6 = gmsh.model.geo.addLine(p6, p1, 6)

pf1 = gmsh.model.geo.addPoint(0.75,    1, 0, hf, 7)
pf2 = gmsh.model.geo.addPoint( 0.5, 0.95, 0, hf, 8)
pf3 = gmsh.model.geo.addPoint(0.25,    1, 0, hf, 9)
lf = gmsh.model.geo.addSpline([p3, pf1, pf2, pf3, p6], 10)

gmsh.model.geo.addCurveLoop([l1, l2, lf, l6], 1)
gmsh.model.geo.addPlaneSurface([1], 0)
gmsh.model.geo.addCurveLoop([l3, l4, l5, -lf], 2)
gmsh.model.geo.addPlaneSurface([2], 1)

gmsh.model.geo.synchronize()

gmsh.model.mesh.generate(dim=2)
gmsh.write(name)
gmsh.finalize()
