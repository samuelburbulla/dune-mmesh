# A rectangle grid with a centered horizontal interface
def create(file, h, hf):
    import pygmsh
    import meshio
    import numpy
    geom = pygmsh.built_in.Geometry()
    p1 = geom.add_point([0.0, 0.0, 0.0], lcar=h)
    p2 = geom.add_point([1.0, 0.0, 0.0], lcar=h)
    p3 = geom.add_point([0.0, 1.0, 0.0], lcar=hf)
    p4 = geom.add_point([1.0, 1.0, 0.0], lcar=hf)
    p5 = geom.add_point([0.0, 2.0, 0.0], lcar=h)
    p6 = geom.add_point([1.0, 2.0, 0.0], lcar=h)
    l0 = geom.add_line(p1, p2)
    l1 = geom.add_line(p2, p4)
    lf = geom.add_line(p4, p3)
    l3 = geom.add_line(p3, p1)
    ll0 = geom.add_line_loop([l0, l1, lf, l3])
    rect0 = geom.add_surface(ll0)
    geom.add_physical(rect0, label=0)
    l4 = geom.add_line(p4, p6)
    l5 = geom.add_line(p6, p5)
    l6 = geom.add_line(p5, p3)
    ll1 = geom.add_line_loop([l4, l5, l6, -lf])
    rect1 = geom.add_surface(ll1)
    geom.add_physical(rect1, label=1)
    geom.add_physical(lf, label=10)
    geom.add_physical(p1, label=11)
    geom.add_physical(p2, label=12)
    mesh = pygmsh.generate_mesh(geom, verbose=False)
    meshio.gmsh.write(filename=file, mesh=mesh, fmt_version="2.2", binary=False)
    return geom
