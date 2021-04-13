import sys

if len(sys.argv) > 1:
    name = sys.argv[1]
else:
    name = 'grid.msh'

tikz = []
readNodes = False
readElements = False
edges = set()

file = open(name, 'r')
for line in file.readlines():
    line = line.strip()
    if line == '$Nodes':
        readNodes = True
    if line == '$EndNodes':
        readNodes = False
    if line == '$Elements':
        readElements = True
    if line == '$EndElements':
        readElements = False

    if readNodes:
        data = line.split(' ')
        if len(data) != 4:
            continue
        i, x, y, z = data
        tikz += ['\coordinate ('+i+') at ('+x+','+y+');']

    if readElements:
        data = line.split(' ')
        # normal edges
        if len(data) == 8:
            _, _, _, _, _, a, b, c = line.split(' ')
            for edge in [(a, b), (a, c), (b, c)]:
                p, q = sorted(edge)
                if not (p, q) in edges:
                    tikz += ['\draw[very thin] ('+p+') -- ('+q+');']
                    edges.add((p, q));

        # interface edges
        if len(data) == 7:
            _, _, _, o, p, a, b = line.split(' ')
            if int(o) >= 10 or int(p) >= 10:
                tikz += ['\draw[very thick] ('+a+') -- ('+b+');']

out = open(name+'.tikz', 'w')
out.write("\n".join(tikz))
out.close()
