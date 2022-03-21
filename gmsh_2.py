import gmsh
import sys
import os
import math

gmsh.initialize()

path = os.path.dirname(os.path.abspath(__file__))

gmsh.merge(os.path.join(path, 'jet.stl'))

angle = 40

forceParametrizablePatches = True

includeBoundary = True

curveAngle = 180


gmsh.model.mesh.classifySurfaces(angle * math.pi / 180., includeBoundary,
                                 forceParametrizablePatches,
                                 curveAngle * math.pi / 180.)

gmsh.model.mesh.createGeometry()

s = gmsh.model.getEntities(2)

l = gmsh.model.geo.addSurfaceLoop([s[i][1] for i in range(len(s))])
gmsh.model.geo.addVolume([l])
gmsh.model.geo.synchronize()

v = gmsh.model.getEntities()

r = 2 #  Меняем размер сетки, более 5 бессмыслено

gmsh.model.mesh.setSize(v, r)

gmsh.model.occ.synchronize()

funny = False
f = gmsh.model.mesh.field.add("MathEval")
if funny:
    gmsh.model.mesh.field.setString(f, "F", "2*Sin((x+y)/5) + 3")
else:
    gmsh.model.mesh.field.setString(f, "F", "4")
gmsh.model.mesh.field.setAsBackgroundMesh(f)

gmsh.model.mesh.generate(3)
gmsh.write('jet.msh')

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
