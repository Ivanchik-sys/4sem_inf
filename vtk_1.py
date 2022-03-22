import vtkmodules.all as vtk
import numpy as np
import gmsh
import math
import os


class CalcMesh:

    def __init__(self, nodes_coords, tetrs_points):
        self.nodes = np.array([nodes_coords[0::3], nodes_coords[1::3], nodes_coords[2::3]])
        self.a_1, self.a_2 = min(nodes_coords[0::3]), max(nodes_coords[0::3])
        self.b_1, self.b_2 = min(nodes_coords[1::3]), max(nodes_coords[1::3])
        self.c_1, self.c_2 = min(nodes_coords[2::3]), max(nodes_coords[2::3])
        #print(a_1, a_2, b_1, b_2, c_1, c_2)
        self.temperature = np.power(np.power(self.b_2 - self.nodes[1, :], 2) + np.power(self.c_2 - self.nodes[2, :], 2), 0.5) + 20.
        self.temperature_c = self.temperature
        n = int(len(nodes_coords) / 3)
        self.velocity = np.array([[0 for i in range(n)], nodes_coords[1::3] - self.b_1, [0 for i in range(n)]])

        self.tetrs = np.array([tetrs_points[0::4], tetrs_points[1::4], tetrs_points[2::4], tetrs_points[3::4]])
        self.tetrs -= 1

    def move(self, tau, i):
        self.temperature = self.temperature_c + tau * i * (15 + 0.1 * (50. -self.nodes[1, :]))
        self.nodes += self.velocity * tau * self.temperature * 0.002

    def snapshot(self, snap_number):
        unstructuredGrid = vtk.vtkUnstructuredGrid()

        points = vtk.vtkPoints()

        temp = vtk.vtkDoubleArray()
        temp.SetName("temp")

        vel = vtk.vtkDoubleArray()
        vel.SetNumberOfComponents(3)
        vel.SetName("vel")

        for i in range(0, len(self.nodes[0])):

            points.InsertNextPoint(self.nodes[0, i], self.nodes[1, i], self.nodes[2, i])

            temp.InsertNextValue(self.temperature[i])

            vel.InsertNextTuple((self.velocity[0, i], self.velocity[1, i], self.velocity[2, i]))


        unstructuredGrid.SetPoints(points)

        unstructuredGrid.GetPointData().AddArray(temp)
        unstructuredGrid.GetPointData().AddArray(vel)


        for i in range(0, len(self.tetrs[0])):
            tetr = vtk.vtkTetra()
            for j in range(0, 4):
                tetr.GetPointIds().SetId(j, self.tetrs[j, i])
            unstructuredGrid.InsertNextCell(tetr.GetCellType(), tetr.GetPointIds())

        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputDataObject(unstructuredGrid)
        writer.SetFileName("jet3d-step-" + str(snap_number) + ".vtu")
        writer.Write()


gmsh.initialize()

try:
    path = os.path.dirname(os.path.abspath(__file__))
    gmsh.merge(os.path.join(path, 'jet.stl'))
except:
    print("Could not load STL mesh: bye!")
    gmsh.finalize()
    exit(-1)

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

r = 2    #  Меняем размер сетки, более 5 бессмыслено

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

nodeTags, nodesCoord, parametricCoord = gmsh.model.mesh.getNodes()

GMSH_TETR_CODE = 4
tetrsNodesTags = None
elementTypes, elementTags, elementNodeTags = gmsh.model.mesh.getElements()
for i in range(0, len(elementTypes)):
    if elementTypes[i] != GMSH_TETR_CODE:
        continue
    tetrsNodesTags = elementNodeTags[i]

if tetrsNodesTags is None:
    print("Can not find tetra data. Exiting.")
    gmsh.finalize()
    exit(-2)

print("The model has %d nodes and %d tetrs" % (len(nodeTags), len(tetrsNodesTags) / 4))

for i in range(0, len(nodeTags)):

    assert (i == nodeTags[i] - 1)

assert(len(tetrsNodesTags) % 4 == 0)

mesh = CalcMesh(nodesCoord, tetrsNodesTags)
mesh.snapshot(0)

tau = 0.01

for i in range(1, 100):
    mesh.move(tau, i)
    mesh.snapshot(i)

gmsh.finalize()
