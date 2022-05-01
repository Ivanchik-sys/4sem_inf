from fenics import *
from ufl import nabla_grad
from ufl import nabla_div
import os
from math import *


path = os.path.dirname(os.path.abspath(__file__))

mesh = Mesh(os.path.join(path, 'prototype_1.xml'))

mesh_points = mesh.coordinates()

nodes_z = [mesh_points[i][2] for i in range(len(mesh_points))]
nodes_y = [mesh_points[i][1] for i in range(len(mesh_points))]
nodes_x = [mesh_points[i][0] for i in range(len(mesh_points))]

min_z, max_z = min(nodes_z), max(nodes_z)
min_y, max_y = min(nodes_y), max(nodes_y)
min_x, max_x = min(nodes_x), max(nodes_x)
print('parameters:', min_x, max_x, min_y, max_y, min_z, max_z)
print('lengths:', -min_x + max_x, -min_y + max_y, -min_z + max_z)

nu = 0.34  # poisson coefficient
E = 70 * 1e6  # module Unga
Lambda = nu * E / ((1 + nu) * (1 - 2*nu))
mu = E / (2 * (1 + nu))  # Lame coefficients aluminium

omega = pi / 4
theta = 0
g = 9.8
rho = 2700

V = VectorFunctionSpace(mesh, 'Lagrange', 1)

domains = MeshFunction("size_t", mesh, mesh.topology().dim())
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)


def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)


def sigma(u):
    return Lambda * nabla_div(u) * Identity(d) + 2*mu*epsilon(u)


def boundary(x, on_boundary):
    return on_boundary


f_g_ex, f_g_ey, f_g_ez = "0", "-g * sin(theta)", "-g * cos(theta)"
#f_c_ex, f_c_ey, f_c_ez = 'abs(x[0]) * sqrt(x[0]*x[0] + x[1]*x[1]) * omega*omega', 'abs(x[1]) * sqrt(x[0]*x[0] + x[1]*x[1]) * omega*omega', '0'
f_ex, f_ey, f_ez = f_g_ex, f_g_ey, f_g_ez

f = Expression(('0', '0', '-5'), degree=1)
print('OK')
#f.omega = omega

bc = DirichletBC(V, Constant((0, 0, 0)), boundary)

u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx

u = Function(V)
solve(a == L, u, bc)

File('catapult/displacement.pvd') << u
