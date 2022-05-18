from fenics import *
from ufl import nabla_grad
from ufl import nabla_div
import os
import numpy as np


path = os.path.dirname(os.path.abspath(__file__))

mesh = Mesh(os.path.join(path, 'Accelerator.xml'))

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
E = 70 * 1e9  # module Unga in SI
Lambda = nu * E / ((1 + nu) * (1 - 2*nu))
mu = E / (2 * (1 + nu))  # Lame coefficients aluminium

V = VectorFunctionSpace(mesh, 'P', 1)

domains = MeshFunction("size_t", mesh, mesh.topology().dim())
boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1)


def ctan(x):
    return 1 / tan(x)


def epsilon(u):
    return 0.5 * (nabla_grad(u) + nabla_grad(u).T)


def sigma(u):
    return Lambda * nabla_div(u) * Identity(d) + 2*mu*epsilon(u)


def clamped_boundary(x, on_boundary):
    return on_boundary and near(x[2], min_z)


bc = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary)

vtkfile_1 = File('catapult/displacement.pvd')
vtkfile_2 = File('catapult/von_mises.pvd')

T = 5400  # seconds
omega_0 = 30.367  # rad / s
eps = omega_0 / T
alpha = 0
beta = 35 * pi / 180
g = 9.8  # in SI
rho = 2700  # kg / m3
u_1, u_2 = Function(V), Function(V)
dt = 54
for ti in range(0, T, dt):

    omega = eps * ti
    gamma = (omega * ti) / 2
    a11 = cos(alpha)*cos(gamma) - cos(beta)*sin(alpha)*sin(gamma)
    a12 = -cos(gamma)*sin(alpha) - cos(alpha)*cos(beta)*sin(gamma)
    a13 = sin(beta)*sin(gamma)
    a21 = cos(beta)*cos(gamma)*sin(alpha) + cos(alpha)*sin(gamma)
    a22 = cos(alpha)*cos(beta)*cos(gamma) - sin(alpha)*sin(gamma)
    a23 = -cos(gamma)*sin(beta)
    a31 = sin(alpha)*sin(beta)
    a32 = cos(alpha)*sin(beta)
    a33 = cos(beta)
    A = np.array([[a11, a12, a13], [a21, a22, a23], [a31, a32, a33]])
    v_x, v_y, v_z = A @ np.array([[1], [0], [0]]), A @ np.array([[0], [1], [0]]), A @ np.array([[0], [0], [1]])
    v_g = np.array([0, 0, -1])
    cosx, cosy, cosz = (v_g @ v_x)[0], (v_g @ v_y)[0], (v_g @ v_z)[0]

    f_g_ex, f_g_ey, f_g_ez = "g * rho * cosx", "g * rho * cosy", "g * rho * cosz"
    f_c_ex, f_c_ey, f_c_ez = "x[0] * omega*omega * rho", "x[1] * omega*omega * rho", "0"
    f_aa_x, f_aa_y, f_aa_z = "eps * x[0] * rho", "eps * x[1] * rho", "0"
    f_ex = f_g_ex + "+" + f_c_ex + "+" + f_aa_x
    f_ey = f_g_ey + "+" + f_c_ey + "+" + f_aa_y
    f_ez = f_g_ez + "+" + f_c_ez + "+" + f_aa_z

    f = Expression((f_ex, f_ey, f_ez), g=g, rho=rho, omega=omega, eps=eps, cosx=cosx, cosy=cosy, cosz=cosz, degree=2)

    u = TrialFunction(V)
    d = u.geometric_dimension()
    v = TestFunction(V)
    T = Expression(("0", "0", "0"), degree=2)
    a = inner(sigma(u), epsilon(v))*dx
    L = dot(f, v)*dx + dot(T, v)*ds

    u = Function(V)
    solve(a == L, u, bc)

    s = sigma(u) - (1. / 3) * tr(sigma(u)) * Identity(d)
    von_Mises = sqrt(1.5 * inner(s, s))
    V1 = FunctionSpace(mesh, 'P', 1)
    von_Mises = project(von_Mises, V1)

    u_1.assign(u)
    u_2.assign(von_Mises)
    vtkfile_1 << (u_1, ti * dt)
    vtkfile_2 << (u_2, ti * dt)

