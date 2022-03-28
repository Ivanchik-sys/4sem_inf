from fenics import *
from ufl import nabla_grad
from ufl import nabla_div
import os


path = os.path.dirname(os.path.abspath(__file__))

mesh = Mesh(os.path.join(path, 'jet.xml'))

L = 1; W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma

V = VectorFunctionSpace(mesh, 'P', 1)

tol = 1E-14


def clamped_boundary(x, on_boundary):
    return on_boundary and x[0] < tol


bc = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary)


def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)


def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)


u = TrialFunction(V)
d = u.geometric_dimension()
v = TestFunction(V)
f = Constant((0, 0, -rho*g))
T = Constant((0, 0, 0))
a = inner(sigma(u), epsilon(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds


u = Function(V)
solve(a == L, u, bc)

s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
von_Mises = project(von_Mises, V)

u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)

File('elasticity/displacement.pvd') << u
File('elasticity/von_mises.pvd') << von_Mises
File('elasticity/magnitude.pvd') << u_magnitude
