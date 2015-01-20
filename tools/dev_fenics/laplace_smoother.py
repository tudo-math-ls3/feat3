#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# vim: ft=python

# Laplace smoother using FEniCS. Original implementation by Steffen Basting,
# modified and enhanced for computing on the deformed domain by Jordi Paul
import sys
from dolfin import *
import mesh_deformation

def gen_mesh(n_gridpoints,dim):
  if dim == 2:
    mesh = RectangleMesh(0,0,1,1,n_gridpoints,n_gridpoints,"left");
  else:
    mesh = BoxMesh(0,0,0,1,1,1,n_gridpoints,n_gridpoints,n_gridpoints);

  return mesh

def mesh_problem(mesh, V, Q):
  delta_t = 0.01

  # Define boundaries
  def everywhere(x,on_boundary): return on_boundary
  def midline(x,on_boundary): return (abs(x[1]-0.5) < 1e-10)

  # If we update the mesh for computation, that means that the
  # computational domain is the domain from the last timestep
  update_mesh = True

  u = TrialFunction(V)
  coords_old = Function(V)
  coords_new = Function(V)
  deform = Function(V)
  v = TestFunction(V)

  # Laplace
  A = assemble(inner(grad(u), grad(v))*dx)
  rhs = assemble(inner(coords_old,v)*dx)

  # Output files
  outfile = File("results_laplace_smoother/mesh.pvd")
  outfile_u = File("results_laplace_smoother/mesh_u.pvd")

  # write out initial mesh
  outfile << mesh

  # We need the mesh coordinates in the FE function coords_old, so we solve
  # Laplace with corresponding Dirichlet boundary values as this gives the
  # identity
  bc = DirichletBC(V, mesh_deformation.my_deform.bc1, everywhere)
  bc.apply(A, rhs)
  solve(A, coords_old.vector(), rhs)

  t = 0
  while t < 1:
    t += delta_t
    print "t = ",t

    # Update the boundary values
    # Assuming the mesh movement is linear in time...
    if update_mesh:
      mesh_deformation.my_deform.bc0.t = delta_t
      mesh_deformation.my_deform.bc1.t = delta_t
    else:
      mesh_deformation.my_deform.bc0.t = t
      mesh_deformation.my_deform.bc1.t = t

    # Inner boundary condition
    #bc0 = DirichletBC(V, mesh_deformation.my_deform.bc0, midline)
    bc1 = DirichletBC(V, mesh_deformation.my_deform.bc1, everywhere)
    #bcs = [bc0, bc1]
    bcs = [bc1]

    for bc in bcs:
        bc.apply(A, rhs)

    solve(A, coords_new.vector(), rhs)

    # Deformation: deform = coords_new - coords_old
    # As long as the system matrix is not reassembled, this has no effect on
    # the equations
    deform.vector().zero()
    deform.vector().axpy(-1., coords_old.vector())
    deform.vector().axpy(1., coords_new.vector())

    # If we update the mesh for computation, save the coordinates and reassemble
    # the matrix on the deformed mesh
    if(update_mesh):
      mesh.move(deform)
      mesh_new = mesh
      A = assemble(inner(grad(u), grad(v))*dx)
      coords_old = coords_new.copy(deepcopy=True)
    else:
      mesh_new = Mesh(mesh)
      mesh_new.move(deform)

    # Output
    outfile_u << deform
    outfile << mesh_new

  return

def main():
  # Optimization options for the form compiler
  parameters["form_compiler"]["cpp_optimize"] = True
  parameters["form_compiler"]["optimize"] = True
  parameters.linear_algebra_backend = "uBLAS"

  dim = 2
  mesh = gen_mesh(32, dim)
  V = VectorFunctionSpace(mesh, "Lagrange", 1)
  Q = FunctionSpace(mesh, "Lagrange", 1)

  mesh_problem(mesh, V, Q)
  sys.exit(0)

if __name__ == '__main__':
        main()
