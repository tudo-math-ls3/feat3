#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# vim: ft=python

# Laplace smoother using FEniCS. Original implementation by Steffen Basting,
# modified and enhanced for computing on the deformed domain by Jordi Paul
import sys
from dolfin import *
import math
import mesh_deformation
from mesh_deformation import gen_mesh, Top, NotTop

def mesh_problem(mesh, V, Q):
  t = 0.
  delta_t = 0.01
  t_end = 0.5

  # If we update the mesh for computation, that means that the
  # computational domain is the domain from the last timestep
  update_mesh = True

  print "delta_t = ", delta_t, " t_end = ", t_end
  print "update_mesh = ", update_mesh

  # Output files
  if(update_mesh):
    outfile = File("results_dudv_smoother_moving/mesh.pvd")
    outfile_u = File("results_dudv_smoother_moving/mesh_u.pvd")
  else:
    outfile = File("results_dudv_smoother/mesh.pvd")
    outfile_u = File("results_dudv_smoother/mesh_u.pvd")

  # write out initial mesh
  outfile << mesh

  # Define boundaries
  boundaries = FacetFunction("uint", mesh)
  top_boundary = Top()
  other_boundaries = NotTop()
  top_boundary.mark(boundaries, 1)
  other_boundaries.mark(boundaries, 2)

  def everywhere(x,on_boundary): return on_boundary

  # Trial and test functions
  u = TrialFunction(V)
  coords_old = Function(V)
  coords_new = Function(V)
  deform = Function(V)
  v = TestFunction(V)

  # Laplace
  A = assemble(0.5*inner(grad(u)+transpose(grad(u)), grad(v)+transpose(grad(v)))*dx)
  rhs = assemble(inner(coords_old,v)*dx)

  # We need the mesh coordinates in the FE function coords_old, so we solve
  # Laplace with corresponding Dirichlet boundary values as this gives the
  # identity
  bc = DirichletBC(V, mesh_deformation.my_deform.bc1, everywhere)
  bc.apply(A, rhs)
  solve(A, coords_old.vector(), rhs)

  while t < t_end:
    t += delta_t
    print "t = ",t

    # Update the boundary values
    mesh_deformation.my_deform.bc0.t = t
    mesh_deformation.my_deform.bc1.t = t

    # The boundary condition should be this:
    #bc0 = DirichletBC(V, mesh_deformation.my_deform.bc1, boundaries, 1)
    # but (due to a bug in FEniCS 1.4?) we have to hard code the expression in this file
    # (the bug leads to the nonlinear solver throwing NaNs in the RumpfSmoother)
    bc0 = DirichletBC(V, dolfin.Expression(("x[0]", "1 + t*sin(2*pi*x[0])"),t=t,pi=math.pi), boundaries, 1)
    bc1 = DirichletBC(V, mesh_deformation.my_deform.bc1, boundaries, 2)
    bcs = [bc0, bc1]

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
      A = assemble(0.5*inner( grad(u)+transpose(grad(u)), grad(v)+transpose(grad(v))) *dx)
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
