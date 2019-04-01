#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.

from dolfin import *
import mesh_deformation
import sys

# This file implements a fake biharmonic mesh smoother from
# B. T. Helenbrook: Mesh deformation using the biharmonic operator,
# International Journal for Numerical Methods in Engineering}, 56(7),
# pp. 1007-1021, 2003
# The author is not very clear about which system is solved with which
# boundary conditions. Implementing it either boils down to just solving a
# Poisson problem or something incompatible, so there might be a fundamental
# misunderstanding on my part.
#
# Implementing a real biharmonic mesh smoother (or solver) for the first BVP
# in a mixed formulation appears to be impossible in FEniCS at the moment as
# we need to test functions from the space with correct Dirichlet boundary
# values with test functions from H^1 (and not H^1_0) and for the underlying
# problem we need to test H^1_0 functions with H^1 test functions.
# This lack of symmetry cannot be expressed in FEniCS 1.4 as it is not
# possible to enforce BCs on test functions, or trial functions only.


def gen_mesh(n_gridpoints,dim):
  if dim == 2:
    mesh = RectangleMesh(0,0,1,1,n_gridpoints,n_gridpoints,"left");
  else:
    mesh = BoxMesh(0,0,0,1,1,1,n_gridpoints,n_gridpoints,n_gridpoints);

  return mesh

# Define Dirichlet boundary
class DirichletBoundary(SubDomain):
  def inside(self, x, on_boundary):
    return on_boundary

def mesh_problem(mesh, V):
  delta_t = 0.01

  # Define boundaries
  def everywhere(x,on_boundary): return on_boundary
  def midline(x,on_boundary): return (abs(x[1]-0.5) < 1e-10)

  # If we update the mesh for computation, that means that the
  # computational domain is the domain from the last timestep
  update_mesh = False

  # Output files
  outfile = File("results_biharmonic_smoother/mesh.pvd")
  outfile_u = File("results_biharmonic_smoother/mesh_u.pvd")

  # write out initial mesh
  outfile << mesh

  u = TrialFunction(V)
  v = TestFunction(V)

  deform = Function(V)
  coords_new = Function(V)
  coords_old = Function(V)

  u0 = Function(V)
  w0 = Function(V)
  f = Function(V)

  C = VectorFunctionSpace(mesh, "R", 0)
  VN = MixedFunctionSpace( [V, C])

  [wn, wc] = TrialFunctions(VN)
  [vn, vc] = TestFunctions(VN)
  w = Function(VN)

  A_D = assemble(inner(grad(u), grad(v))*dx)
  A_N = assemble((inner(grad(wn), grad(vn)) + dot(wn, vc) + dot(wc, vn) )*dx)
  M = assemble(inner(u,v)*dx)
  n = FacetNormal(mesh)
  # We need the mesh coordinates in the FE function coords_old, so we solve
  # Laplace with corresponding Dirichlet boundary values as this gives the
  # identity
  rhs = assemble(inner(coords_old,v)*dx)
  bc = DirichletBC(V, mesh_deformation.my_deform.bc1, everywhere)
  bc.apply(A_D, rhs)
  solve(A_D, coords_old.vector(), rhs)

  coords_0 = coords_old.copy(deepcopy=True)

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
    bc0 = DirichletBC(V, mesh_deformation.my_deform.bc0, midline)
    bc1 = DirichletBC(V, mesh_deformation.my_deform.bc1, everywhere)
    bcs = [bc1]

    # On the undeformed mesh with the new bounday
    tmp = coords_old.copy(deepcopy=True)
    for bc in bcs:
      bc.apply(tmp.vector())

    # ... assuming the normal flux is 0
    rhs = assemble( inner(grad(tmp), grad(v))*dx)
    # Complete equation
    #rhs = assemble( inner(grad(tmp), grad(v))*dx - dot(grad(tmp)*n,v)*ds )

    # S is calculated
    solve(M, w0.vector(), rhs)
    #print "Mean value of w0 = ", assemble(w0.sub(0)*dx(mesh)), assemble(w0.sub(1)*dx(mesh))
    w0_m = Constant((assemble(w0.sub(0)*dx), assemble(w0.sub(1)*dx)))
    #w0_m = Constant( (0.5, 0.5) )
    w0_mean = interpolate(w0_m, C)
    #DEBUG
    #print "Mean value of w0 = ", assemble(w0.sub(0)*dx), assemble(w0.sub(1)*dx)
    #outfile_u << w0

    # The error in the solution is calculated
    # rhs = assemble(inner(grad(w0), grad(v))*dx - dot(grad(w0)*n,v)*ds )
    rhs = assemble(inner(grad(w0), grad(v))*dx - dot(grad(w0)*n,v)*ds )
    solve(M, f.vector(), rhs)

    # The negative of the error is treated as a source term
    #rhs = assemble( - dot(f,vn)*dx + dot(grad(w0)*n,vn)*ds )
    rhs = assemble( -dot(f,vn)*dx - dot(grad(w0)*n,vn)*ds +dot(w0_mean,vc)*dx)
    solve(A_N, w.vector(), rhs)

    outfile_u << w.sub(1)
    rhs = assemble(dot(w.sub(0),v)*dx)
    for bc in bcs:
        bc.apply(A_D, rhs)

    solve(A_D, coords_new.vector(), rhs)

    # Deformation: deform = coords_new - coords_old
    deform.vector().zero()
    deform.vector().axpy(-1., coords_old.vector())
    deform.vector().axpy(1., coords_new.vector())

    # If we update the mesh for computation, save the coordinates and reassemble
    # the matrix on the deformed mesh
    if(update_mesh):
      mesh.move(deform)
      mesh_new = mesh
      A_D = assemble(inner(grad(u), grad(v))*dx)
      A_N = assemble((inner(grad(wn), grad(vn)) + dot(wn, vc) + dot(wc, vn) )*dx)
      M = assemble(inner(u,v)*dx)
      n = FacetNormal(mesh)
      coords_old = coords_new.copy(deepcopy=True)
    else:
      mesh_new = Mesh(mesh)
      mesh_new.move(deform)

    # Output
    #outfile_u << deform
    outfile << mesh_new


def main():
  # Optimization options for the form compiler
  parameters["form_compiler"]["cpp_optimize"] = True
  parameters["form_compiler"]["optimize"] = True
  parameters.linear_algebra_backend = "uBLAS"

  dim = 2
  mesh = gen_mesh(16, dim)
  V = VectorFunctionSpace(mesh, "Lagrange", 1)

  mesh_problem(mesh, V)
  sys.exit(0)

if __name__ == '__main__':
        main()
