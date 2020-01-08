#!/usr/bin/env python
# FEAT3: Finite Element Analysis Toolbox, Version 3
# Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
# FEAT3 is released under the GNU General Public License version 3,
# see the file 'copyright.txt' in the top level directory for details.

# Rumpf smoother using FEniCS. Original implementation by Steffen Basting,
# modified and enhanced by levelset based r-adaptivity by Jordi Paul

import math
import sys
from dolfin import *
import mesh_deformation
from mesh_deformation import gen_mesh, Top, NotTop
import analytic_lvlset

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"]=8
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

# Code for C++ evaluation of Rumpf trafo
rumpf_trafo_code = """

class RumpfTrafo : public Expression
{
public:

  // Create expression with 4 components
  RumpfTrafo() : Expression(4) {}

  // Function for evaluating expression on each cell
  void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
  {
    const uint D = cell.topological_dimension;
    const uint cell_index = cell.index;
    values[0] = (*c11)[cell_index];
    values[1] = (*c12)[cell_index];
    values[2] = (*c21)[cell_index];
    values[3] = (*c22)[cell_index];
  }

  // The data stored in mesh functions
  std::shared_ptr<MeshFunction<double> > c11;
  std::shared_ptr<MeshFunction<double> > c12;
  std::shared_ptr<MeshFunction<double> > c21;
  std::shared_ptr<MeshFunction<double> > c22;

};
"""


# Computes the sum of the determinants of the transformations of the standard
# Rumpf simplex to the physical cells
def compute_sum_det(mesh):
  sum_det = 0.
  for cell in cells(mesh):
    x = cell.get_vertex_coordinates()
    x1 = x[0]
    x2 = x[2]
    x3 = x[4]
    y1 = x[1]
    y2 = x[3]
    y3 = x[5]
    r11 = x2-x1
    r12 = x3-x1
    r21 = y2-y1
    r22 = y3-y1
    a11 = r11
    a12 = (-r11+2.*r12)/sqrt(3.)
    a21 = r21
    a22 = (-r21+2.*r22)/sqrt(3.)
    sum_det += abs(a11*a22-a12*a21)
  return sum_det

# Concentration function for mesh density
def conc_function(dist):
  alpha = 1e-1
  p = 0.5

  return pow(alpha + abs(dist), p);

# Computes the mesh concentration conc for all elements in the mesh according to
# the given levelset function
def compute_conc(conc, mesh, lvlset_function):
  my_dof_map = lvlset_function.function_space().dofmap()

  for cell in cells(mesh):
    lvlset_val = 0.
    for dof in my_dof_map.cell_dofs(cell.index()) :
      lvlset_val += lvlset_function.vector()[dof]

    conc[cell] = conc_function(lvlset_val)

  return

# Computes the local target cell sizes rumpf_lambda according to a uniform distribution
def compute_lambda_uniform(rumpf_lambda, mesh):
  vol = 0.
  for cell in cells(mesh):
    vol += cell.volume()

  ncells = 1./mesh.num_cells()

  for cell in cells(mesh):
    rumpf_lambda.array()[cell.index()] = ncells*vol

  return

# Computes the local target cell sizes rumpf_lambda according to a nonuniform
# distribution defined by the concentration function
def compute_lambda_nonuniform(rumpf_lambda, mesh, lvlset_function):
  conc = MeshFunction("double", mesh, 2)
  compute_conc(conc, mesh, lvlset_function)

  sum_conc = 0.
  for cell in cells(mesh):
    sum_conc += conc[cell]

  sum_conc = 1./sum_conc
  for cell in cells(mesh):
    rumpf_lambda.array()[cell.index()] = sum_conc*conc.array()[cell.index()]

  return

# Computes the local optimal scales h given by the local target cell sizes
# rumpf_lambda
def compute_h(h, dim, mesh, rumpf_lambda):
  sum_det = compute_sum_det(mesh)
  exponent = 1./dim
  for cell in cells(mesh):
    h.array()[cell.index()] = pow(rumpf_lambda.array()[cell.index()]*sum_det,exponent)

  return

# Evaluates the analytic levelset function at time t and computes its FE
# interpolant lvlset_function
def prepare_lvlset(lvlset_function, analytic_lvlset, t):
  analytic_lvlset.my_lvlset_function.value.t = t
  lvlset_function.vector()[:] = (interpolate(analytic_lvlset.my_lvlset_function.value,\
	  lvlset_function.function_space())).vector()
  return

# Prepares the functional for evaluation and minimisation at time instant t
# by computing the local Rumpf transformation matrix c, the local optimal
# scales h, the local target cell sizes rumpf_lambda according to the given
# mesh
def prepare(c, h, rumpf_lambda, mesh, lvlset_function, use_r_adaptivity, use_rumpf_trafo):
  if(use_r_adaptivity==1):
    compute_lambda_nonuniform(rumpf_lambda, mesh, lvlset_function)
    compute_h(h, 2, mesh, rumpf_lambda)
  else:
    compute_lambda_uniform(rumpf_lambda, mesh)
    compute_h(h, 2, mesh, rumpf_lambda)

  if(use_rumpf_trafo):
    for cell in cells(mesh):
      x = cell.get_vertex_coordinates()
      x1 = x[0]
      x2 = x[2]
      x3 = x[4]
      y1 = x[1]
      y2 = x[3]
      y3 = x[5]
      r11 = x2-x1
      r12 = x3-x1
      r21 = y2-y1
      r22 = y3-y1
      c.c11.array()[cell.index()] = r11/h.array()[cell.index()]
      c.c12.array()[cell.index()] = (-r11+2.*r12)/sqrt(3.)/h.array()[cell.index()]
      c.c21.array()[cell.index()] = r21/h.array()[cell.index()]
      c.c22.array()[cell.index()] = (-r21+2.*r22)/sqrt(3.)/h.array()[cell.index()]
  else:
    for cell in cells(mesh):
      c.c11.array()[cell.index()] = 1.0
      c.c12.array()[cell.index()] = 0.0
      c.c21.array()[cell.index()] = 0.0
      c.c22.array()[cell.index()] = 1.0

def assemble_operator(fac_norm, fac_det, fac_reg, p, C, u, v):

  # Factor for the 1/det term so that the identity mapping minimises the local
  # functional
  if(p==1):
    # det + 1 / det
    fac_rec_det = fac_det*(1. + fac_reg*fac_reg + sqrt(1. + fac_reg*fac_reg))
  elif(p==2):
    # det^2 + 1 / det (obsolete)
    #fac_rec_det = fac_det*( sqrt(1.+fac_reg*fac_reg) * (1. + sqrt(1.+fac_reg*fac_reg))*(1. + sqrt(1.+fac_reg*fac_reg)))
    # det^2 + 1 / det^2
    fac_rec_det = fac_det*(2.*sqrt(fac_reg*fac_reg + 1.) + 2.*fac_reg*fac_reg + 2. + sqrt(fac_reg*fac_reg + 1)\
                  *fac_reg*fac_reg)
  else:
    print "Invalid p:", p
    sys.exit(1)

  # Rumpf functional including transformation to the local Rumpf reference
  # cells. Multiplication by C from the right.
  det_rumpf = det(dolfin.grad(u))*abs(det(C))
  F = fac_norm*4.*inner((inner(dolfin.grad(u)*(C),dolfin.grad(u)*(C))-2.)*dolfin.grad(u)*(C),dolfin.grad(v)*(C))*dx
  if(p==1):
    # det + 1/det
    F += inner((fac_det*det_rumpf - fac_rec_det/ \
	((sqrt(det_rumpf*det_rumpf+fac_reg*fac_reg)+det_rumpf)*sqrt(det_rumpf*det_rumpf+fac_reg*fac_reg)))\
	*inv(transpose(dolfin.grad(u)*C)), dolfin.grad(v)*C)*dx
  else :
    # det^2 + 1/det^2
    F += 2.*inner((fac_det*det_rumpf*det_rumpf - fac_rec_det*det_rumpf \
         /( pow(sqrt(det_rumpf*det_rumpf+fac_reg*fac_reg)+det_rumpf,2.0) \
         *(sqrt(det_rumpf*det_rumpf+fac_reg*fac_reg))))*inv(transpose(dolfin.grad(u)*(C))), dolfin.grad(v)*(C))*dx
  return F

# Mesh optimisation code
def mesh_problem(mesh, function_space_family, function_space_parameter):
  use_r_adaptivity = 0
  # Regularisation parameter for the 1/det term in the Rumpf functional
  fac_reg = (1e-8)
  # Factor for the Frobenius norm term in the Rumpf functional
  fac_norm = (1e-2)
  # Factor for the det term in the Rumpf functional
  fac_det = (1e0)
  # The det term is fac_det*det^exponent_det + fac_rec_det/(det + sqrt(det^2 + fac_reg^2))^exponent_det
  exponent_det = 2
  # If we update the mesh for computation, that means that the computational
  # domain is the domain from the last timestep
  update_mesh = False
  # Use the transformation to the Rumpf reference simplex?
  use_rumpf_trafo = True

  filename="hyperelasticity"
  if(use_rumpf_trafo):
    filename+="_trafo"
  else:
    filename+="_notrafo"

  if(update_mesh):
    filename+="_moving"
  else:
    filename+="_fixed"

  filename+="_fac_norm_"+str(fac_norm)
  filename+="_d"+str(exponent_det)

  # Output file for the deformed mesh
  outfile =  File(filename+"/mesh.pvd")
  # Output file for the deformation plotted on the original mesh
  outfile_u = File(filename+"/mesh_u.pvd")

  # Starting time
  t = 0.
  # Timestep size
  delta_t = 0.01
  # End time
  t_end = 0.5

  print "delta_t = ", delta_t, " t_end = ", t_end
  print "use_r_adaptivity = ", use_r_adaptivity, ", update_mesh = ", update_mesh
  print "fac_norm = ", fac_norm, "fac_det = ", fac_det, "fac_reg = ", fac_reg, "exponent_det = ", exponent_det
  print "writing to directory ", filename
  # Create a copy of the original mesh that we can deform
  deformed_mesh = Mesh(mesh)
  # Create a copy of the original mesh for out possible changing computational domain
  mesh_comp = Mesh(mesh)

  # Define boundaries
  boundaries = FacetFunction("uint", mesh_comp)
  top_boundary = Top()
  other_boundaries = NotTop()
  top_boundary.mark(boundaries, 1)
  other_boundaries.mark(boundaries, 2)

  def everywhere(x,on_boundary): return on_boundary

  mesh_deformation.my_deform.bc0.t = t
  mesh_deformation.my_deform.bc1.t = t

  # Function space for computations on the original mesh
  V_0 = VectorFunctionSpace(mesh, function_space_family, function_space_parameter)
  # Function space for the deformation that transforms the old deformed mesh
  # into the new deformed mesh
  V = VectorFunctionSpace(mesh_comp, function_space_family, function_space_parameter)

  u = TrialFunction(V)
  v = TestFunction(V)
  coords_0 = Function(V)
  # Deformation from the old deformed mesh to the new deformed mesh
  deform = Function(V)
  # Deformation from the original undeformed mesh to the new deformed mesh
  deform_0 = Function(V_0)

  # Levelset function stuff
  lvlset_space = FunctionSpace(mesh_comp, "Lagrange", 1)
  lvlset_function = Function(lvlset_space,name='lvlset')

  scaling = -1.
  displacement = 0.15
  analytic_lvlset.my_lvlset_function.value.displacement = displacement
  analytic_lvlset.my_lvlset_function.value.scaling = scaling

  prepare_lvlset(lvlset_function, analytic_lvlset, t)

  # Target cell sizes. Get computed on the deformed mesh, but have to be
  # copied to ...
  rumpf_lambda = MeshFunction("double", mesh_comp, 2)
  # ... the local optimal scales, defined on the original mesh, as the
  # computations are carried out there
  h = MeshFunction("double", mesh_comp, 2)

  # Local Rumpf transformation matrix
  c = Expression(cppcode=rumpf_trafo_code)
  c.c11 = MeshFunction("double", mesh_comp, 2)
  c.c12 = MeshFunction("double", mesh_comp, 2)
  c.c21 = MeshFunction("double", mesh_comp, 2)
  c.c22 = MeshFunction("double", mesh_comp, 2)

  C = as_matrix(((c[0], c[1]), (c[2], c[3])))

  # We need the coefficients of the FE function representing the original
  #mesh, so we compute a uniform rumpf_lambda and solve the system once
  compute_lambda_uniform(rumpf_lambda, mesh)
  compute_h(h, 2, mesh, rumpf_lambda)

  if(use_rumpf_trafo):
    for cell in cells(mesh_comp):
      x = cell.get_vertex_coordinates()
      x1 = x[0]
      x2 = x[2]
      x3 = x[4]
      y1 = x[1]
      y2 = x[3]
      y3 = x[5]
      r11 = x2-x1
      r12 = x3-x1
      r21 = y2-y1
      r22 = y3-y1
      c.c11.array()[cell.index()] = r11/h.array()[cell.index()]
      c.c12.array()[cell.index()] = (-r11+2.*r12)/sqrt(3.)/h.array()[cell.index()]
      c.c21.array()[cell.index()] = r21/h.array()[cell.index()]
      c.c22.array()[cell.index()] = (-r21+2.*r22)/sqrt(3.)/h.array()[cell.index()]
      #print cell.index(), c.c11.array()[cell.index()], c.c12.array()[cell.index()], \
      #	      c.c21.array()[cell.index()], c.c22.array()[cell.index()]
  else:
    for cell in cells(mesh_comp):
      c.c11.array()[cell.index()] = 1.0
      c.c12.array()[cell.index()] = 0.0
      c.c21.array()[cell.index()] = 0.0
      c.c22.array()[cell.index()] = 1.0

  F = assemble_operator(fac_norm, fac_det, fac_reg, exponent_det, C, u, v)

  coords = interpolate(mesh_deformation.my_deform.bc0, V)
  coords_old = Function(V, name = "coords_old")

  F = action(F, coords)
  J = derivative(F, coords, u)

  # bc1 is Dirichlet bcs on the boundary
  bc = DirichletBC(V, mesh_deformation.my_deform.bc1, everywhere)

  problem = NonlinearVariationalProblem(F, coords, bc, J)
  solver  = NonlinearVariationalSolver(problem)

  prm = solver.parameters
  info(prm, True)
  prm['nonlinear_solver'] = 'newton'
  #prm['newton_solver_parameters']['line_search'] = 'basic'
  prm['newton_solver']['linear_solver']= 'petsc'
  prm['newton_solver']['maximum_iterations'] = 300
  solver.solve()
  # Save original vertex coordinates
  coords_0.vector()[:] = coords.vector()

  deform_0.vector().zero()

  # Now prepare everything again with (potentially) nonuniform rumpf_lambda
  prepare(c, h, rumpf_lambda, mesh_comp, lvlset_function, use_r_adaptivity, use_rumpf_trafo)

  F = assemble_operator(fac_norm, fac_det, fac_reg, exponent_det, C, u, v)

  outfile << rumpf_lambda
  outfile_u << deform_0

  F = action(F, coords)
  J = derivative(F, coords, u)

  # Time loop
  while (t < t_end):
    t += delta_t
    print "t = ", t

    mesh_deformation.my_deform.bc0.t = t
    mesh_deformation.my_deform.bc1.t = t

    # Save old vertex coordinates
    coords_old.vector()[:] = coords.vector()

    # Select boundary conditions here
    if use_r_adaptivity :
      # t = 0 means Picard iteration at t = 0
      prepare_lvlset(lvlset_function, analytic_lvlset, 0)
      prepare(c, h, rumpf_lambda, mesh, lvlset_function, use_r_adaptivity, use_rumpf_trafo)
      bc0 = DirichletBC(V, mesh_deformation.my_deform.bc1, everywhere)
      bcs = [bc0]
    else :
      # The boundary condition should be this:
      #bc0 = DirichletBC(V, mesh_deformation.my_deform.bc1, boundaries, 1)
      # but (due to a bug in FEniCS 1.4?) we have to hard code the expression in this file
      bc0 = DirichletBC(V, dolfin.Expression(("x[0]", "1 + t*sin(2*pi*x[0])"),t=t,pi=math.pi), boundaries, 1)
      bc1 = DirichletBC(V, mesh_deformation.my_deform.bc1, boundaries, 2)
      bcs = [bc0, bc1]

    # If we update the mesh for computation, save the coordinates and reassemble
    # the matrix on the deformed mesh
    if(update_mesh):
      mesh_comp.move(deform)
      F = assemble_operator(fac_norm, fac_det, fac_reg, exponent_det, C, u, v)
      F = action(F, coords)
      J = derivative(F, coords, u)

    problem = NonlinearVariationalProblem(F, coords, bcs, J)
    solver  = NonlinearVariationalSolver(problem)

    prm = solver.parameters
    prm['nonlinear_solver'] = 'newton'
    #prm['newton_solver']['line_search'] = 'basic'
    prm['newton_solver']['linear_solver']= 'petsc'
    prm['newton_solver']['krylov_solver']['gmres']['restart'] = 30
    prm['newton_solver']['maximum_iterations'] = 25
    # When using r-adaptivity, Newton has to be damped for the Krylov solvers
    # to converge
    if(use_r_adaptivity):
      prm['newton_solver']['relaxation_parameter']= 0.1
      prm['newton_solver']['krylov_solver']['gmres']['restart'] = 100

    solver.solve()

    # Compute deformation function from the original mesh to the new one
    deform_0.vector().zero()
    deform_0.vector().axpy(-1., coords_0.vector())
    deform_0.vector().axpy(1., coords.vector())

    # Compute deformation function from the old mesh to the new one
    deform.vector().zero()
    deform.vector().axpy(-1., coords_old.vector())
    deform.vector().axpy(1., coords.vector())

    # Update deformed_mesh for writing out
    deformed_mesh = mesh_comp
    deformed_mesh.move(deform)

    outfile << rumpf_lambda #deformed_mesh
    #outfile << lvlset_function
    outfile_u << deform_0

  return

def main():
  print list_linear_solver_methods()
  dim = 2
  mesh = gen_mesh(32, dim)

  mesh_problem(mesh, "Lagrange", 1)
  sys.exit(0)

if __name__ == '__main__':
  main()
