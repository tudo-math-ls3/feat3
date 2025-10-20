// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// ====================================================================================================================
// k-epsilon Turbulence Model Simulation (Kuzmin-Mierka) for Backward-Facing Step
// --------------------------------------------------------------------------------------------------------------------
// This application implements a monolithic unsteady incompressible Navier-Stokes solver coupled with the
// k-epsilon turbulence model according to "On the implementation of the k − epsilon turbulence model in
// incompressible ﬂow solvers based on a ﬁnite element discretization" by Kuzmin and Mierka \cite{kuzmin2007implementation}.
//
// The simulation is pre-configured to automatically run the backward-facing step benchmark (backfacestep.xml)
// on refinement level 2. The target Reynolds number is 47625, with a fixed inflow velocity of 1.0.
//
// WARNING: These fixed parameters (level, Re, velocity) are hard-coded for this specific benchmark.
//          Users should only modify them if they know exactly what they are doing.
//
// The primary command line arguments to specify are:
//
// --problem <type>
// Specifies the near-wall boundary condition handling and simulation strategy.
//   - 'dirichlet': Applies Dirichlet boundary conditions and solves directly for Re = 47625.
//   - 'neumann':   Applies Neumann boundary conditions. The simulation starts with a solution for Re = 1000
//                  and then iteratively increases the Reynolds number up to 47625. This gradual
//                  increase is necessary to ensure convergence with this boundary condition type.
//
// --mesh <mesh>
// Specifies the mesh file and has to be ~/feat3.git/data/meshes/backfacestep.xml otherwise the program will not run
//
// To obtain a list of all other command line arguments, which are currently supported by this application, simply
// compile and run this application without any command line arguments and utilize the visual sensors in your skull
// to parse the program output.
//
// The solver must be configured either with SuperLU-MPI for parallel runs or with UMFPACK for serial runs.
// By default, a VTK file named backfacestep_.vtu* is written every 100 time steps.
//
// \author Pia Ritter

#include "base.hpp"
#include "parameter_settings.hpp"
#include "boundary_info.hpp"
#include "stokes_solver.hpp"
#include "time_stepping.hpp"
#include "check_point.hpp"
#include "vtk_writer.hpp"
#include "burgers_assembly_job.hpp"
#include "turb_trace_assembler.hpp"
#include "turb_assembly_full.hpp"

#include <kernel/util/math.hpp>
#include <control/scalar_basic.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/solver/direct_sparse_solver.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/richardson.hpp>

/**
 * \brief Runs the k-epsilon turbulence model simulation for a given mesh configuration.
 *
 * This templated function executes the full simulation of the incompressible Navier-Stokes equations
 * coupled with a k-epsilon turbulence model. It includes a laminar initialization phase and a turbulent phase,
 * where wall functions for k and epsilon are enforced via Dirichlet or Neumann boundary conditions.
 * The function is instantiated with a mesh-specific BoundaryInfoType and is called from main().
 *
 * \tparam BoundaryInfoType Provides boundary-specific data and functions for the selected mesh.
 * \param domain_control          Mesh and domain configuration.
 * \param params                  User-defined simulation parameters.
 * \param comm                    MPI communicator.
 * \param stokes_solver           Navier-Stokes solver infrastructure.
 * \param time_stepping           Time stepping controller for turbulent phase.
 * \param time_stepping_laminar   Time stepping controller for laminar phase.
 * \param check_point             Checkpoint manager for restart capability.
 * \param vtk_writer              VTK output writer for turbulent phase.
 * \param write_laminar           Enables VTK output for laminar phase.
 * \param vtk_writer_laminar      VTK writer used if laminar output is enabled.
 */
template<typename BoundaryInfoType>
void run_k_epsilon_model(
  Turb::DomainControl& domain_control,
  Turb::ParameterSettings& params,
  FEAT::Dist::Comm& comm,
  Turb::StokesSolver& stokes_solver,
  Turb::TimeStepping& time_stepping,
  Turb::TimeStepping& time_stepping_laminar,
  Turb::CheckPoint& check_point,
  Turb::VtkWriter& vtk_writer,
  bool write_laminar,
  Turb::VtkWriter& vtk_writer_laminar
)
{
  using namespace Turb;

  // ================================================================================================
  // Boundary Info Object
  // ------------------------------------------------------------------------------------------------

  typedef Geometry::MeshPart<MeshType> MeshPartType;
  auto* mesh_node = domain_control.front()->get_mesh_node();

  // Create BoundaryInfo Object
  BoundaryInfoType bnd_info(*mesh_node);
  bnd_info.u_mean = params.u_mean;
  bnd_info.C_mu = params.C_mu;
  bnd_info.c_bc = params.c_bc;
  bnd_info.l_0 = params.l_0;

  // Set nu corresponding to Reynolds number and inflow velocity
  DataType nu = params.u_mean * bnd_info.get_l() / params.Re_laminar;
  stokes_solver.nu = nu;

  // Set wall_distance dependent on mesh and level of refinement
  DataType wall_distance = bnd_info.get_wall_dist() / Math::pow(2.0, DataType(domain_control.max_level_index()));

  // streamline diffusion stabilization parameter
  stokes_solver.upsam = 0.3;

  // ================================================================================================
  // Infos to set the inflow boundary condition
  // ------------------------------------------------------------------------------------------------

  // find out whether this process is adjacent to the left inflow boundary
  bool have_bnd_l = (mesh_node->find_mesh_part("bnd:l") != nullptr);

  // ================================================================================================
  // Creating Stokes Solver and Assembling Boundary Conditions
  // ------------------------------------------------------------------------------------------------

  // create levels
  stokes_solver.create_levels();

  // get a reference to the finest stokes level
  StokesLevel& the_stokes_level = *stokes_solver.stokes_levels.front();

  // assemble corresponding boundary conditions depending on which dimension we're in
  stokes_solver.assemble_boundary_conditions(bnd_info.no_flow_parts, bnd_info.slip_parts, bnd_info.flow_parts, bnd_info.flow_func_2d(), false);

  // compile local systems for Vanka smoother
  stokes_solver.compile_local_systems();

  // create multigrid stokes_solver
  stokes_solver.create_multigrid_solver();

  // enable or disable full plot mode
  stokes_solver.plot_nonlinear = stokes_solver.plot_nonlinear_header = time_stepping.full_plot;

  // ================================================================================================
  // Creating Stokes Vectors
  // ------------------------------------------------------------------------------------------------

  // create RHS vector
  GlobalStokesVector vec_stokes_rhs = stokes_solver.create_rhs_vector();

  // create initial solution vector
  GlobalStokesVector vec_stokes_sol  = stokes_solver.create_sol_vector();

  // we need up to four additional stokes vectors for our time stepping scheme
  GlobalStokesVector vec_stokes_sol_1 = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);
  GlobalStokesVector vec_stokes_sol_2 = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);
  GlobalStokesVector vec_stokes_sol_3 = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);
  GlobalStokesVector vec_stokes_tmp = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);

  // register Stokes solution vector with VTK writer
  vtk_writer.register_stokes_vector("stokes", vec_stokes_sol);
  if (write_laminar)
    vtk_writer_laminar.register_stokes_vector("stokes", vec_stokes_sol);

  // register all our Stokes vectors to the checkpoint
  check_point.register_stokes_vector("stokes_sol_1", vec_stokes_sol_1);
  check_point.register_stokes_vector("stokes_sol_2", vec_stokes_sol_2);
  check_point.register_stokes_vector("stokes_sol_3", vec_stokes_sol_3);

  // restart from checkpoint or load solution vector?
  if(check_point.load_checkpoint(time_stepping))
  {
    vec_stokes_sol.copy(vec_stokes_sol_1);
  }
  else if(!stokes_solver.load_joined_sol_vector(vec_stokes_sol))
  {
    // in any other case, we start by solving a Stokes system
    comm.print("\nSolving steady-state Stokes system...\n");

    // solve Stokes system
    if(!stokes_solver.solve_linear(vec_stokes_sol, vec_stokes_rhs))
    {
      comm.print("\nStokes Solver didn't succeed!\n");
      return;
    }
    vec_stokes_sol_1.copy(vec_stokes_sol);
  }

  // ================================================================================================
  // Creating Turb Level
  // ------------------------------------------------------------------------------------------------
  typedef Control::ScalarUnitFilterSystemLevel<DataType, IndexType> ScalarLevelType;
  typedef typename ScalarLevelType::GlobalSystemFilter TurbFilterType;

  ScalarLevelType the_turb_level;
  the_turb_level.assemble_gate(domain_control.front(), domain_control.front()->space_turb);
  the_turb_level.symbolic_assembly_std1(domain_control.front()->space_turb);
  //Control::Asm::asm_splitter(domain_control.front(), [](const auto& dl){return &dl.space_turb;}, the_turb_level.base_splitter_sys);

  // ================================================================================================
  // Vector and Matrix Initialization for TurbAssembler
  // ------------------------------------------------------------------------------------------------

  // Typedefs
  typedef typename ScalarLevelType::GlobalSystemMatrix GlobalTurbMatrix;
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> MatrixTypeTurb;
  typedef LAFEM::DenseVector<DataType, IndexType> VectorTypeTurb;

  // Initialize Matrices
  GlobalTurbMatrix matrix_k = the_turb_level.matrix_sys.clone(LAFEM::CloneMode::Shallow);
  GlobalTurbMatrix matrix_e = the_turb_level.matrix_sys.clone(LAFEM::CloneMode::Layout);
  GlobalTurbMatrix matrix_left_k = the_turb_level.matrix_sys.clone(LAFEM::CloneMode::Layout);
  GlobalTurbMatrix matrix_left_e = the_turb_level.matrix_sys.clone(LAFEM::CloneMode::Layout);
  GlobalTurbMatrix matrix_mass_k = the_turb_level.matrix_sys.clone(LAFEM::CloneMode::Layout);
  GlobalTurbMatrix matrix_mass_e = the_turb_level.matrix_sys.clone(LAFEM::CloneMode::Layout);

  // Initialize Vectors
  GlobalScalarVectorType vec_rhs_k = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_f_k = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_rhs_e = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_f_e = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_sol_k = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_sol_e = matrix_k.create_vector_r();
  GlobalScalarVectorType nu_vec = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_gamma = matrix_k.create_vector_r();
  GlobalStokesVector surface_int_23 = vec_stokes_sol.clone(LAFEM::CloneMode::Deep);
  GlobalScalarVectorType vec_p_k = matrix_k.create_vector_r();
  GlobalStokesVector vec_rhs_u_neumann = stokes_solver.create_rhs_vector();
  GlobalScalarVectorType vec_sol_k_old = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_sol_e_old = matrix_k.create_vector_r();
  GlobalScalarVectorType res_defect_k = matrix_k.create_vector_r();
  GlobalScalarVectorType res_defect_e = matrix_k.create_vector_r();
  GlobalScalarVectorType res_defect_k_primal = matrix_k.create_vector_r();
  GlobalScalarVectorType res_defect_e_primal = matrix_k.create_vector_r();
  GlobalScalarVectorType surface_int_27 = matrix_e.create_vector_r();
  GlobalScalarVectorType vec_p_k_bnd = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_p_k_neumann = matrix_k.create_vector_r();
  GlobalScalarVectorType vec_rhs_e_neumann = matrix_e.create_vector_r();
  GlobalScalarVectorType reference_vec_1 = matrix_k.create_vector_r();
  GlobalScalarVectorType reference_vec_0 = matrix_k.create_vector_r();
  GlobalScalarVectorType left_wall_val_k = matrix_k.create_vector_r();
  GlobalScalarVectorType left_wall_val_e = matrix_k.create_vector_r();

  // Variables for VTK Export
  vtk_writer.register_scalar_vector("k", vec_sol_k);
  vtk_writer.register_scalar_vector("e", vec_sol_e);
  vtk_writer.register_scalar_vector("rhs_k", vec_rhs_k);
  vtk_writer.register_scalar_vector("rhs_e", vec_rhs_e);
  vtk_writer.register_scalar_vector("nu_t", nu_vec);
  vtk_writer.register_scalar_vector("gamma", vec_gamma);

  // specific variable output for Dirichlet/Neumann case
  if(params.problem == "dirichlet")
  {
    vtk_writer.register_stokes_vector("surface_int_23", surface_int_23);
    vtk_writer.register_scalar_vector("p_k", vec_p_k);
  }
  vtk_writer.register_scalar_vector("defect_k_dual", res_defect_k);
  vtk_writer.register_scalar_vector("defect_e_dual", res_defect_e);
  vtk_writer.register_scalar_vector("defect_k_primal", res_defect_k_primal);
  vtk_writer.register_scalar_vector("defect_e_primal", res_defect_e_primal);
  if(params.problem == "neumann")
  {
    vtk_writer.register_scalar_vector("p_k", vec_p_k);
    vtk_writer.register_scalar_vector("p_k_bnd", vec_p_k_bnd);
    vtk_writer.register_scalar_vector("p_k_neumann", vec_p_k_neumann);
    vtk_writer.register_stokes_vector("surface_int_23", surface_int_23);
    vtk_writer.register_scalar_vector("surface_int_27", surface_int_27);
    vtk_writer.register_scalar_vector("ref", reference_vec_0);
  }

  // Format matrices and vectors
  matrix_k.local().format();
  matrix_e.local().format();
  vec_rhs_k.format();
  vec_f_k.format();
  vec_rhs_e.format();
  vec_f_e.format();
  vec_sol_k.format(0.0);
  vec_sol_e.format(0.0);
  nu_vec.format(0.0);
  vec_gamma.format(0.0);
  surface_int_23.format(0.0);
  vec_p_k.format(0.0);
  vec_rhs_u_neumann.format();
  vec_rhs_e_neumann.format();
  vec_sol_k_old.format(0.0);
  vec_sol_e_old.format(0.0);
  res_defect_k.format(0.0);
  res_defect_e.format(0.0);
  res_defect_k_primal.format(0.0);
  res_defect_e_primal.format(0.0);
  left_wall_val_k.format(0.0);
  left_wall_val_e.format(0.0);
  vec_p_k.format(0.0);
  vec_p_k_bnd.format(0.0);
  vec_p_k_neumann.format(0.0);
  surface_int_23.format(0.0);
  surface_int_27.format(0.0);
  reference_vec_1.format(1.0);
  reference_vec_0.format(0.0);

  // Create TurbAssemblyJob
  TurbSystemAssemblyJob<MatrixTypeTurb, VectorTypeTurb, LocalVeloVector, SpaceTurbType, SpaceVeloType> turb_asm_job;
  turb_asm_job.matrix_k = &matrix_k.local();
  turb_asm_job.matrix_e = &matrix_e.local();
  turb_asm_job.vec_rhs_k = &vec_f_k.local();
  turb_asm_job.vec_rhs_e = &vec_f_e.local();
  turb_asm_job.vec_sol_k = &vec_sol_k.local();
  turb_asm_job.vec_sol_e = &vec_sol_e.local();
  turb_asm_job.vec_sol_v = &vec_stokes_sol.local().first();
  if (params.use_previous_velo)
    turb_asm_job.vec_sol_v = &vec_stokes_sol_1.local().first();
  turb_asm_job.space_turb = &domain_control.front()->space_turb;
  turb_asm_job.space_velo = &domain_control.front()->space_velo;
  turb_asm_job.sigma_k = params.sigma_k;
  turb_asm_job.sigma_e = params.sigma_e;
  turb_asm_job.C_mu = params.C_mu;
  turb_asm_job.l_max = bnd_info.get_l_max();
  turb_asm_job.nu_min = params.nu_min;
  turb_asm_job.theta_matrix_e = params.C_2;
  turb_asm_job.theta_matrix_k = DataType(1.0);
  turb_asm_job.theta_vec_rhs = params.C_1;
  turb_asm_job.scal_diff = params.scal_diff;
  turb_asm_job.scal_conv = params.scal_conv;
  turb_asm_job.scal_reac = params.scal_reac;
  turb_asm_job.scal_p_k = params.scal_rhs;
  turb_asm_job.problem = params.problem;
  turb_asm_job.nu = params.u_mean*bnd_info.get_l() / params.Re_turbulent;
  turb_asm_job.ref = &reference_vec_0.local();
  turb_asm_job.p_k = &vec_p_k.local();
  turb_asm_job.p_k_bnd = &vec_p_k_bnd.local();
  turb_asm_job.p_k_neumann = &vec_p_k_neumann.local();

  // start with a smaller Re for Neumann
  if(params.problem == "neumann")
  {
    params.Re_turbulent = 1000;
    nu = params.u_mean*bnd_info.get_l() / params.Re_turbulent;
    turb_asm_job.nu = params.u_mean*bnd_info.get_l() / params.Re_turbulent;
    turb_asm_job.ref = &reference_vec_0.local();
    turb_asm_job.p_k = &vec_p_k.local();
    turb_asm_job.p_k_bnd = &vec_p_k_bnd.local();
    turb_asm_job.p_k_neumann = &vec_p_k_neumann.local();
  }

  // ================================================================================================
  // UnitFilter for k and e for InflowBoundary
  // ------------------------------------------------------------------------------------------------
  // Initialize and create Filter and FilterAssembler

  // filter for eq (18)
  TurbFilterType filter_k_wall;
  TurbFilterType filter_e_wall;

  // filter for eq (27)
  TurbFilterType filter_k_neumann;
  TurbFilterType filter_e_neumann;

  // our unit assembler
  Assembly::UnitFilterAssembler<MeshType> unit_asm;
  Assembly::UnitFilterAssembler<MeshType> unit_asm_wall;

  // Find and add mesh part
  MeshPartType* bnd_left = mesh_node->find_mesh_part("bnd:l");
  if(bnd_left != nullptr)
    unit_asm.add_mesh_part(*bnd_left);

  if(params.problem == "dirichlet")
  {
    for (const auto* part : bnd_info.trace_mesh_parts)
    {
      if(part != nullptr)
      {
        unit_asm_wall.add_mesh_part(*part);
      }
    }
  }

  // Assemble unit filter
  if(params.problem == "neumann")
  {
    unit_asm.assemble(filter_k_neumann.local(), domain_control.front()->space_turb, bnd_info.k_sol_func());
    unit_asm.assemble(filter_e_neumann.local(), domain_control.front()->space_turb, bnd_info.e_sol_func());
  }

  // create wall values
  GlobalScalarVectorType k_wall = matrix_k.create_vector_r();
  GlobalScalarVectorType e_wall = matrix_k.create_vector_r();
  DataType u_tau = 0;

  if(params.problem == "dirichlet")
  {
    // compute u_tau, k_wall and e_wall according to eq (18)
    for(Index i = 0; i < k_wall.local().size(); ++i)
    {
      u_tau = Math::max(Math::pow(params.C_mu, 0.25) * Math::sqrt(vec_sol_k.local()(i)), vec_stokes_sol.local().first()(i).norm_euclid()/11.06);
      k_wall.local()(i, pow(u_tau, 2.0) / Math::sqrt(params.C_mu));
      e_wall.local()(i, pow(u_tau, 3.0) / (0.41 * wall_distance));
    }

    // assemble wall filter
    unit_asm_wall.assemble(filter_k_wall.local(), domain_control.front()->space_turb, k_wall.local());
    unit_asm_wall.assemble(filter_e_wall.local(), domain_control.front()->space_turb, e_wall.local());

    // evaluate k_sol_func and e_sol_func into a vector
    Assembly::Interpolator interpolator;
    interpolator.project(left_wall_val_k.local(), bnd_info.k_sol_func(), domain_control.front()->space_turb);
    interpolator.project(left_wall_val_e.local(), bnd_info.e_sol_func(), domain_control.front()->space_turb);

    // loop over left wall indices and set filter value manually because we have a different value at the inflow
    if(have_bnd_l)
    {
      const FEAT::Geometry::TargetSet& vertices_bnd_l = mesh_node->find_mesh_part("bnd:l")->get_target_set<0>();

      for(Index i = 0; i < vertices_bnd_l.get_num_entities(); ++i)
      {
        filter_k_wall.local().add(vertices_bnd_l[i], left_wall_val_k.local()(vertices_bnd_l[i]));
        filter_e_wall.local().add(vertices_bnd_l[i], left_wall_val_e.local()(vertices_bnd_l[i]));
      }
    }

    // export k_wall and e_wall vector to make sure boundary values are enforced
    vtk_writer.register_scalar_vector("k wall", k_wall);
    vtk_writer.register_scalar_vector("e wall", e_wall);
  }

  // Apply filter
  if(params.problem == "dirichlet")
    filter_k_wall.filter_sol(vec_sol_k);
  else if (params.problem == "neumann")
    filter_k_neumann.filter_sol(vec_sol_k);

  if(params.problem == "dirichlet")
    filter_e_wall.filter_sol(vec_sol_e);
  else if (params.problem == "neumann")
    filter_e_neumann.filter_sol(vec_sol_e);

  // ================================================================================================
  // Adding mesh parts to TraceAssembler for u, k and e
  // ------------------------------------------------------------------------------------------------

  // Creating TurbTraceAssembler
  TurbTraceAssembler<decltype(domain_control.front()->trafo)> turb_trace_assembler(domain_control.front()->trafo, wall_distance);

  if(params.problem == "neumann")
  {
    // Create reference filter and assembler (for eqs (24) and (28))
    TurbFilterType reference_filter;
    Assembly::UnitFilterAssembler<MeshType> reference_asm;

    // Add wall-meshparts to turb_trace_assembler and reference filter
    for (const auto* part : bnd_info.trace_mesh_parts)
    {
      if(part != nullptr)
      {
        turb_trace_assembler.add_mesh_part(*part);
        reference_asm.add_mesh_part(*part);
      }
    }
    // compile TurbTraceAssembler
    turb_trace_assembler.compile();

    // Assemble and filter reference vectors
    reference_asm.assemble(reference_filter.local(), domain_control.front()->space_turb, reference_vec_1.local());
    reference_filter.filter_sol(reference_vec_0);
  }

  if(params.problem == "dirichlet")
  {
    for (const auto* part : bnd_info.trace_mesh_parts)
    {
      if(part != nullptr)
        turb_trace_assembler.add_mesh_part(*part);
    }
    turb_trace_assembler.compile();
  }

  // ================================================================================================
  // Matrix and Vector Assembly
  // ------------------------------------------------------------------------------------------------

  comm.print("Assembling system matrix...\n");

  // Symbolic Assembly matrix_left_* and matrix_mass_*
  Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_left_k.local(), domain_control.front()->space_turb);
  Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_left_e.local(), domain_control.front()->space_turb);
  matrix_left_k.local().format();
  matrix_left_e.local().format();
  matrix_mass_k = matrix_left_k.clone(LAFEM::CloneMode::Weak);
  matrix_mass_e = matrix_left_e.clone(LAFEM::CloneMode::Weak);

  // Assemble mass matrices
  Assembly::Common::IdentityOperator identity_operator;
  Assembly::assemble_bilinear_operator_matrix_1(
    domain_control.front()->domain_asm,
    matrix_mass_k.local(),
    identity_operator,
    domain_control.front()->space_turb,
    params.cubature_name
  );

  Assembly::assemble_bilinear_operator_matrix_1(
    domain_control.front()->domain_asm,
    matrix_mass_e.local(),
    identity_operator,
    domain_control.front()->space_turb,
    params.cubature_name
  );

  // assemble TurbAssemblyJob
  domain_control.front()->domain_asm.assemble(turb_asm_job);

  // Create vec_mass_l for Computation of primal vector P_k
  auto vec_mass_l = matrix_mass_k.lump_rows();
  vec_mass_l.component_invert(vec_mass_l, 1.0);

  // ================================================================================================
  // Initial Solution
  // ------------------------------------------------------------------------------------------------

  // equation (10)
  vec_sol_k.format(params.k_0);
  vec_sol_e.format(params.e_0);

  // equation (23)
  // convert to obtain a type-0 vector
  surface_int_23.from_1_to_0();
  turb_trace_assembler.assemble_functional_vector_velocity(
    surface_int_23.local().first(),
    vec_sol_k.local(),
    vec_sol_e.local(),
    vec_stokes_sol.local().first(),
    domain_control.front()->space_turb,
    domain_control.front()->space_velo,
    params.cubature_name,
    DataType(1.0)
  );

  // synchronize to obtain a type-1 vector
  surface_int_23.sync_0();

  if(params.problem == "neumann")
  {
    // equation (27)
    // convert to obtain a type-0 vector
    surface_int_27.from_1_to_0();
    turb_trace_assembler.assemble_functional_vector_dissipation(
      surface_int_27.local(),
      vec_sol_k.local(),
      vec_sol_e.local(),
      vec_stokes_sol.local().first(),
      domain_control.front()->space_turb,
      domain_control.front()->space_velo,
      params.cubature_name,
      DataType(1.0)
    );

    // synchronize to obtain a type-1 vector
    surface_int_27.sync_0();
  }

  // ================================================================================================
  // Preconditioner and Solver
  // ------------------------------------------------------------------------------------------------

  // Create Preconditioner and Solver
  auto& solver_filter_k = (params.problem == "neumann" ? filter_k_neumann : filter_k_wall);
  auto& solver_filter_e = (params.problem == "neumann" ? filter_e_neumann : filter_e_wall);

  // these solvers require a suitable third-party direct solver (UMFPACK, SuperLU, MKL-DSS or cuDSS)
  auto sparse_sol_k = Solver::new_direct_sparse_solver(matrix_left_k, solver_filter_k);
  auto sparse_sol_e = Solver::new_direct_sparse_solver(matrix_left_e, solver_filter_e);
  auto solver_k = Solver::new_richardson(matrix_left_k, solver_filter_k, 1, sparse_sol_k);
  auto solver_e = Solver::new_richardson(matrix_left_e, solver_filter_e, 1, sparse_sol_e);

  // Enable convergence plot
  solver_k->set_plot_mode(Solver::PlotMode::summary);
  solver_k->set_max_iter(10000);
  solver_k->set_tol_rel(1e-10);

  solver_e->set_plot_mode(Solver::PlotMode::summary);
  solver_e->set_max_iter(10000);
  solver_e->set_tol_rel(1e-10);

  // Initialize the solver
  solver_k->init_symbolic();
  solver_e->init_symbolic();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Laminar Phase
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  if(params.do_laminar)
  {
    // Print message
    comm.print(String("\n") + String(100u, '=') + "\n");
    comm.print("Starting laminar phase...\n\n");
    comm.print(
      "Laminar Reynolds Information:\n\n  nu_laminar          : " + stringify(nu)
      + "\n  Re_laminar          : " + stringify(params.Re_laminar)
      + "\n  u_mean (laminar)    : " + stringify(params.u_mean) + "\n\n"
    );

    // StopWatches for time information
    StopWatch laminar_extrapolate, laminar_filter_sol, laminar_rhs, laminar_filter_rhs, laminar_solve;

    // begin time stepping loop
    time_stepping_laminar.begin_loop();

    while(time_stepping_laminar.advance_step())
    {
      comm.print("Current time step " + stringify(time_stepping_laminar.time_step));

      // ============================================================================================
      // Step 1: extrapolate previous Navier-Stokes solution to obtain an initial guess
      // --------------------------------------------------------------------------------------------

      // start watch
      laminar_extrapolate.start();
      the_stokes_level.extrapolate_sol(
        vec_stokes_sol,
        time_stepping_laminar.expolc,
        vec_stokes_sol_1,
        vec_stokes_sol_2,
        vec_stokes_sol_3
      );

      // stop watch
      laminar_extrapolate.stop();

      // start watch
      laminar_filter_sol.start();

      // apply filter
      the_stokes_level.filter_sys.filter_sol(vec_stokes_sol);

      // stop watch
      laminar_filter_sol.stop();

      // ============================================================================================
      // Step 2: assemble RHS vector for Navier-Stokes
      // --------------------------------------------------------------------------------------------

      // start watch
      laminar_rhs.start();

      // format RHS
      vec_stokes_rhs.format();

      // add all terms from the time stepping scheme
      the_stokes_level.update_unsteady_rhs(
        vec_stokes_rhs,
        time_stepping_laminar.theta,
        vec_stokes_sol_1,
        vec_stokes_sol_2
      );

      // synchronize to obtain a type-1 vector
      vec_stokes_rhs.sync_0();

      // stop watch
      laminar_rhs.stop();

      // start watch
      laminar_filter_rhs.start();

      // apply filter
      the_stokes_level.filter_sys.filter_rhs(vec_stokes_rhs);

      // stop watch
      laminar_filter_rhs.stop();

      // ============================================================================================
      // Step 3: solve Navier-Stokes equations
      // --------------------------------------------------------------------------------------------

      // set theta for Burgers matrix assembly
      stokes_solver.theta = time_stepping_laminar.theta[0];

      // start watch
      laminar_solve.start();

      // solve the actual Navier-Stokes system for this time-step
      if(!stokes_solver.solve_nonlinear_turbulent(vec_stokes_sol, vec_stokes_rhs, vec_sol_k.local(), vec_sol_e.local()))
      {
        comm.print("\nStokes Solver didn't succeed!\n");
        return;
      }

      // stop watch
      laminar_solve.stop();

      // ============================================================================================
      // Step 4: write VTK file if desired
      // --------------------------------------------------------------------------------------------

      // write VTK file if desired
      if (write_laminar)
        vtk_writer_laminar.write_registered(time_stepping_laminar);

      // ============================================================================================
      // Step 5: update Stokes vectors
      // --------------------------------------------------------------------------------------------

      // shift the Stokes solution vectors
      vec_stokes_sol_3.copy(vec_stokes_sol_2); // u_{k-3}
      vec_stokes_sol_2.copy(vec_stokes_sol_1); // u_{k-2}
      vec_stokes_sol_1.copy(vec_stokes_sol);   // u_{k-1}
      check_point.save_checkpoint(time_stepping_laminar);
    }
    time_stepping_laminar.finish_loop();

    // Print out times
    comm.print(
      stringify("\nLaminar phase timings (s):\n\n  Extrapolate Sol:        ") + stringify(laminar_extrapolate.elapsed())
      + "\n  Filter Sol:             " + stringify(laminar_filter_sol.elapsed())
      + "\n  RHS Assembly:           " + stringify(laminar_rhs.elapsed())
      + "\n  Filter RHS:             " + stringify(laminar_filter_rhs.elapsed())
      + "\n  Solve Navier-Stokes:    " + stringify(laminar_solve.elapsed()) + "\n"
    );
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Turbulent Phase
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Print message
  comm.print(String("\n") + String(100u, '=') + "\n");
  comm.print("Solving turbulent Navier-Stokes system...\n");

  // Update Reynolds number for k-epsilon
  nu = params.u_mean*bnd_info.get_l() / params.Re_turbulent;
  stokes_solver.nu = nu;

  comm.print(
    stringify("\nTurbulent Reynolds Information:\n\n")
    + "  nu_turb         : " + stringify(nu) + "\n"
    + "  Re_turb         : " + stringify(params.Re_turbulent) + "\n"
    + "  u_mean (turb)   : " + stringify(params.u_mean) + "\n\n"
  );

  // StopWatches for time information
  StopWatch turb_extrapolate, turb_filter_sol;
  StopWatch turb_rhs, turb_neumann, turb_filter_rhs, turb_solve;
  StopWatch turb_k_solve, turb_e_solve;

  //         "       1      0.1000000000 [  1.00%]       0.000   5: 1.060e-01 > 3.824e-11 :   15;"
  comm.print("    Step       Time           Done       Runtime  NL  Def-Init    Def-Final     MG");
  comm.print("----------------------------------------------------------------------------------");

  time_stepping.begin_loop();
  vtk_writer.write_registered(time_stepping);
  stokes_solver.enable_turbulence = true;

  // create values for defect calculation
  DataType defect_old_k = 0.0;
  DataType defect_new_k = 0.0;
  DataType defect_old_e = 0.0;
  DataType defect_new_e = 0.0;

  // for Neumann we need to start with a smaller Re and then increase it after solving
  // for Dirichlet it is not necessary
  IndexType num_inc = 1;
  if(params.problem =="neumann")
    num_inc = 51;

  // increasing Re loop
  for(IndexType m = 0; m < num_inc; m++)
  {
    // read in old solutions (with smaller Re) for Neumann problem
    /*if (m > 0)
    {
      the_stokes_level.base_splitter_sys.split_read_from(vec_stokes_sol, "vec.u");
      the_turb_level.base_splitter_sys.split_read_from(vec_sol_k, "vec.k");
      the_turb_level.base_splitter_sys.split_read_from(vec_sol_e, "vec.e");
    }*/

    while(time_stepping.advance_step())
    {
      comm.print(String("\n") + String(100u, '=') + "\n");

      // Counter for setting k, e
      IndexType count_e = 0;
      IndexType count_k = 0;
      IndexType count_nu = 0;

      vec_sol_k_old.copy(vec_sol_k);
      vec_sol_e_old.copy(vec_sol_e);

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Solve u
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // ============================================================================================
      // Step 0: extrapolate previous Navier-Stokes solution to obtain an initial guess
      // --------------------------------------------------------------------------------------------

      // start watch
      turb_extrapolate.start();

      // initialize solution for this time-step by extrapolating the previous ones
      the_stokes_level.extrapolate_sol(
        vec_stokes_sol,
        time_stepping.expolc,
        vec_stokes_sol_1,
        vec_stokes_sol_2,
        vec_stokes_sol_3
      );

      // stop watch
      turb_extrapolate.stop();

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Outer coupling loop
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      for(IndexType k = 0; k < params.coupling_steps; k++)
      {
        comm.print("Timestep: " + stringify(time_stepping.time_step) + "\n");
        comm.print("Current coupling iteration: " + stringify(k + 1) + "\n");

        // ============================================================================================
        // Step 1: Filter velocity, k and e
        // --------------------------------------------------------------------------------------------

        // apply filter
        the_stokes_level.filter_sys.filter_sol(vec_stokes_sol);

        if(params.problem == "dirichlet")
        {
          // clear filter
          filter_k_wall.local().clear();
          filter_e_wall.local().clear();

          // compute u_tau, k_wall and e_wall values according to eq (18)
          for(Index i = 0; i < k_wall.local().size(); ++i)
          {
            u_tau = Math::max(Math::pow(params.C_mu, 0.25) * Math::sqrt(vec_sol_k.local()(i)), vec_stokes_sol.local().first()(i).norm_euclid()/11.06);
            k_wall.local()(i, pow(u_tau, 2.0) / Math::sqrt(params.C_mu));
            e_wall.local()(i, pow(u_tau, 3.0) / (0.41 * wall_distance));
          }

          // assemble wall filter
          unit_asm_wall.assemble(filter_k_wall.local(), domain_control.front()->space_turb, k_wall.local());
          unit_asm_wall.assemble(filter_e_wall.local(), domain_control.front()->space_turb, e_wall.local());

          // format left wall vectors
          left_wall_val_k.format();
          left_wall_val_e.format();

          // evaluate k_sol_func and e_sol_func into a vector
          Assembly::Interpolator::project(left_wall_val_k.local(), bnd_info.k_sol_func(), domain_control.front()->space_turb);
          Assembly::Interpolator::project(left_wall_val_e.local(), bnd_info.e_sol_func(), domain_control.front()->space_turb);

          // loop over left wall indices and set filter value manually because we have a different value at the inflow
          if(have_bnd_l)
          {
            const FEAT::Geometry::TargetSet& vertices_bnd_l = mesh_node->find_mesh_part("bnd:l")->get_target_set<0>();

            for(Index i = 0; i < vertices_bnd_l.get_num_entities(); ++i)
            {
              filter_k_wall.local().add(vertices_bnd_l[i], left_wall_val_k.local()(vertices_bnd_l[i]));
              filter_e_wall.local().add(vertices_bnd_l[i], left_wall_val_e.local()(vertices_bnd_l[i]));
            }
          }
        }

        if(params.problem == "neumann")
        {
          // clear filter
          filter_k_neumann.local().clear();
          filter_e_neumann.local().clear();

          // assemble inflow filter
          unit_asm.assemble(filter_k_neumann.local(), domain_control.front()->space_turb, bnd_info.k_sol_func());
          unit_asm.assemble(filter_e_neumann.local(), domain_control.front()->space_turb, bnd_info.e_sol_func());
        }

        // ============================================================================================
        // Step 2: assemble RHS vector for Navier-Stokes
        // --------------------------------------------------------------------------------------------

        // start watch
        turb_rhs.start();

        // format RHS (type 0)
        vec_stokes_rhs.format();
        the_stokes_level.update_unsteady_rhs(
          vec_stokes_rhs,
          time_stepping.theta,
          vec_stokes_sol_1,
          vec_stokes_sol_2
        );

        // synchronize to obtain a type-1 vector
        vec_stokes_rhs.sync_0();

        // start watch
        turb_neumann.start();

        // Neumann contributions
        vec_rhs_u_neumann.format();
        turb_trace_assembler.assemble_functional_vector_velocity(
          vec_rhs_u_neumann.local().first(),
          vec_sol_k.local(),
          vec_sol_e.local(),
          vec_stokes_sol.local().first(),
          domain_control.front()->space_turb,
          domain_control.front()->space_velo,
          params.cubature_name,
          1.0
        );

        // synchronize to obtain a type-1 vector
        vec_rhs_u_neumann.sync_0();

        // add to RHS
        vec_stokes_rhs.local().first().axpy(vec_rhs_u_neumann.local().first());

        // stop watch
        turb_neumann.stop();

        // start watch
        turb_filter_rhs.start();

        // apply filter
        the_stokes_level.filter_sys.filter_rhs(vec_stokes_rhs);

        // stop watches
        turb_filter_rhs.stop();
        turb_rhs.stop();

        // ============================================================================================
        // Step 3: solve Navier-Stokes equations
        // --------------------------------------------------------------------------------------------

        // set theta for Burgers matrix assembly
        stokes_solver.theta = time_stepping.theta[0];

        // start watch
        turb_solve.start();

        // solve the actual Navier-Stokes system for this time-step
        if(!stokes_solver.solve_nonlinear_turbulent(
          vec_stokes_sol,
          vec_stokes_rhs,
          vec_sol_k.local(),
          vec_sol_e.local()))
        {
          comm.print("\nERROR: Stokes solver failed!\n");
          return;
        }

        // stop watch
        turb_solve.stop();

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Solve k
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Implicit Euler-Scheme:
        // (M + delta_t * A_{t}) * k_{t} = M * k_{t-1} + delta_t * f_k_{t}

        // Format vectors and matrices
        vec_rhs_k.local().format();
        vec_f_k.format();
        matrix_k.local().format();
        vec_rhs_e.format();
        vec_f_e.local().format();
        matrix_e.local().format();
        vec_p_k.format(0.0);
        vec_p_k_bnd.format(0.0);
        vec_p_k_neumann.format(0.0);

        // Compute: vec_rhs_k = M * k_{t-1}
        matrix_mass_k.local().apply(vec_rhs_k.local(), vec_sol_k_old.local());

        // Update matrix_k, vec_f_k
        domain_control.front()->domain_asm.assemble(turb_asm_job);

        // Compute: vec_rhs_k = M * k_{t-1} + delta_t * f_k_{t}
        vec_rhs_k.axpy(vec_f_k, time_stepping.delta_t);

        // synchronize to obtain a type-1 vector
        vec_rhs_k.sync_0();

        // Compute: matrix_left_k = M + delta_t*matrix_k
        matrix_left_k.local().copy(matrix_mass_k.local());

        // add to RHS
        matrix_left_k.local().axpy(matrix_k.local(), time_stepping.delta_t);

        // apply filter
        if(params.problem == "dirichlet")
        {
          filter_k_wall.filter_sol(vec_sol_k);
          filter_k_wall.filter_rhs(vec_rhs_k);
        }
        if(params.problem == "neumann")
        {
          filter_k_neumann.filter_sol(vec_sol_k);
          filter_k_neumann.filter_rhs(vec_rhs_k);
        }

        // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        // start watch
        turb_k_solve.start();

        // Solve the system for k
        solver_k->init_numeric();
        if(params.problem == "dirichlet")
        {
          auto status_k = Solver::solve(*solver_k, vec_sol_k, vec_rhs_k, matrix_left_k, filter_k_wall);
          if(!Solver::status_success(status_k))
          {
            comm.print("Solver k didn't succeed!\n");
            break;
          }
        }
        if(params.problem == "neumann")
        {
          auto status_k = Solver::solve(*solver_k, vec_sol_k, vec_rhs_k, matrix_left_k, filter_k_neumann);
          if(!Solver::status_success(status_k))
          {
            comm.print("Solver k didn't succeed!\n");
            comm.print("Final Reynoldsnumber " + stringify(params.Re_turbulent) + "\n");
            break;
          }
        }

        // release the solver
        solver_k->done_numeric();

        // stop watch
        turb_k_solve.stop();

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Solve e
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // Implicit Euler-Scheme:
        // (M + delta_t * A_{t}) * e_{t} = M * e_{t-1} + delta_t * f_e_{t}
        // All matrices and vectors have already been assembled

        // Compute: vec_rhs_e = M * e_{t-1}
        matrix_mass_e.apply(vec_rhs_e, vec_sol_e_old);

        // synchronize to obtain a type-0 vector
        vec_rhs_e.from_1_to_0();

        // Compute: vec_rhs_e = M * e_{t-1} + delta_t * f_e_{t}
        vec_rhs_e.axpy(vec_f_e, time_stepping.delta_t);

        if(params.problem == "neumann")
        {
          // Compute: vec_rhs_e = M * e_{t-1} + delta_t * f_e_{t} + delta_t * eq (27)
          vec_rhs_e_neumann.format();

          // synchronize to obtain a type-0 vector
          vec_rhs_e_neumann.from_1_to_0();
          turb_trace_assembler.assemble_functional_vector_dissipation(
            vec_rhs_e_neumann.local(),
            vec_sol_k.local(),
            vec_sol_e.local(),
            vec_stokes_sol.local().first(),
            domain_control.front()->space_turb,
            domain_control.front()->space_velo,
            params.cubature_name,
            DataType(time_stepping.delta_t)
          );

          // synchronize to obtain a type-1 vector
          vec_rhs_e_neumann.sync_0();

          // add to RHS
          vec_rhs_e.axpy(vec_rhs_e_neumann);
        }

        // synchronize to obtain a type-1 vector
        vec_rhs_e.sync_0();

        // Compute: matrix_left_e = M + delta_t * matrix_e
        matrix_left_e.local().copy(matrix_mass_e.local());
        matrix_left_e.local().axpy(matrix_e.local(), time_stepping.delta_t);

        // apply filter
        if(params.problem =="dirichlet")
        {
          filter_e_wall.filter_sol(vec_sol_e);
          filter_e_wall.filter_rhs(vec_rhs_e);
        }
        if(params.problem =="neumann")
        {
          filter_e_neumann.filter_sol(vec_sol_e);
          filter_e_neumann.filter_rhs(vec_rhs_e);
        }

        // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

        // start watch
        turb_e_solve.start();

        // solve for e
        solver_e->init_numeric();
        if(params.problem == "dirichlet")
        {
          auto status_e = Solver::solve(*solver_e, vec_sol_e, vec_rhs_e, matrix_left_e, filter_e_wall);
          if(!Solver::status_success(status_e))
          {
            comm.print("Solver e didn't succeed!\n");
            break;
          }
        }
        if(params.problem == "neumann")
        {
          auto status_e = Solver::solve(*solver_e, vec_sol_e, vec_rhs_e, matrix_left_e, filter_e_neumann);
          if(!Solver::status_success(status_e))
          {
            comm.print("Solver e didn't succeed!\n");
            comm.print("Final Reynoldsnumber " + stringify(params.Re_turbulent) + "\n");
            break;
          }
        }

        // release the solver
        solver_e->done_numeric();

        // stop watch
        turb_e_solve.stop();

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Set values of e, k to 1e-10 if they become negative
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        DataType* data_ptr_e = vec_sol_e.local().elements();
        DataType* data_ptr_k = vec_sol_k.local().elements();
        DataType* data_ptr_nu = nu_vec.local().elements();
        DataType* data_ptr_gamma = vec_gamma.local().elements();
        IndexType data_size = vec_sol_e.local().size();
        DataType l_star;
        count_e = 0;
        count_k = 0;
        count_nu = 0;

        for (IndexType i = 0; i < data_size; ++i)
        {
          // Set e positive
          if (data_ptr_e[i] < 1e-10)
          {
            data_ptr_e[i] = 1e-10;
            ++count_e;
          }

          // Set k positive
          if (data_ptr_k[i] < 1e-10)
          {
            data_ptr_k[i] = 1e-10;
            ++count_k;
          }

          // Compute nu_T according to eq (8)
          if (params.C_mu * data_ptr_k[i] * Math::sqrt(data_ptr_k[i]) < data_ptr_e[i] * bnd_info.get_l_max())
          {
            l_star = params.C_mu * data_ptr_k[i] * Math::sqrt(data_ptr_k[i]) / data_ptr_e[i];
          }
          else
          {
            l_star = bnd_info.get_l_max();
            ++count_nu;
          }
          data_ptr_nu[i] = Math::max(params.nu_min, l_star * Math::sqrt(data_ptr_k[i]));

          // Compute gamma according to eq (9)
          data_ptr_gamma[i] = params.C_mu * data_ptr_k[i] / data_ptr_nu[i];
        }

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Calculate defect
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        res_defect_k.format();
        res_defect_e.format();

        //res_defect_k = matrix_k * vec_sol_k
        matrix_k.apply(res_defect_k, vec_sol_k);
        matrix_e.apply(res_defect_e, vec_sol_e);

        // res_defect_k -= vec_rhs_k
        // synchronize to obtain a type-0 vector
        res_defect_k.from_1_to_0();
        res_defect_e.from_1_to_0();

        // subtract from res_defect_*
        res_defect_k.axpy(vec_f_k, DataType(-1));
        res_defect_e.axpy(vec_f_e, DataType(-1));

        // synchronize to obtain a type-1 vector
        res_defect_k.sync_0();
        res_defect_e.sync_0();

        //filter defect
        if(params.problem == "dirichlet")
        {
          filter_k_wall.filter_def(res_defect_k);
          filter_e_wall.filter_def(res_defect_e);
        }
        if(params.problem == "neumann")
        {
          filter_k_neumann.filter_def(res_defect_k);
          filter_e_neumann.filter_def(res_defect_e);
        }

        // dual -> primal
        res_defect_k_primal.format();
        res_defect_e_primal.format();
        res_defect_k_primal.component_product(vec_mass_l, res_defect_k);
        res_defect_e_primal.component_product(vec_mass_l, res_defect_e);

        comm.print("\nNormalized primal defects:\n");
        comm.print("  defect_k_primal : " + stringify(res_defect_k_primal.norm2() / Math::sqrt(DataType(domain_control.front()->space_turb.get_num_dofs()))));
        comm.print("  defect_e_primal : " + stringify(res_defect_e_primal.norm2() / Math::sqrt(DataType(domain_control.front()->space_turb.get_num_dofs()))));

        defect_new_k = res_defect_k.norm2() / Math::sqrt(DataType(domain_control.front()->space_turb.get_num_dofs()));
        defect_new_e = res_defect_e.norm2() / Math::sqrt(DataType(domain_control.front()->space_turb.get_num_dofs()));

        // stop when both defects stagnated
        if(Math::abs(defect_new_k - defect_old_k) < 1e-12 && Math::abs(defect_new_e - defect_old_e) < 1e-12)
        {
          comm.print("Defect k and e stagnated");
          return;
        }

        // only print defect in last coupling loop
        if(k == params.coupling_steps - 1)
        {
          comm.print("\nFinal defect report:\n");
          comm.print("  defect_k (new) : " + stringify(defect_new_k));
          comm.print("  defect_k (old) : " + stringify(defect_old_k));
          comm.print("");
          comm.print("  defect_e (new) : " + stringify(defect_new_e));
          comm.print("  defect_e (old) : " + stringify(defect_old_e));
          comm.print("");
          defect_old_k = defect_new_k;
          defect_old_e = defect_new_e;
        }

        // stop when defects are small enough
        if(defect_new_k < 1e-12 && defect_new_e < 1e-12)
        {
          comm.print("Defect k and e are sufficiently small — convergence reached.");
        }
      } // outer coupling loop

      // ============================================================================================
      // Step 4: update Stokes vectors
      // --------------------------------------------------------------------------------------------

      // shift the Stokes solution vectors
      vec_stokes_sol_3.copy(vec_stokes_sol_2); // u_{k-3}
      vec_stokes_sol_2.copy(vec_stokes_sol_1); // u_{k-2}
      vec_stokes_sol_1.copy(vec_stokes_sol);   // u_{k-1}

      // ============================================================================================
      // Step 5: write checkpoint if desired
      // --------------------------------------------------------------------------------------------

      // save checkpoint if desired
      check_point.save_checkpoint(time_stepping);

      // ============================================================================================
      // Step 6: print plot line
      // --------------------------------------------------------------------------------------------

      // print plot line if full plot mode is disabled
      if(!time_stepping.full_plot)
      {
        String plot_line = time_stepping.plot_line;
        plot_line += stokes_solver.plot_line;
        plot_line += check_point.plot_line;
        plot_line += vtk_writer.plot_line;
        comm.print(plot_line);
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Print time and solver information
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      comm.print(
        stringify(time_stepping.time_step).pad_front(5)
        + stringify("k").pad_front(10)
        + stringify_fp_fix(time, 5).pad_front(10)
        + stringify(solver_k->get_num_iter()).pad_front(5)
      );

      // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      comm.print(
        stringify(time_stepping.time_step).pad_front(5)
        + stringify("e").pad_front(10)
        + stringify_fp_fix(time, 5).pad_front(10)
        + stringify(solver_e->get_num_iter()).pad_front(5)
      );

      // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      // Print out counter
      comm.print(
        "\n  k was set to a positive value " + stringify(count_k) + " times.\n"
        + "  e was set to a positive value " + stringify(count_e) + " times.\n"
        + "  l_* was set to l_max = " + stringify(bnd_info.get_l_max()) + " : " + stringify(count_nu) + " times.\n"
      );

      // -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      // Print out times
      comm.print(
        stringify("\nTurbulent phase timings (s):\n\n")
        + "  Extrapolate Sol:        " + stringify(turb_extrapolate.elapsed()) + "\n"
        + "  RHS Assembly:           " + stringify(turb_rhs.elapsed()) + "\n"
        + "  Neumann Assembly:       " + stringify(turb_neumann.elapsed()) + "\n"
        + "  Filter RHS:             " + stringify(turb_filter_rhs.elapsed()) + "\n"
        + "  Solve Navier-Stokes:    " + stringify(turb_solve.elapsed()) + "\n"
        + "  E Solve:                " + stringify(turb_e_solve.elapsed()) + "\n"
        + "  K Solve:                " + stringify(turb_k_solve.elapsed()) + "\n"
      );

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Calculate primal vector P_k
      if(params.problem == "dirichlet")
      {
        vec_p_k.format(0.0);
        vec_p_k.component_product(vec_mass_l, vec_f_k);
      }
      if(params.problem == "neumann")
      {
        vec_p_k.component_product(vec_mass_l, vec_p_k);
        vec_p_k_bnd.component_product(vec_mass_l, vec_p_k_bnd);
        vec_p_k_neumann.component_product(vec_mass_l, vec_p_k_neumann);
      }

      // equation (23)
      surface_int_23.format();

      // synchronize to obtain a type-0 vector
      surface_int_23.from_1_to_0();
      turb_trace_assembler.assemble_functional_vector_velocity(
        surface_int_23.local().first(),
        vec_sol_k.local(),
        vec_sol_e.local(),
        vec_stokes_sol.local().first(),
        domain_control.front()->space_turb,
        domain_control.front()->space_velo,
        params.cubature_name,
        DataType(1.0)
      );

      // synchronize to obtain a type-1 vector
      surface_int_23.sync_0();

      if(params.problem == "neumann")
      {
        // equation (27)
        surface_int_27.format();

        // synchronize to obtain a type-0 vector
        surface_int_27.from_1_to_0();
        turb_trace_assembler.assemble_functional_vector_dissipation(
          surface_int_27.local(),
          vec_sol_k.local(),
          vec_sol_e.local(),
          vec_stokes_sol.local().first(),
          domain_control.front()->space_turb,
          domain_control.front()->space_velo,
          params.cubature_name,
          DataType(1.0)
        );

        // synchronize to obtain a type-1 vector
        surface_int_27.sync_0();
      }

      // Vtk Export
      if(params.problem == "dirichlet")
      {
        filter_k_wall.filter_sol(vec_sol_k);
        filter_e_wall.filter_sol(vec_sol_e);
      }
      if(params.problem =="neumann")
      {
        filter_k_neumann.filter_sol(vec_sol_k);
        filter_e_neumann.filter_sol(vec_sol_e);
      }

      // write VTK
      vtk_writer.write_registered(time_stepping);

      // flush cout buffer
      comm.print_flush();

      if(params.problem == "neumann")
      {
        // increase Re by 1000 after every 1000 steps (for Neumann)
        if (time_stepping.time_step % 1000 == 0)
        {
          // if (time_stepping.time_step % 1000 == 0 && params.Re_turbulent != 47625 && params.Re_turbulent != 47000)
          if(params.Re_turbulent < 46999.0)
          {
            /*the_stokes_level.base_splitter_sys.join_write_out(vec_stokes_sol, "vec.u");
            the_turb_level.base_splitter_sys.join_write_out(vec_sol_k, "vec.k");
            the_turb_level.base_splitter_sys.join_write_out(vec_sol_e, "vec.e");*/
            params.Re_turbulent += DataType(1000);
            nu = params.u_mean*bnd_info.get_l() / params.Re_turbulent;
            stokes_solver.nu = nu;
            comm.print("Increasing Reynoldsnumber to " + stringify(params.Re_turbulent) + "\n");
            break;
          }
          // make sure the last Re is 47625
          // else if (time_stepping.time_step % 1000 == 0 && m > 0 && params.Re_turbulent == 47000)
          else if(m > 0)
          {
            /*the_stokes_level.base_splitter_sys.join_write_out(vec_stokes_sol, "vec.u");
            the_turb_level.base_splitter_sys.join_write_out(vec_sol_k, "vec.k");
            the_turb_level.base_splitter_sys.join_write_out(vec_sol_e, "vec.e");*/
            params.Re_turbulent = DataType(47625);
            nu = params.u_mean*bnd_info.get_l() / params.Re_turbulent;
            stokes_solver.nu = nu;
            comm.print("Increasing Reynoldsnumber to " + stringify(params.Re_turbulent) + "\n");
            break;
          }
        }
      }
    } // time stepping loop
  } // increasing Reynoldsnumber loop

  // print Re information
  comm.print("Final Reynoldsnumber " + stringify(params.Re_turbulent) + "\n");
  time_stepping.finish_loop();

  // And release the solver
  solver_k->done_symbolic();
  solver_e->done_symbolic();

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // write joined solution vector if desired
  stokes_solver.save_joined_sol_vector(vec_stokes_sol);

  // release the stokes_solver
  stokes_solver.release();
} // run_k_epsilon_model


// Here's our main method
int main(int argc, char* argv[])
{
  using namespace Turb;

  static_assert(dim == 2, "\nThis application currently supports only 2D simulations.\n");

  Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // create world communicator
  Dist::Comm comm(Dist::Comm::world());

  // print number of processes
  print_pad(comm, "Number of MPI Processes", stringify(comm.size()));

  // ================================================================================================
  // SimpleArgParser
  // ------------------------------------------------------------------------------------------------

  // create argument parser
  SimpleArgParser args(argc, argv);

  // check command line arguments
  DomainControl::add_supported_args(args);
  StokesSolver::add_supported_args(args);
  TimeStepping::add_supported_args(args);
  CheckPoint::add_supported_args(args);
  VtkWriter::add_supported_args(args);
  ParameterSettings::add_supported_args(args);

  // no arguments given?
  if(args.num_args() <= 1)
  {
    comm.print("\nNo arguments given; supported arguments:\n");
    comm.print(args.get_supported_help());
    comm.print("\nAt least the problem has to be specified");
    return 0;
  }

  // check for unsupported options
  {
    auto unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      // print all unsupported options
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      // abort
      return -1;
    }
  }

  // ================================================================================================
  // Parsing
  // ------------------------------------------------------------------------------------------------

  // a stop watch for our total runtime measurement
  StopWatch watch_total;

  // start watch
  watch_total.start();

  // create domain control object and print info
  DomainControl domain_control(comm);
  domain_control.create_domain(args);
  domain_control.print_info();

  // create parameter setting object and print info
  ParameterSettings params(comm);
  params.parse_args(args);
  params.print_config();

  // create Stokes solver, set default parameters, parse arguments and print configuration
  StokesSolver stokes_solver(domain_control);
  stokes_solver.C_mu = params.C_mu;
  stokes_solver.parse_args(args);
  stokes_solver.print_config();

  // create time stepping, set default parameters, parse arguments and print configuration
  TimeStepping time_stepping(comm);
  time_stepping.t_max = params.t_max_turb;
  time_stepping.delta_t = params.delta_t_turb;
  time_stepping.bdf_type = 1;
  time_stepping.parse_args(args);
  time_stepping.print_config();

  TimeStepping time_stepping_laminar(comm);
  time_stepping_laminar.t_max = params.t_max_laminar;
  time_stepping_laminar.delta_t = params.delta_t_laminar;
  time_stepping_laminar.bdf_type = 1;
  time_stepping_laminar.parse_args(args);
  time_stepping_laminar.print_config();

  // create checkpoint, parse arguments and print configuration
  CheckPoint check_point(comm);
  check_point.parse_args(args);
  check_point.print_config();

  // Mesh name
  String mesh_file_path = domain_control.mesh_file_names.front();
  String mesh_file =
  mesh_file_path.split_by_charset("/\\").back().split_by_string(".").front();

  // create VTK writer, set default parameters, parse arguments and
  // print configuration
  VtkWriter vtk_writer(domain_control);
  vtk_writer.name_prefix = mesh_file + "_" + stringify(dim) + "d";
  vtk_writer.parse_args(args);
  vtk_writer.print_config();
  vtk_writer.stepping = params.vtk_stepping;

  bool write_laminar = true;
  VtkWriter vtk_writer_laminar(domain_control);
  if (write_laminar)
  {
    vtk_writer_laminar.name_prefix = mesh_file + "_" + stringify(dim) +
    "d" + "_laminar";
    vtk_writer_laminar.parse_args(args);
    vtk_writer_laminar.print_config();
    vtk_writer_laminar.stepping = params.vtk_stepping;
  }

  // ================================================================================================
  // Boundary Info Object
  // ------------------------------------------------------------------------------------------------

  if (mesh_file == "backfacestep")
  {
    run_k_epsilon_model<BackWardConstant<typename DomainControl::MeshNodeType>>(
      domain_control,
      params,
      comm,
      stokes_solver,
      time_stepping,
      time_stepping_laminar,
      check_point,
      vtk_writer,
      write_laminar,
      vtk_writer_laminar
    );
  }
  else
  {
    comm.print("ERROR: Unknown mesh: " + mesh_file);
    Runtime::abort();
  }

  // make sure everyone is finished
  comm.barrier();

  // stop watch
  watch_total.stop();

  // print final statistics
  comm.print(String("\n") + String(100u, '=') + "\n");
  print_pad(comm, "Total Runtime", watch_total.elapsed_string().pad_front(10) + " sec {Balance} [Fraction]");
  domain_control.print_runtime(watch_total.elapsed());
  time_stepping.print_runtime(watch_total.elapsed());
  stokes_solver.print_runtime(watch_total.elapsed());
  check_point.print_runtime(watch_total.elapsed());
  vtk_writer.print_runtime(watch_total.elapsed());
  print_memory_usage(comm);

  // okay, we're done
  return 0;
} // END main
