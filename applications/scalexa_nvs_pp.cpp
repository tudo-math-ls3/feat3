// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// PURPOSE
// -------
// This application implements a 2D/3D unsteady incompressible Navier-Stokes solver, which uses
// Crank-Nicolson time discretization and a van Kan pressure projection scheme, which separates
// the solution of the coupled system into a disjoint velocity Burgers system solution step and
// a pressure Poisson system solution steps. The velocity Burgers system is by default solved by
// a BiCGStab-Jacobi solver and the pressure Poisson system is solved by a PCG-BoomerAMG solver.
//
// The benchmark that is solved by this application is the 2D/3D flow-around-a-cylinder problem
// with unsteady inflow boundary conditions on the time interval [0,8], which is described here:
//
// https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/flow/dfg_benchmark3_re100.html
//
// USAGE
// -----
// run this application with MPI with an arbitrary number of MPI ranks
//
// command line options:
//
// --mesh <path-to-meshfile>
// MANDATORY: This command line option specifies the path to the mesh file to be read in.
// This application was written to support pretty much every mesh file that discretizes the
// 2D/3D flow-around-a-cylinder/square domains, especially the following mesh-files:
// * flowbench_c2d_00_quad_130.xml: 2D flow-around-a-circle, 130 quads
// * flowbench_c2d_02_quad_48.xml:  2D flow-around-a-circle, 48 quads
// * flowbench_c2d_03_quad_64.xml:  2D flow-around-a-circle, 64 quads
// * flowbench_c2d_04_quad_144.xml: 2D flow-around-a-circle, 144 quads
// * flowbench_q2d_02_quad_48.xml:  2D flow-around-a-square, 48 quads
// * flowbench_q2d_03_quad_64.xml:  2D flow-around-a-square, 64 quads
// * flowbench_c3d_02_hexa_192.xml: 3D flow-around-a-cylinder, 192 hexahedra
// * flowbench_c3d_03_hexa_256.xml: 3D flow-around-a-cylinder, 256 hexahedra
// * flowbench_c3d_04_hexa_384.xml: 3D flow-around-a-cylinder, 384 hexahedra
// * flowbench_q3d_03_hexa_256.xml: 3D flow-around-a-cuboid, 256 hexahedra
//
// --level <level-max>
// MANDATORY: This command line option specifies the maximum refinement level <level-max>, which
// defines how many times the input mesh has to be refined regularly to obtain the final mesh that
// the PDE is to be solved on. It is recommended to set the refinement level to 3 or higher
// (for any of the aforementioned meshes) to obtain useful results.
//
// --steps <N>
// MANDATORY: Specifies the total number of time steps to perform for the entire time interval.
// It is recommended to set the number of time steps to 2000 or higher to obtain useful results.
//
// --time-max <T>
// Optional: Specifies the end point T of the simulation time interval (0,T).
// If not given, then T is set to 8.
//
// --nu <nu>
// Optional: Specifies the viscosity parameter nu.
// If not given, then nu is set to 0.001.
//
// --v-max <vmax>
// Optional: Specifies the maximum inflow velocity.
// If not given, then vmax is set to 1.5 in 2D and to 2.25 in 3D, which corresponds to Re=100.
//
// --inflow <meshpart-name>
// Optional: Specifies the name(s) of the inflow boundary meshpart(s).
// If not given, then 'bnd:l' is used as the inflow boundary.
//
// --outflow <meshpart-name>
// Optional: Specifies the name(s) of the outflow boundary meshpart(s).
// If not given, then 'bnd:r' is used as the outflow boundary.
//
// --obstacle <meshpart-name>
// Optional: Specifies the name(s) of the obstacle boundary meshpart(s).
// If not given, then 'bnd:c' is used as the obstacle boundary.
//
// --sigma <sigma>
// Optional: Specifies the scaling factor for the "reactive" pressure preconditioner.
// If not given, then sigma is set to 0.
//
// --upsam <uspam>
// Optional: Specifies the streamline diffusion stabilization parameter.
// If not given, then upsam is set to 0.01.
//
// --vtk [<filename> [<stepping>]]
// Optional: Specifies that the application should write out PVTU/VTU files during the simulation
// to visualize the solutions. If given, the optional <stepping> parameter specifies that VTU files
// should only be written out every <stepping> time steps to reduce the total number of VTU files
// generated. If <stepping> is not given, then VTU files are written out every 10 time steps.
//
// --dump
// Optional: Enables additional debug output about the partitioning.
//

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/global/alg_dof_parti.hpp>
#include <kernel/global/pmdcdsc_matrix.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/burgers_assembly_job.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/hypre.hpp>
#include <kernel/solver/superlu.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/ilu_precond.hpp>
#include <kernel/solver/frosch.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/fgmres.hpp>

#include <control/domain/unit_cube_domain_control.hpp>
#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/statistics.hpp>

#include <vector>

#ifndef FEAT_HAVE_MPI
#error This application must be configured and compiled with MPI support enabled.
#endif

#ifndef SCALEXA_DIM
#error You have to define the dimension of your app
#define SCALEXA_DIM 2
#endif

using DataType = double;
using IndexType = FEAT::Index;

namespace ScalexaNavierStokes
{
  using namespace FEAT;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Use the following typedefs to define the shape type of the mesh as well as the finite element space

  // select the element shape to be used; must be one of:
#if SCALEXA_DIM == 3
  typedef Shape::Hexahedron ShapeType;  // 3D Hexahedral Elements
#else
  typedef Shape::Quadrilateral ShapeType; // 2D Quadrilateral Elements
#endif

  // define the mesh type; this is always
  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  // define the trafo type; this is usually
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  // select an inf-sup stable finite element space pair to be used:
  typedef Space::Lagrange2::Element<TrafoType>       SpaceVeloType; // P2 or Q2 finite elements
  typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType; // P1dc finite elements

  // define the domain level type; for this application it is typically
  typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;

  // define our arch types
  typedef double DataType; // usually double
  typedef Index IndexType; // usually 64 bit unsigned int

  // get the dimension
  static constexpr int dim = ShapeType::dimension;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // this main function is called by the actual main function at the end of this source file
  void main(int argc, char* argv[])
  {
    // stop watch for time measurement
    StopWatch watch_total;
    watch_total.start();

    // create world communicator (this is our MPI_Comm wrapper object)
    Dist::Comm comm(Dist::Comm::world());
    comm.print(String(100u, '*'));

    // dump system call
    {
      String s("Arguments: ");
      s.append(argv[0]);
      for(int i(1); i < argc; ++i)
        s.append(" ").append(argv[i]);
      comm.print(s);
    }

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("mesh");
    args.support("level");
    args.support("steps");
    args.support("time-end");
    args.support("inflow");
    args.support("outflow");
    args.support("obstacle");
    args.support("v-max");
    args.support("sigma");
    args.support("upsam");
    args.support("dump");
    args.support("vtk");

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // FROSch
    Solver::Trilinos::add_supported_fpl_args(args);

    args.support("krylov_dim");
    args.support("krylov_tol_rel");
    args.support("krylov_tol_abs");
    args.support("krylov_inner_res_scale");
    args.support("krylov_maxit");
    args.support("krylov_minit");
    args.support("krylov_info");
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      // print all unsupported options to cerr
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());

      // abort
      FEAT::Runtime::abort();
    }

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh <meshfile>' is missing!");
      FEAT::Runtime::abort();
    }
    if(args.check("level") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
      FEAT::Runtime::abort();
    }
    if(args.check("steps") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--steps <steps>' is missing!");
      FEAT::Runtime::abort();
    }

    // viscosity parameter nu = 1/RE
    DataType nu = 1E-3;

    // theta scheme parameter (1.0 = backward Euler, 0.5 = Crank-Nicolson)
    DataType theta = 0.5;

    // scaling factor for diffusion Schur-complement preconditioner (pressure mass solver)
    DataType sigma = 0.0;

    // streamline diffusion stabilization parameter
    DataType upsam = 0.01;

    // maximum inflow velocity: 1.5 in 2D, 2.25 in 3D
    DataType max_inflow_velo = DataType(0.75) * DataType(dim);

    // simulation time end point
    DataType time_end = 8.0;
    if(args.parse("time-end", time_end) < 0)
    {
      comm.print(std::cerr, "ERROR: Failed to parse number of maximum time!");
      FEAT::Runtime::abort();
    }

    // number of time-steps
    Index num_steps = 1000;
    if(args.parse("steps", num_steps) < 0)
    {
      comm.print(std::cerr, "ERROR: Failed to parse number of time-steps!");
      FEAT::Runtime::abort();
    }

    args.parse("sigma", sigma);
    args.parse("upsam", upsam);
    args.parse("v-max", max_inflow_velo);

    Index vtk_step(10u);
    String vtk_name = "scalexa_nvs_pp_02";
    if(args.check("vtk") < 0)
      vtk_step = Index(0);
    if(args.parse("vtk", vtk_name, vtk_step) < 0)
    {
      comm.print(std::cerr, "ERROR: Failed to parse VTK filename or stepping!");
      FEAT::Runtime::abort();
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    std::shared_ptr<Solver::Trilinos::FROSchParameterList> params = std::make_shared<Solver::Trilinos::FROSchParameterList>(comm, dim, Solver::Trilinos::FROSchParameterList::PRESSUREPOISSON);
    params->parse_args(args);
    params->create_core();
    params->print();

    Index krylov_dim = 10;
    DataType tol_rel = 1E-8;
    DataType tol_abs = 1E-8;
    DataType inner_res_scale = 0;
    Index krylov_maxit = 1000;
    Index krylov_minit = 10;
    bool krylov_info = false;
    args.parse("krylov_dim", krylov_dim);
    args.parse("krylov_tol_rel", tol_rel);
    args.parse("krylov_tol_abs", tol_abs);
    args.parse("krylov_inner_res_scale", inner_res_scale);
    args.parse("krylov_maxit", krylov_maxit);
    args.parse("krylov_minit", krylov_minit);
    args.parse("krylov_info", krylov_info);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    StopWatch watch_domain_setup;
    watch_domain_setup.start();

    // create our domain control
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, false);

    // parse command line arguments
    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);

    // create domain partitioning; this is where the magic happens
    domain.create(args.query("mesh")->second);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // make sure that the desired level was actually chosen
    if(domain.get_desired_level_max() != domain.front()->get_level_index())
    {
      comm.print("ERROR: Domain control could not provide a valid partitioning for desired level " + domain.format_desired_levels());
      Runtime::abort();
    }

    // dump domain info if desired
    if(args.check("dump") >= 0)
    {
      comm.print("\nDomain Layers:");
      comm.allprint(domain.dump_layers());
      comm.print("\nDomain Layer Levels:");
      comm.allprint(domain.dump_layer_levels());
      comm.print("\nDomain Virtual Levels:");
      comm.allprint(domain.dump_virt_levels());
    }

    // in case that we need to write VTK files, we'll also need a refined mesh for the output
    std::unique_ptr<typename DomainLevelType::MeshNodeType> refined_mesh_node;
    if(vtk_step > Index(0))
      refined_mesh_node = domain.front()->get_mesh_node()->refine_unique();

    watch_domain_setup.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    StopWatch watch_system_asm;
    watch_system_asm.start();

    // define our system level
    typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim, DataType, IndexType> SystemLevelType;

    // get the finest domain level
    DomainLevelType& the_domain_level = *domain.front();

    // get the list of all boundary mesh-parts in the mesh
    std::set<String> all_meshpart_names, bnd_meshpart_names;
    for(const String& name : the_domain_level.get_mesh_node()->get_mesh_part_names(true))
    {
      all_meshpart_names.insert(name);
      if(name.starts_with("bnd:"))
        bnd_meshpart_names.insert(name);
    }

    // determine which boundary mesh-parts are inflow and which are noflow
    std::set<String> boundary_inflow, boundary_noflow, boundary_outflow, boundary_obstacle;

    // mixed Dirichlet-Neumann boundary conditions
    std::deque<String> desired_inflow_partnames;
    if(args.check("inflow") > 0)
      desired_inflow_partnames = args.query("inflow")->second;
    else
      desired_inflow_partnames.push_back("bnd:l");
    for(const String& name : desired_inflow_partnames)
    {
      if(all_meshpart_names.find(name) == all_meshpart_names.end())
      {
        comm.print(String("ERROR: Inflow meshpart '") + name + String("' not found!"));
        comm.print(String("Meshpart names for the given mesh are: '") + stringify_join(all_meshpart_names, "', '") + String("'"));
        Runtime::abort();
      }
      boundary_inflow.insert(name);
      auto it = bnd_meshpart_names.find(name);
      if(it != bnd_meshpart_names.end())
        bnd_meshpart_names.erase(it);
    }

    std::deque<String> desired_outflow_partnames;
    if(args.check("outflow") > 0)
      desired_outflow_partnames = args.query("outflow")->second;
    else
      desired_outflow_partnames.push_back("bnd:r");
    for(const String& name : desired_outflow_partnames)
    {
      if(all_meshpart_names.find(name) == all_meshpart_names.end())
      {
        comm.print(String("ERROR: Outflow meshpart '") + name + String("' not found!"));
        comm.print(String("Meshpart names for the given mesh are: '") + stringify_join(all_meshpart_names, "', '") + String("'"));
        Runtime::abort();
      }
      boundary_outflow.insert(name);
      auto it = bnd_meshpart_names.find(name);
      if(it != bnd_meshpart_names.end())
        bnd_meshpart_names.erase(it);
    }

    // all remaining boundaries are no-flow
    boundary_noflow = bnd_meshpart_names;

    // obstacle boundary: this is only required for body forces computation
    std::deque<String> desired_obstacle_partnames;
    if(args.check("obstacle") > 0)
      desired_obstacle_partnames = args.query("obstacle")->second;
    else
      desired_obstacle_partnames.push_back("bnd:c");
    for(const String& name : desired_obstacle_partnames)
    {
      if(all_meshpart_names.find(name) == all_meshpart_names.end())
      {
        comm.print(String("ERROR: Obstacle meshpart '") + name + String("' not found!"));
        comm.print(String("Meshpart names for the given mesh are: '") + stringify_join(all_meshpart_names, "', '") + String("'"));
        Runtime::abort();
      }
      boundary_obstacle.insert(name);
    }

    // get number of elements on finest mesh and sum up over all MPI processes
    Index num_elems_total = domain.front()->get_mesh().get_num_elements();
    comm.allreduce(&num_elems_total, &num_elems_total, 1u, Dist::op_sum);

    comm.print("\nProblem Configuration Summary:");
    comm.print("Dimension..................: " + stringify(dim));
    comm.print("Input Mesh File(s).........: " + stringify_join(args.query("mesh")->second, ", "));
    comm.print("Fine Mesh Level............: " + domain.format_chosen_levels());
    comm.print("Fine Mesh Elements.........: " + stringify(num_elems_total));
    comm.print("Maximum Inflow Velocity....: " + stringify_fp_fix(max_inflow_velo));
    comm.print("Viscosity Parameter nu.....: " + stringify_fp_sci(nu));
    comm.print("Simulation Time End........: " + stringify_fp_fix(time_end));
    comm.print("Number of Time Steps.......: " + stringify(num_steps));
    comm.print("Time Stepping Theta........: " + stringify_fp_fix(theta));
    comm.print("Diffusive Pressure Sigma...: " + stringify_fp_fix(sigma));
    comm.print("Streamline Diffusion Upsam.: " + stringify_fp_fix(upsam));


    if(boundary_inflow.empty())
    {
      comm.print("ERROR: no inflow boundary was specifies via --inflow <meshpart>");
      Runtime::abort();
    }
    else
      comm.print(String("Inflow  Boundaries.........: ") + stringify_join(boundary_inflow, " "));

    if(boundary_outflow.empty())
    {
      comm.print("ERROR: no outflow boundary was specifies via --outflow <meshpart>");
      Runtime::abort();
    }
    else
      comm.print(String("OutFlow Boundaries.........: ") + stringify_join(boundary_outflow, " "));

    if(boundary_noflow.empty())
      comm.print(String("No-Flow Boundaries.........: -N/A-"));
    else
      comm.print(String("No-Flow Boundaries.........: ") + stringify_join(boundary_noflow, " "));

    if(boundary_obstacle.empty())
      comm.print(String("Obstacle Boundary..........: -N/A-"));
    else
      comm.print(String("Obstacle Boundary..........: ") + stringify_join(boundary_obstacle, " "));

    // create a system level for the finest domain level
    SystemLevelType the_system_level;

    // define our cubature rule
    const String cubature("gauss-legendre:5");// + stringify(SpaceVeloType::local_degree+1));

    // get our FE spaces
    SpaceVeloType& space_velo = the_domain_level.space_velo;
    SpacePresType& space_pres = the_domain_level.space_pres;

    // prepare domain assembler
    the_domain_level.domain_asm.compile_all_elements();

    // assemble gate; this is the main communication structure
    the_system_level.assemble_gates(domain.front());

    // assemble velocity mass matrix
    the_system_level.assemble_velo_struct(space_velo);

    // assemble pressure gradient/velocity divergence matrices
    the_system_level.assemble_grad_div_matrices(space_velo, space_pres, cubature);

    // get the global matrix types for B and D (=B^T)
    typedef typename SystemLevelType::LocalMatrixBlockA LocalMatrixBlockA;
    typedef typename SystemLevelType::GlobalMatrixBlockA GlobalMatrixBlockA;
    typedef typename SystemLevelType::GlobalMatrixBlockB GlobalMatrixBlockB;
    typedef typename SystemLevelType::GlobalMatrixBlockD GlobalMatrixBlockD;
    typedef typename SystemLevelType::GlobalVeloFilter GlobalVeloFilter;
    typedef typename SystemLevelType::GlobalPresFilter GlobalPresFilter;
    typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    typedef typename SystemLevelType::LocalVeloVector LocalVeloVector;
    typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;
    typedef typename SystemLevelType::LocalPresVector LocalPresVector;
    typedef typename SystemLevelType::LocalScalarMatrix LocalScalarMatrix;
    typedef typename SystemLevelType::LocalScalarVector LocalScalarVector;

    // get our matrices
    GlobalMatrixBlockA& matrix_a = the_system_level.matrix_a; // velocity burgers matrix
    GlobalMatrixBlockB& matrix_b = the_system_level.matrix_b; // pressure gradient matrix
    GlobalMatrixBlockD& matrix_d = the_system_level.matrix_d; // velocity divergence matrix
    GlobalMatrixBlockA  matrix_m = matrix_a.clone();          // velocity mass matrix
    GlobalMatrixBlockA  matrix_l = matrix_a.clone();          // velocity laplace matrix
    GlobalMatrixBlockA  matrix_k = matrix_a.clone();          // velocity convection matrix
    LocalScalarMatrix   matrix_o;                             // pressure mass matrix

    // create four new vectors
    GlobalVeloVector vec_sol_v = matrix_d.create_vector_r(); // v_k
    GlobalVeloVector vec_pre_v = matrix_d.create_vector_r(); // v_{k-1}
    GlobalVeloVector vec_tmp_v = matrix_d.create_vector_r();
    GlobalVeloVector vec_rhs_v = matrix_d.create_vector_r(); // r_k

    GlobalPresVector vec_sol_p = matrix_b.create_vector_r(); // p_k
    GlobalPresVector vec_cor_p = matrix_b.create_vector_r();
    GlobalPresVector vec_rhs_p = matrix_b.create_vector_r();
    GlobalPresVector vec_tmp_p = matrix_b.create_vector_r();
    GlobalPresVector vec_div_v = matrix_b.create_vector_r();

    LocalVeloVector lvec_mod_v(space_velo.get_num_dofs(), DataType(0)); // (nu*L + K(v_{k-1}))*v_{k-1}
    LocalVeloVector lvec_mod_p(space_velo.get_num_dofs(), DataType(0)); // B*p_{k-1}
    LocalVeloVector lvec_bdfd(space_velo.get_num_dofs(), DataType(0)); // body forces defect vector

    LocalScalarVector lvec_obst(space_velo.get_num_dofs(), DataType(0)); // obstacle boundary characteristic vector

    // format vectors to zero
    vec_sol_v.format();
    vec_pre_v.format();
    vec_tmp_v.format();
    vec_rhs_v.format();
    vec_sol_p.format();
    vec_rhs_p.format();

    // get our filters
    GlobalVeloFilter& filter_v = the_system_level.filter_velo; // velocity unit filter
    GlobalPresFilter& filter_p = the_system_level.filter_pres; // pressure none filter (does nothing, actually)

    // loop over all inflow and noflow boundary parts
    Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow, unit_asm_noflow;
    for(const String& name : boundary_inflow)
    {
      const auto* mpart = the_domain_level.get_mesh_node()->find_mesh_part(name);
      if(mpart)
        unit_asm_inflow.add_mesh_part(*mpart);
    }
    for(const String& name : boundary_noflow)
    {
      const auto* mpart = the_domain_level.get_mesh_node()->find_mesh_part(name);
      if(mpart)
        unit_asm_noflow.add_mesh_part(*mpart);
    }
    unit_asm_inflow.assemble(filter_v.local(), space_velo);
    unit_asm_noflow.assemble(filter_v.local(), space_velo);

    // assemble temporary obstacle filter and filter obstacle characteristic vector
    {
      Assembly::UnitFilterAssembler<MeshType> unit_asm_obstacle;
      for(const String& name : boundary_obstacle)
      {
        const auto* mpart = the_domain_level.get_mesh_node()->find_mesh_part(name);
        if(mpart)
          unit_asm_obstacle.add_mesh_part(*mpart);
      }
      LAFEM::UnitFilter<DataType, IndexType> filter_obstacle_v(space_velo.get_num_dofs());
      unit_asm_obstacle.assemble(filter_obstacle_v, space_velo);
      filter_obstacle_v.get_filter_vector().format(DataType(1));
      filter_obstacle_v.filter_sol(lvec_obst);
    }

    // assemble velocity mass and laplace matrices
    {
      matrix_m.local().format();
      matrix_l.local().format();
      Assembly::Common::IdentityOperatorBlocked<dim> identity_op;
      Assembly::Common::LaplaceOperatorBlocked<dim> laplace_op;
      Assembly::assemble_bilinear_operator_matrix_1(the_domain_level.domain_asm, matrix_m.local(), identity_op, space_velo, cubature);
      Assembly::assemble_bilinear_operator_matrix_1(the_domain_level.domain_asm, matrix_l.local(), laplace_op, space_velo, cubature);
    }

    // assemble pressure mass matrix
    {
      Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_o, space_pres);
      Assembly::Common::IdentityOperator identity_op;
      matrix_o.format();
      Assembly::assemble_bilinear_operator_matrix_1(the_domain_level.domain_asm, matrix_o, identity_op, space_pres, cubature);
    }

    // compute global lumped velocity mass matrix
    GlobalVeloVector lumped_mass_v = matrix_m.lump_rows(true);
    //comm.print("\nVelocity Lumped Mass Min/Max: " + stringify_fp_sci(lumped_mass_v.min_abs_element()) + " / " + stringify_fp_sci(lumped_mass_v.max_abs_element()));

    // invert lumped velocity mass matrix
    GlobalVeloVector lumped_mass_inv_v = lumped_mass_v.clone();
    lumped_mass_inv_v.component_invert(lumped_mass_v, DataType(1));

    // apply velocity filters onto the inverse lumped mass matrix
    filter_v.filter_def(lumped_mass_inv_v);

    // unmap pressure evaluation points p_a and p_e
    Trafo::InverseMappingData<DataType, dim> point_iv_a, point_iv_e;
    {
      typedef Trafo::InverseMapping<TrafoType, DataType> InvMappingType;
      InvMappingType inv_mapping(space_velo.get_trafo());

      // reference pressure points
      typename InvMappingType::ImagePointType v_a, v_e;
      if(dim == 2)
      {
        v_a[0] = DataType(0.15);
        v_e[0] = DataType(0.25);
        v_a[1] = v_e[1] = DataType(0.2);
      }
      else
      {
        v_a[0] = DataType(0.45);
        v_e[0] = DataType(0.55);
        v_a[1] = v_e[1] = DataType(0.2);
        v_a[2] = v_e[2] = DataType(0.205);
      }

      // unmap points
      point_iv_a = inv_mapping.unmap_point(v_a, true);
      point_iv_e = inv_mapping.unmap_point(v_e, true);
    }


    watch_system_asm.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    comm.print("\nInitializing Global Pre-Multiplied Discontinuous Diagonal Schur-Complement Matrix...");

    StopWatch watch_pmdcdsc_init;
    watch_pmdcdsc_init.start();

    // define a global Pre-Multiplied Discontinuous Diagonal Schur-Complement Matrix
    Global::PMDCDSCMatrix<GlobalMatrixBlockB, GlobalMatrixBlockD> matrix_s(lumped_mass_inv_v, matrix_b, matrix_d);
    matrix_s.init();

    watch_pmdcdsc_init.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    comm.print("Initializing Velocity Burgers Solver...");

    StopWatch watch_solver_v_setup, watch_solver_v_apply;
    watch_solver_v_setup.start();

    // create BiCGStab-Jacobi solver and initialize it symbolically
    auto solver_v = Solver::new_bicgstab(matrix_a, filter_v, Solver::new_jacobi_precond(matrix_a, filter_v));

    // solver_v->set_tol_rel(1e-5);
    // solver_v->set_tol_abs(1e-10);
    solver_v->set_tol_rel(tol_rel);
    solver_v->set_tol_abs(tol_abs);
    solver_v->set_tol_abs_low(1e-14);
    solver_v->set_max_iter(1000);
    solver_v->init_symbolic();

    watch_solver_v_setup.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    comm.print("Initializing Pressure Poisson Solver...");

    StopWatch watch_solver_p_setup, watch_solver_p_apply;
    watch_solver_p_setup.start();

    // create BoomerAMG solver and initialize it
    //auto solver_p = Solver::new_boomeramg(matrix_s, filter_p);
    //auto solver_p = Solver::new_richardson(matrix_s, filter_p, 1.0, Solver::new_boomeramg(matrix_s, filter_p));
    //auto solver_p = Solver::new_pcg(matrix_s, filter_p, Solver::new_boomeramg(matrix_s, filter_p));
    //auto solver_p = Solver::new_pcg(matrix_s, filter_p, Solver::new_jacobi_precond(matrix_s, filter_p));
    auto frosch_precond = Solver::new_frosch(matrix_s, filter_p, *params);
    // auto solver_p  = Solver::new_gmres(matrix_s,
    //                                    filter_p,
    //                                    krylov_dim,
    //                                    inner_res_scale,
    //                                    frosch_precond);
    auto solver_p  = Solver::new_fgmres(matrix_s,
                                       filter_p,
                                       krylov_dim,
                                       inner_res_scale,
                                       frosch_precond);
    solver_p->set_tol_rel(tol_rel);
    solver_p->set_tol_abs(tol_abs);
    solver_p->set_tol_abs_low(1e-14);
    solver_p->set_min_iter(krylov_minit);
    solver_p->set_max_iter(krylov_maxit);

    // auto solver_p = Solver::new_bicgstab(matrix_s, filter_p, frosch_precond);
    // solver_p->set_tol_rel(tol_rel);
    // solver_p->set_tol_abs(tol_abs);
    // solver_p->set_tol_abs_low(1e-14);
    // solver_p->set_max_iter(krylov_maxit);

    if(krylov_info)
      solver_p->set_plot_mode(Solver::PlotMode::iter);

    // solver_p->set_tol_rel(1e-5);
    // solver_p->set_tol_abs(1e-10);
    // solver_p->set_tol_abs_low(1e-14);
    // solver_p->set_max_iter(10000);

    solver_p->init();

    watch_solver_p_setup.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    StopWatch watch_solver_o_setup, watch_solver_o_apply;

    comm.print("Initializing Pressure Mass Solver...");

    watch_solver_o_setup.start();

    // create a pressure mass solver: ILU(0) is an exact solver in our case
    auto solver_o = Solver::new_ilu_precond(PreferredBackend::generic, matrix_o, filter_p.local());

    if(sigma > DataType(0))
      solver_o->init();

    watch_solver_o_setup.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // some constants we'll need
    const DataType pi = Math::pi<DataType>();
    const DataType channel_d = DataType(0.41); // benchmark channel diameter

    comm.print("\nStarting Time-Stepping loop...");

    //         "     1   0.00500     6: 1.419e-08 > 1.132e-13 [   0.131] |   209: 3.139e-03 > 9.149e-11 [   1.787] | 7.84692e-06 |   0.1435862  -0.0000059   0.0848067"
    comm.print("  Step      Time  #VBS  Def-Init    Def-Final   Runtime  |  #PPS  Def-Init    Def-Final   Runtime  | Divergence  |   Drag        Lift        P-Diff");
    comm.print("------------------------------------------------------------------------------------------------------------------------------------------------------");

    StopWatch watch_time_loop, watch_compute_defect, watch_body_forces, watch_vtk_export;
    watch_time_loop.start();

    // beginning of time step loop
    bool quit_loop = false;
    for(Index time_step = 1; time_step <= num_steps; ++time_step)
    {
      // compute time step size and current time
      const DataType delta_t = time_end / DataType(num_steps);
      const DataType cur_time = delta_t * DataType(time_step);

      // our output line for this step
      String line, line_extra;
      line += stringify(time_step).pad_front(6);
      line += stringify_fp_fix(cur_time, 5, 10);

      // ==============================================================================================================

      // STEP 1: re-assemble velocity inflow filter of corresponding dimension
      watch_system_asm.start();
      if constexpr(dim == 2)
      {
        const DataType iota = max_inflow_velo * DataType(4) * Math::sin(pi * cur_time * DataType(0.125)) / Math::sqr(channel_d);
        auto inflow_func_2d = Analytic::create_lambda_function_vector_2d(
          [iota, channel_d] (DataType, DataType y) {return iota * y * (channel_d - y);},
          [] (DataType, DataType) {return DataType(0);});
        unit_asm_inflow.assemble(filter_v.local(), the_domain_level.space_velo, inflow_func_2d);
      }
      if constexpr(dim == 3)
      {
        const DataType iota = max_inflow_velo * DataType(16) * Math::sin(pi * cur_time * DataType(0.125)) / Math::sqr(Math::sqr(channel_d));
        auto inflow_func_3d = Analytic::create_lambda_function_vector_3d(
          [iota, channel_d] (DataType, DataType y, DataType z) {return iota * y * (channel_d - y) * z * (channel_d - z);},
          [] (DataType, DataType, DataType) {return DataType(0);},
          [] (DataType, DataType, DataType) {return DataType(0);});
        unit_asm_inflow.assemble(filter_v.local(), the_domain_level.space_velo, inflow_func_3d);
      }
      watch_system_asm.stop();

      // ==============================================================================================================

      // STEP 2: extrapolate velocity vector: w_k := 2*v_{k-1} - v_{k-2}
      if(time_step > Index(1))
      {
        vec_tmp_v.copy(vec_pre_v); // v_{k-2}
        vec_pre_v.copy(vec_sol_v); // v_{k-1}
        vec_sol_v.scale(vec_pre_v, DataType(2));
        vec_sol_v.axpy(vec_tmp_v, -DataType(1));
      }
      else
      {
        vec_pre_v.copy(vec_sol_v);
      }

      // apply filter to impose boundary conditions of current step
      filter_v.filter_sol(vec_sol_v);

      // ==============================================================================================================

      // STEP 3: assemble Burgers matrix: A := nu*L + theta*M + K(w_k)

      watch_system_asm.start();

      // A := M + theta*delta_t*(nu*L + K(w_k))
      matrix_a.local().copy(matrix_m.local());
      matrix_a.local().axpy(matrix_l.local(), nu * theta * delta_t);

      // assemble convective term K(v) by using Burgers assembly job
      Assembly::BurgersBlockedMatrixAssemblyJob<LocalMatrixBlockA, SpaceVeloType> burgers_matrix_job(
        matrix_a.local(), vec_sol_v.local(), space_velo, cubature);
      burgers_matrix_job.beta = theta * delta_t;
      if(upsam > DataType(0))
      {
        burgers_matrix_job.sd_nu = nu;
        burgers_matrix_job.sd_delta = upsam;
        burgers_matrix_job.set_sd_v_norm(vec_sol_v);
      }
      the_domain_level.domain_asm.assemble(burgers_matrix_job);

      watch_system_asm.stop();

      // ==============================================================================================================

      // STEP 4: build right hand side for momentum equation

      // r_k := M*v_{k-1} - (1 - theta)*delta_t*(nu*L + K(v_{k-1}))*v_{k-1} - delta_t*B*p_{k-1}
      //      = M*v_{k-1} - (1 - theta)*delta_t*mod_v_k - delta_t*mod_p_k
      matrix_m.local().apply(vec_rhs_v.local(), vec_pre_v.local());
      vec_rhs_v.local().axpy(lvec_mod_v, -(1.0-theta)*delta_t);
      vec_rhs_v.local().axpy(lvec_mod_p, -delta_t);
      vec_rhs_v.sync_0();

      filter_v.filter_rhs(vec_rhs_v);

      // ==============================================================================================================

      // STEP 5: solve velocity Burgers problem

      StopWatch watch_v;
      watch_v.start();

      // initialize solver numerically
      watch_solver_v_setup.start();
      solver_v->init_numeric();
      watch_solver_v_setup.stop();

      // apply velocity Burgers solver
      watch_solver_v_apply.start();
      Solver::Status status_v = solver_v->correct(vec_sol_v, vec_rhs_v);
      watch_solver_v_apply.stop();

      // release solver numerically
      solver_v->done_numeric();

      // print number of iterations and final defect
      line += stringify(solver_v->get_num_iter()).pad_front(6);
      line += ":";
      line += stringify_fp_sci(solver_v->get_def_initial(), 3, 10);
      line += " >";
      line += stringify_fp_sci(solver_v->get_def_final(), 3, 10);

      // analyze solver status
      switch(status_v)
      {
      case Solver::Status::success:
        // Yay!
        break;

      case Solver::Status::diverged:
        line += "\n\nERROR: velocity solver diverged; aborting!\n";
        quit_loop = true;
        break;

      case Solver::Status::aborted:
        line += "\n\nERROR: velocity solver aborted; aborting!\n";
        quit_loop = true;
        break;

      case Solver::Status::max_iter:
        // not necessarily an error; so we'll just issue a warning later
        line_extra += " ; velocity solver reached maximum iterations";
        break;

      case Solver::Status::stagnated:
        // not necessarily an error; so we'll just issue a warning later
        line_extra += " ; velocity solver stagnated";
        break;

      default:
        // we should never arrive here
        line += "\n\nERROR: velocity solver returned unknown status; aborting!\n";
        quit_loop = true;
        break;
      }

      // is this the end?
      if(quit_loop)
      {
        if(!line_extra.empty())
          line += line_extra;
        comm.print(line);
        break;
      }

      watch_v.stop();

      line += " [";
      line += watch_v.elapsed_string().pad_front(8);
      line += "]";

      // ==============================================================================================================

      // STEP 6: solve pressure Poisson problem

      StopWatch watch_p;
      watch_p.start();

      watch_compute_defect.start();

      // compute divergence defect
      vec_rhs_p.format();
      matrix_d.apply(vec_rhs_p, vec_sol_v, vec_rhs_p, DataType(1) / delta_t);

      // filter divergence defect (doesn't actually do anything in our case)
      filter_p.filter_def(vec_rhs_p);

      // compute defect norm
      watch_compute_defect.stop();

      if(sigma > DataType(0))
      {
        // apply pressure Mass solver
        watch_solver_o_apply.start();
        Solver::Status status_o = solver_o->apply(vec_cor_p.local(), vec_rhs_p.local());
        watch_solver_o_apply.stop();
        // this one should never fail
        XASSERT(status_o == Solver::Status::success);

        // update pressure solution
        vec_sol_p.axpy(vec_cor_p, nu * delta_t * sigma);
      }

      // apply pressure Poisson solver
      watch_solver_p_apply.start();
      Solver::Status status_p = solver_p->apply(vec_cor_p, vec_rhs_p);
      watch_solver_p_apply.stop();

      // print number of iterations and final defect
      line += " |";
      line += stringify(solver_p->get_num_iter()).pad_front(6);
      line += ":";
      line += stringify_fp_sci(solver_p->get_def_initial(), 3, 10);
      line += " >";
      line += stringify_fp_sci(solver_p->get_def_final(), 3, 10);

      // analyze solver status
      switch(status_p)
      {
      case Solver::Status::success:
        // Yay!
        break;

      case Solver::Status::diverged:
        line += "\n\nERROR: pressure solver diverged; aborting!\n";
        quit_loop = true;
        break;

      case Solver::Status::aborted:
        line += "\n\nERROR: pressure solver aborted; aborting!\n";
        quit_loop = true;
        break;

      case Solver::Status::max_iter:
        // not necessarily an error; so we'll just issue a warning later
        line_extra += " ; pressure solver reached maximum iterations";
        break;

      case Solver::Status::stagnated:
        // not necessarily an error; so we'll just issue a warning later
        line_extra += " ; pressure solver stagnated";
        break;

      default:
        // we should never arrive here
        line += "\n\nERROR: pressure solver return unknown status; aborting!\n";
        quit_loop = true;
        break;
      }

      // is this the end?
      if(quit_loop)
      {
        if(!line_extra.empty())
          line += line_extra;
        comm.print(line);
        break;
      }

      // update pressure solution
      vec_sol_p.axpy(vec_cor_p);

      watch_p.stop();

      line += " [";
      line += watch_p.elapsed_string().pad_front(8);
      line += "]";

      // ==============================================================================================================

      // STEP 7: update velocity with new pressure and pre-compute defect

      // v_k := w_k - theta*delta_t*M_l^{-1} * B * q_k
      matrix_b.apply(vec_tmp_v, vec_cor_p);
      vec_tmp_v.component_product(vec_tmp_v, lumped_mass_inv_v);
      filter_v.filter_cor(vec_tmp_v);
      vec_sol_v.axpy(vec_tmp_v, -theta*delta_t);

      // ==============================================================================================================

      // STEP 8: assemble final momentum defect of this time step (without the mass)

      watch_system_asm.start();

      // mod_v_k = nu*L*v_k + K(v_k)*v_k
      lvec_mod_v.format();
      matrix_l.local().apply(lvec_mod_v, vec_sol_v.local(), lvec_mod_v, nu);
      Assembly::BurgersBlockedVectorAssemblyJob<LocalVeloVector, SpaceVeloType> burgers_vector_job(
        lvec_mod_v, vec_sol_v.local(), vec_sol_v.local(), space_velo, cubature);
      burgers_vector_job.beta = DataType(1);
      the_domain_level.domain_asm.assemble(burgers_vector_job);

      // mod_p_k := B*p_k = B*p_{k-1} + B*q_k = mod_p_{k-1} + B*q_k
      matrix_b.local().apply(lvec_mod_p, vec_sol_p.local());

      watch_system_asm.stop();

      // ==============================================================================================================

      // STEP 9: compute divergence of velocity field

      watch_compute_defect.start();

      matrix_d.apply(vec_div_v, vec_sol_v);
      filter_p.filter_def(vec_div_v);

      const DataType div_v = vec_div_v.norm2();

      line += " |";
      line += stringify_fp_sci(div_v, 5, 12);

      watch_compute_defect.stop();

      // ==============================================================================================================

      // STEP 10: post-processing: compute drag & lift forces

      watch_body_forces.start();

      // compute body forces defect vector
      // lvec_bdfd.axpy(lvec_mod_v, lvec_mod_p);
      lvec_bdfd.copy(lvec_mod_p);
      lvec_bdfd.axpy(lvec_mod_v);

      // compute body forces
      Tiny::Vector<DataType, 3> body_forces(DataType(0));

      {
        Tiny::Vector<DataType, dim> body_frc_d(DataType(0));
        const Index n = vec_tmp_v.local().size();
        const Tiny::Vector<DataType, dim>* vt = lvec_bdfd.elements();
        const DataType* vc = lvec_obst.elements();
        for(Index i(0); i < n; ++i)
          body_frc_d.axpy(vc[i], vt[i]);
        comm.allreduce(body_forces.v, body_forces.v, 3u, Dist::op_sum);
        body_forces.template copy_n<dim>(body_frc_d);
      }

      const DataType bf_scale = (dim == 2 ?
        DataType(2) / (DataType(0.100)*Math::sqr(max_inflow_velo*(DataType(2)/DataType(3)))) : // = 2 / (rho * U^2 * D)
        DataType(2) / (DataType(0.041)*Math::sqr(max_inflow_velo*(DataType(4)/DataType(9))))); // = 2 / (rho * U^2 * D * H)

      // compute pressure difference at obstacle
      DataType pval_a = Assembly::DiscreteEvaluator::eval_fe_function(point_iv_a, vec_sol_p.local(), space_pres).mean_value_dist(comm);
      DataType pval_e = Assembly::DiscreteEvaluator::eval_fe_function(point_iv_e, vec_sol_p.local(), space_pres).mean_value_dist(comm);

      // print drag and lift coefficients
      line += " |";
      line += stringify_fp_fix(-bf_scale * body_forces[0], 7, 12); // drag coefficient
      line += stringify_fp_fix(-bf_scale * body_forces[1], 7, 12); // lift coefficient
      line += stringify_fp_fix(pval_a - pval_e, 7, 12);            // pressure difference

      watch_body_forces.stop();

      // ==============================================================================================================

      // STEP 11: write VTK file if desired

      if(!vtk_name.empty() && (vtk_step > Index(0)) && ((time_step % vtk_step) == 0))
      {
        watch_vtk_export.start();

        XASSERT(bool(refined_mesh_node));
        String vtk_filename = vtk_name + "." + stringify(time_step).pad_front(6, '0');
        line_extra += " ; writing VTK file '" + vtk_filename + "'";

        Geometry::ExportVTK<MeshType> vtk_exp(*refined_mesh_node->get_mesh());
        vtk_exp.add_vertex_vector("velocity", vec_sol_v.local());

        // compute time derivative of velocity
        // vec_tmp_v.axpy(vec_pre_v, vec_sol_v, -1.0);
        vec_tmp_v.copy(vec_sol_v);
        vec_tmp_v.axpy(vec_pre_v, -1.0);
        vec_tmp_v.scale(vec_tmp_v, delta_t);
        vtk_exp.add_vertex_vector("velocity_dt", vec_tmp_v.local());

        // project pressure onto refined mesh
        LocalPresVector prj_p, prj_q;
        Assembly::DiscreteCellProjector::project_refined(prj_p, vec_sol_p.local(), space_pres);
        Assembly::DiscreteCellProjector::project_refined(prj_q, vec_div_v.local(), space_pres);
        vtk_exp.add_cell_scalar("pressure", prj_p.elements());
        vtk_exp.add_cell_scalar("divergence", prj_q.elements());

        vtk_exp.write(vtk_filename, comm);

        watch_vtk_export.stop();
      }

      // ==============================================================================================================

      // print line for this time step
      if(!line_extra.empty())
        line += line_extra;
      comm.print(line);
    } // time stepping loop

    watch_time_loop.stop();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // release solvers
    if(sigma > DataType(0))
      solver_o->done();
    solver_p->done();
    solver_v->done_symbolic();

    // print elapsed runtime
    watch_total.stop();
    const double time_total = watch_total.elapsed();
    comm.print("\nSimulation Runtime Summary:");
    comm.print("Total Runtime.........................: " + watch_total.elapsed_string().pad_front(10) + " seconds");
    comm.print("Time for Domain Setup.................: " + watch_domain_setup.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_domain_setup.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for System Assembly..............: " + watch_system_asm.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_system_asm.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for PMDCDSC Matrix setup.........: " + watch_pmdcdsc_init.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_pmdcdsc_init.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for Velocity Burgers solver setup: " + watch_solver_v_setup.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_solver_v_setup.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for Velocity Burgers solver apply: " + watch_solver_v_apply.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_solver_v_apply.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for Pressure Poisson solver setup: " + watch_solver_p_setup.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_solver_p_setup.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for Pressure Poisson solver apply: " + watch_solver_p_apply.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_solver_p_apply.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for Pressure Mass solver setup...: " + watch_solver_o_setup.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_solver_o_setup.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for Pressure Mass solver apply...: " + watch_solver_o_apply.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_solver_o_apply.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for Body Forces computation......: " + watch_body_forces.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_body_forces.elapsed() / time_total, 2, 6) + "%]");
    comm.print("Time for VTK export...................: " + watch_vtk_export.elapsed_string().pad_front(10) + " seconds [" +
      stringify_fp_fix(100.0 * watch_vtk_export.elapsed() / time_total, 2, 6) + "%]");
  }
} // namespace ScalexaNavierStokes

int main(int argc, char* argv [])
{
  FEAT::Runtime::ScopeGuard scope_guard(argc, argv);
  try
  {
    ScalexaNavierStokes::main(argc, argv);
  }
  catch (const std::exception& exc)
  {
    std::cerr << "ERROR: unhandled exception: " << exc.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return 0;
}
