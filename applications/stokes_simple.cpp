// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// ====================================================================================================================
// Simple steady-state Stokes Example Application with fully-grown MPI-parallel Geometric Multigrid Solver
// --------------------------------------------------------------------------------------------------------------------
// This application is a simple steady-state Stokes solver, which utilized a fully MPI-parallel geometric multigrid
// with a GMRES(4)-AmaVanka smoother as a Richardson iteration. This application is configured to use a conforming
// Q2/P1dc discretization defined on a unstructured 2D quadrilateral mesh, however, both the Shape-type of the
// underlying mesh as well as the finite-element space can be easily changed by changing the corresponding typedefs in
// the code below. However, all the pre-defined problems are only defined in 2D, so you will have to come up with a
// different problem type if you want to change this application to a 3D code.
// The entire application functionality is contained in the main function below.
//
// The documentation of this application assumes that you are familiar with the basic concepts and classes of FEAT,
// which have been covered in the tutorials.
//
// To run this application, you have to specify at least the following three mandatory command line
// arguments:
//
// --problem <problem>
// This application defines a set of Poisson problems and this command line option specifies which of the predefined
// problems is to be solved; must be one of the following:
// * "dricav": solve Stokes on unit-square domain in lid-driven cavity configuration; no analytic solution
// * "taygre": solve Stokes on unit-square domain with the Taylor-Green-vortex analytic solution:
//             2D: u(x,y) = [ sin(pi*x)*cos(pi*y), -cos(pi*x)*sin(pi*y) ]
//                 p(x,y) = (cos(pi*x)^2 + cos(pi*y)^2 -1) / 2
//             Note: this solution has slip BCs on [0,1]^d, but it can also be solved with Dirichlet BCs
// * "parpro": solve Stokes on unit-square domain with the parabolic profile Poiseuille-Flow analytic solution:
//             2D: u(x,y) = [ 4*y*(1-y), 0 ]
//                 p(x,y) = 8*(1 - x)
// * "ribovo": solve Stokes on unit-circle domain with the "rigid-body-vortex" analytic solution:
//             2D: u(x,y) = [ -y, x ]
//                 p(x,y) = 0
//             Note: Dirichlet BCs on the entire boundary are required for this problem
//
// --mesh <meshfile>
// Specifies the mesh file to read. In principle, all of the above problems can be solved on any domain, as long as
// the mesh type matches this application's configuration, but naturally one would use a mesh discretizing the
// unit-circle domain for the 'ribovo' problem and the unit-square domain [0,1]^2 for all other of the above problems.
//
// --level <level-max> [[<levels-med...>] <level-min>]
// Specifies the refinement levels for the multigrid hierarchy. The mandatory first parameter <level-max> specifies the
// maximum refinement level on which the system is actually solved. The optional final parameter <level-min> specifies
// the minimum refinement level, which represents the coarse mesh level for the multigrid solver. The optional
// intermediary parameters <level-mid...> specify a sequence of partition change levels of the form 'level:ranks' where
// 'level' is the refinement level index and 'ranks' specifies the number of MPI ranks that should be used from this
// level to the next coarser partitioning change level or the coarse level, if no further partition change levels were
// given.
//
//
// In addition to the above mandatory parameters, this application also supports a small set of optional parameters:
//
// --slip
// Specifies that the application should assemble slip boundary conditions instead of no-slip boundary conditions,
// where ever applicable for the selected problem type.
//
// --deform
// Use deformation tensor Du:Dv instead of the gradient tensor grad(u):grad(v) as the bilinear operator.
//
// --vtk [<filename>]
// Specifies that the application should export the final solution to a (partitioned) VTU file. If no filename is
// given, then a default filename starting with 'poisson-simple-blocked' will be used.
//
// --no-err
// Specifies that the error computation of the final solution against the analytical reference solution is to be
// skipped.
//
// \author Peter Zajac
// ====================================================================================================================

#include <kernel/runtime.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/distance_function.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/multigrid.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>

// we're using the FEAT namespace here
using namespace FEAT;

int main(int argc, char* argv [])
{
  // always the very first step: create a runtime scope guard
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // create world communicator
  Dist::Comm comm(Dist::Comm::world());

  // print number of processes
  comm.print("Number of Processes: " + stringify(comm.size()));

  // create a stop watch
  StopWatch watch_total;
  watch_total.start();

  // create simple argument parser
  SimpleArgParser args(argc, argv);

  // add all supported command line arguments to parser
  Control::Domain::add_supported_pdc_args(args);
  args.support("mesh");
  args.support("level");
  args.support("deform");
  args.support("problem");
  args.support("no-err");
  args.support("slip");
  args.support("vtk");

  // check for unsupported options
  auto unsupported = args.query_unsupported();
  if (!unsupported.empty())
  {
    // print all unsupported options to cerr
    for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
      comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

    comm.print(std::cerr, "Supported Options are:");
    comm.print(std::cerr, args.get_supported_help());

    // abort
    FEAT::Runtime::abort();
  }

  // we always need a mesh file
  if(args.check("mesh") < 1)
  {
    comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
    FEAT::Runtime::abort();
  }

  // we always need refinement levels
  if(args.check("level") < 1)
  {
    comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
    FEAT::Runtime::abort();
  }

  // define our desired shape type and mesh type
  // 2D quadrilateral conforming mesh
  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;

  // Q2/P1dc finite element space pair defined on a standard first order transformation
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceTypeVelo;
  typedef Space::Discontinuous::ElementP1<TrafoType> SpaceTypePres;

  // get our shape-dimension
  static constexpr int shape_dim = ShapeType::dimension;

  // parse our problem type
  String problem;
  if(args.parse("problem", problem) < 1)
  {
    comm.print(std::cerr, "ERROR: Mandatory option '--problem <problem>' is missing!");
    FEAT::Runtime::abort();
  }
  if(!problem.is_one_of("taygre ribovo parpro dricav"))
  {
    comm.print(std::cerr, "ERROR: Invalid problem type: '" + problem + "'");
    FEAT::Runtime::abort();
  }
  if constexpr(shape_dim != 2)
  {
    if(problem.is_one_of("taygre ribovo parpro dricav"))
    {
      comm.print(std::cerr, "ERROR: Problem types 'taygre', 'ribovo', 'parpro' and 'dricav' are only available in 2D");
      FEAT::Runtime::abort();
    }
  }

  // use deformation tensor?
  const bool deformation_tensor = (args.check("deform") >= 0);

  StopWatch watch_domain_setup;
  watch_domain_setup.start();

  // create our domain control
  typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceTypeVelo, SpaceTypePres> DomainLevelType;
  typedef Control::Domain::PartiDomainControl<DomainLevelType> DomainControlType;

  // create the domain control, set the domain arguments and level and create the mesh hierarchy
  DomainControlType domain(comm, true);
  domain.parse_args(args);
  domain.set_desired_levels(args.query("level")->second);
  domain.create(args.query("mesh")->second);
  domain.add_trafo_mesh_part_charts();

  watch_domain_setup.stop();

  // plot our levels and partition info
  comm.print("Desired Levels: " + domain.format_desired_levels());
  comm.print("Chosen  Levels: " + domain.format_chosen_levels());
  comm.print("Partition Info: " + domain.get_chosen_parti_info());

  // get the number of levels on this process
  const Index num_levels = domain.size_physical();

  // ==================================================================================================================
  // System Assembly Phase
  // ==================================================================================================================

  // print selected problem
  if(problem == "ribovo")
    comm.print("Selected Problem: Rigid-Body-Vortex");
  else if(problem == "taygre")
    comm.print("Selected Problem: Taylor-Green-Vortex");
  else if(problem == "dricav")
    comm.print("Selected Problem: Lid-Driven Cavity");
  else if(problem == "parpro")
    comm.print("Selected Problem: Parabolic Profile Poiseuille-Flow");
  comm.print("Selected Tensor.: " + String(deformation_tensor ? "deformation tensor" : "gradient tensor"));

  comm.print("Assembling system...");

  // define our desired data and index types
  typedef double DataType;
  typedef Index IndexType;

  // define a set of analytic solutions

  // Taylor-Green vortex velocity and pressure
  Analytic::Common::TaylorGreenVortexVelo2D<DataType> func_taygre_velo;
  Analytic::Common::TaylorGreenVortexPres2D<DataType> func_taygre_pres;

  // rigid body vortex velocity and pressure
  Analytic::Common::RigidBodyVortexVelo2D<DataType> func_ribovo_velo;
  Analytic::Common::ConstantFunction<2, DataType> func_ribovo_pres(0.0);

  // parabolic profile velocity and pressure
  Analytic::Common::ParProfileVector<DataType> func_parpro_velo(0.0, 0.0, 0.0, 1.0, 1.0);
  auto func_parpro_pres = Analytic::create_lambda_function_scalar_2d(
    [](DataType x, DataType) {return DataType(8)*(1.0-x);},
    [](DataType, DataType) {return DataType(-8);},
    [](DataType, DataType) {return DataType(0);});

  // create lid-driven cavity boundary function; no analytical solution given for thsi problem
  auto func_dricav_velo = Analytic::create_lambda_function_vector_2d(
    [](DataType x, DataType y) {return ((y > 0.9999) && (x > 0.0001) && (x < 0.9999) ? 1.0 : 0.0); },
    [](DataType  , DataType  ) {return 0.0; });

  // define our system level
  typedef Control::StokesBlockedCombinedSystemLevel<shape_dim, DataType, IndexType> SystemLevelType;

  // create system levels deque
  std::deque<std::shared_ptr<SystemLevelType>> system_levels;
  for (Index i(0); i < num_levels; ++i)
    system_levels.push_back(std::make_shared<SystemLevelType>());

  StopWatch watch_system_setup;
  watch_system_setup.start();

  // define our cubature formula
  const String cubature("auto-degree:" + stringify(2*SpaceTypeVelo::local_degree+1));

  // assemble gate and compile domain assemblers
  for (Index i(0); i < num_levels; ++i)
  {
    domain.at(i)->domain_asm.compile_all_elements();
    system_levels.at(i)->assemble_gates(domain.at(i));
  }

  // assemble muxers and transfers
  for (Index i(0); (i < domain.size_physical()) && ((i+1) < domain.size_virtual()); ++i)
  {
    system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
    if((i+1) < domain.size_physical())
      system_levels.at(i)->assemble_transfers(*system_levels.at(i+1), domain.at(i), domain.at(i+1), cubature);
    else
      system_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature);
  }

  // assemble Laplace matrices with either gradient or deformation tensor
  Assembly::Common::LaplaceOperatorBlocked<shape_dim> laplace_operator; // for gradient tensor
  Assembly::Common::DuDvOperatorBlocked<shape_dim> dudv_operator;       // for deformation tensor

  for (Index i(0); i < num_levels; ++i)
  {
    // assemble matrix structures
    system_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
    system_levels.at(i)->assemble_pres_struct(domain.at(i)->space_pres);
    system_levels.at(i)->matrix_sys.local().format();

    // assemble diffusion matrix
    if(deformation_tensor)
      Assembly::assemble_bilinear_operator_matrix_1(domain.at(i)->domain_asm, system_levels.at(i)->matrix_a.local(),
        dudv_operator, domain.at(i)->space_velo, cubature);
    else
      Assembly::assemble_bilinear_operator_matrix_1(domain.at(i)->domain_asm, system_levels.at(i)->matrix_a.local(),
        laplace_operator, domain.at(i)->space_velo, cubature);

    // assemble gradient and divergence matrices
    system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->domain_asm, domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);

    // compile the system matrix
    system_levels.at(i)->compile_system_matrix();
  }

  // do we have to assemble slip instead of no-slip?
  const bool want_slip = (args.check("slip") >= 0);

  // assemble unit filter on all levels
  for (Index i(0); i < num_levels; ++i)
  {
    if(problem == "parpro")
    {
      // Parabolic inflow Profile
      system_levels.at(i)->assemble_velocity_unit_filter(*domain.at(i), domain.at(i)->space_velo, "inflow", "bnd:l", func_parpro_velo);
      if(want_slip)
      {
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:b", "bnd:b");
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:t", "bnd:t");
      }
      else
        system_levels.at(i)->assemble_velocity_unit_filter(*domain.at(i), domain.at(i)->space_velo, "noslip", "bnd:b bnd:t");
    }
    else if(problem == "ribovo")
    {
      // Rigid-Body-Vortex
      if(want_slip)
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "outer", "bnd:o");
      else
        system_levels.at(i)->assemble_velocity_unit_filter(*domain.at(i), domain.at(i)->space_velo, "outer", "bnd:o", func_ribovo_velo);
    }
    else if(problem == "dricav")
    {
      // Lid-Driven Cavity
      system_levels.at(i)->assemble_velocity_unit_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:t", "bnd:t", func_dricav_velo);
      if(want_slip)
      {
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:b", "bnd:b");
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:l", "bnd:l");
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:r", "bnd:r");
      }
      else
        system_levels.at(i)->assemble_velocity_unit_filter(*domain.at(i), domain.at(i)->space_velo, "noslip", "bnd:b bnd:l bnd:r");
    }
    else if(problem == "taygre")
    {
      // Taylor-Green vortex
      if(want_slip)
      {
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:b", "bnd:b");
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:t", "bnd:t");
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:l", "bnd:l");
        system_levels.at(i)->assemble_velocity_slip_filter(*domain.at(i), domain.at(i)->space_velo, "bnd:r", "bnd:r");
      }
      else
        system_levels.at(i)->assemble_velocity_unit_filter(*domain.at(i), domain.at(i)->space_velo, "dirichlet", "bnd:b bnd:t bnd:l bnd:r", func_taygre_velo);
    }

    // If there is no do-nothing (outflow) boundary, we'll need an pressure integral mean filter
    if(problem.is_one_of("taygre ribovo dricav"))
    {
      system_levels.at(i)->assemble_pressure_mean_filter(domain.at(i)->space_pres, "gauss-legendre:1");
    }

    // compile the system filters
    system_levels.at(i)->compile_system_filter();
  }

  // get our assembled vector types
  typedef typename SystemLevelType::LocalVeloVector LocalVeloVector;
  typedef typename SystemLevelType::LocalPresVector LocalPresVector;
  typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

  // fetch our finest levels
  DomainLevelType& the_domain_level = *domain.front();
  SystemLevelType& the_system_level = *system_levels.front();

  // create new vectors and format them
  GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
  GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

  vec_sol.format();
  vec_rhs.format();

  // assemble the right-hand-side force functional
  if(problem.is_one_of("ribovo dricav parpro"))
  {
    // rhs = 0
  }
  else if(problem == "taygre")
  {
    // assemble RHS for Taylor-Green
    const DataType pi = Math::pi<DataType>();
    auto force_taygre = Analytic::create_lambda_function_vector_2d(
      [pi](DataType x, DataType y) {return (-Math::cos(pi*x) + 2.0*pi*Math::cos(pi*y))*pi*Math::sin(pi*x);},
      [pi](DataType x, DataType y) {return (-Math::cos(pi*y) - 2.0*pi*Math::cos(pi*x))*pi*Math::sin(pi*y);});
    Assembly::assemble_force_function_vector(the_domain_level.domain_asm, vec_rhs.local().at<0>(), force_taygre, the_domain_level.space_velo, cubature);
  }

  // synchronize the vector
  vec_rhs.sync_0();

  // filter the vectors
  the_system_level.filter_sys.filter_sol(vec_sol);
  the_system_level.filter_sys.filter_rhs(vec_rhs);

  watch_system_setup.stop();

  // ==================================================================================================================
  // Multigrid Solver Setup Phase
  // ==================================================================================================================

  comm.print("Setting up solver...");

  StopWatch watch_solver_setup;
  watch_solver_setup.start();

  // Richardson ( MG-VCycle ( S: Richardson ( AmaVanka )  / C: Richardson ( AmaVanka )  )  )
  auto multigrid_hierarchy = std::make_shared<
    Solver::MultiGridHierarchy<
    typename SystemLevelType::GlobalSystemMatrix,
    typename SystemLevelType::GlobalSystemFilter,
    typename SystemLevelType::GlobalSystemTransfer
      > >(domain.size_virtual());

  // loop over all multigrid levels from the finest to the coarsest
  for (Index i(0); i < num_levels; ++i)
  {
    SystemLevelType& lvl = *system_levels.at(i);

    // compile local type-1 system matrix for our AmaVanka smoother
    lvl.compile_local_matrix_sys_type1();

    // create a damped AmaVanka smoother
    auto amavanka = Solver::new_amavanka(lvl.local_matrix_sys_type1, lvl.filter_sys.local());

    // skip singular elements; this may happen on the coarse level
    amavanka->set_skip_singular(true);
    auto schwarz  = Solver::new_schwarz_precond(amavanka, lvl.filter_sys);

    // stuff the AmaVanka-smoother into a FGMRES(4) solver
    auto smoother = Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, 4, 0.0, schwarz);
    smoother->set_min_iter(4);
    smoother->set_max_iter(4);

    if((i+1) < domain.size_virtual())
    {
      // this is an intermediate level
      multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
    }
    else
    {
      // this is the coarse grid level
      multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, smoother);
    }
  }

  // create a multigrid preconditioner with V-cycle
  auto mgv = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

  // create a CG solver preconditioned with the multigrid
  auto solver = Solver::new_richardson(the_system_level.matrix_sys, the_system_level.filter_sys, 1.0, mgv);

  // enable plotting
  solver->set_plot_mode(Solver::PlotMode::iter);

  // set tolerance and maximum iterations
  solver->set_tol_rel(1E-8);
  solver->set_max_iter(1000);

  // initialize
  multigrid_hierarchy->init();
  solver->init();

  watch_solver_setup.stop();

  // ==================================================================================================================
  // System Solution Phase
  // ==================================================================================================================

  comm.print("\nSolving...");

  StopWatch watch_solver_apply;
  watch_solver_apply.start();

  // solve
  auto result = Solver::solve(*solver, vec_sol, vec_rhs, the_system_level.matrix_sys, the_system_level.filter_sys);

  // check whether the solver was successful
  if (!Solver::status_success(result))
  {
    comm.print("Solver execution FAILED, with status: " + stringify(result));
  }

  // release solver
  solver->done();
  multigrid_hierarchy->done();

  watch_solver_apply.stop();

  // ==================================================================================================================
  // Post-Processing: Error Analysis Phase
  // ==================================================================================================================

  if((args.check("no-err") < 0) && problem.is_one_of("taygre parpro ribovo"))
  {
    // use a cubature formula of higher order for error computation
    const String cubature_error("auto-degree:" + stringify(2*SpaceTypeVelo::local_degree+2));

    // create errors result structures
    Assembly::DiscreteFunctionIntegral<LocalVeloVector, SpaceTypeVelo>::Type errors_v;
    Assembly::DiscreteFunctionIntegral<LocalPresVector, SpaceTypePres>::Type errors_p;

    if(problem == "taygre")
    {
      errors_v = Assembly::integrate_error_function<2>(the_domain_level.domain_asm, func_taygre_velo,
        vec_sol.local().at<0>(), the_domain_level.space_velo, cubature_error);
      errors_p = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_taygre_pres,
        vec_sol.local().at<1>(), the_domain_level.space_pres, cubature_error);
    }
    else if(problem == "parpro")
    {
      // Note that the parabolic profile velocity and pressure functions are contained in our trial space Q2/P1,so the
      // errors primarily depend on the solver stopping criterion and will in general increase with finer mesh levels!
      errors_v = Assembly::integrate_error_function<2>(the_domain_level.domain_asm, func_parpro_velo,
        vec_sol.local().at<0>(), the_domain_level.space_velo, cubature_error);
      errors_p = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_parpro_pres,
        vec_sol.local().at<1>(), the_domain_level.space_pres, cubature_error);
    }
    else if(problem == "ribovo")
    {
      // Note that the rigid body vortex velocity and pressure functions are contained in our trial space Q2/P1,so the
      // errors primarily depend on the solver stopping criterion and will in general increase with finer mesh levels!
      errors_v = Assembly::integrate_error_function<2>(the_domain_level.domain_asm, func_ribovo_velo,
        vec_sol.local().at<0>(), the_domain_level.space_velo, cubature_error);
      errors_p = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_ribovo_pres,
        vec_sol.local().at<1>(), the_domain_level.space_pres, cubature_error);
    }

    // synchronize errors over all processes
    errors_v.synchronize(comm);
    errors_p.synchronize(comm);

    // print error information
    comm.print("\nVelocity Error Analysis:\n" + errors_v.print_norms());
    comm.print("\nPressure Error Analysis:\n" + errors_p.print_norms());
  }

  // ==================================================================================================================
  // Post-Processing: VTK writing phase
  // ==================================================================================================================

  if (args.check("vtk") >= 0)
  {
    // build VTK name
    String vtk_name = String("./stokes-simple");
    vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
    vtk_name += "-n" + stringify(comm.size());
    args.parse("vtk", vtk_name);

    comm.print("\nWriting VTK file to'" + vtk_name + ".[p]vtu'...");

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

    // project velocity to the vertices
    LocalVeloVector vtx_velo;
    Assembly::DiscreteVertexProjector::project(vtx_velo, vec_sol.local().at<0>(), the_domain_level.space_velo);

    // project pressure to the cells
    LocalPresVector vtx_pres;
    Assembly::DiscreteCellProjector::project(vtx_pres, vec_sol.local().at<1>(), the_domain_level.space_pres);

    // write velocity
    exporter.add_vertex_vector("v", vtx_velo);
    exporter.add_cell_scalar("p", vtx_pres.elements());

    // finally, write the VTK file
    exporter.write(vtk_name, comm);
  }

  // print elapsed runtime and memory usage
  watch_total.stop();
  comm.print("\nTotal Runtime....: " + watch_total.elapsed_string().pad_front(10));
  comm.print(  "Domain Setup Time: " + watch_domain_setup.elapsed_string().pad_front(10) + watch_domain_setup.percent_of(watch_total));
  comm.print(  "System Setup Time: " + watch_system_setup.elapsed_string().pad_front(10) + watch_system_setup.percent_of(watch_total));
  comm.print(  "Solver Setup Time: " + watch_solver_setup.elapsed_string().pad_front(10) + watch_solver_setup.percent_of(watch_total));
  comm.print(  "Solver Apply Time: " + watch_solver_apply.elapsed_string().pad_front(10) + watch_solver_apply.percent_of(watch_total));
  comm.print(  "Peak Memory Usage: " + MemoryUsage::format_peak_physical_usage(comm));
}
