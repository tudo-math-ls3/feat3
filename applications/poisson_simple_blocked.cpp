// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// ====================================================================================================================
// Simple Vector-Valued Poisson Example Application with fully-grown MPI-parallel Geometric Multigrid Solver
// --------------------------------------------------------------------------------------------------------------------
// This application is a simple vector-valued Poisson solver, which utilized a fully MPI-parallel geometric multigrid
// as a preconditioner for a parallel PCG solver. This application is configured to use a conforming Q1 discretization
// defined on a unstructured 2D quadrilateral mesh. However, both the Shape-type of the underlying mesh as well as
// the finite-element space can be easily changed by changing the corresponding typedefs in the code below, however,
// all the pre-defined problems are only defined in 2D, so you will have to come up with a different problem type if
// you want to change this application to a 3D code.
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
// * "ribovo": solve -Laplace(u) = f with the "rigid-body-vortex" analytic solution:
//             2D: u(x,y)   =  [ -y, x ]
//             Note: Dirichlet BCs on the entire boundary are required for this problem
// * "taygre": solve -Laplace(u) = f with the Taylor-Green-vortex analytic solution:
//             2D: u(x,y)   = [ sin(pi*x)*cos(pi*y), -cos(pi*x)*sin(pi*y) ]
//             Note: this solution has slip BCs on [0,1]^d
//
// --mesh <meshfile>
// Specifies the mesh file to read. In principle, all of the above problems can be solved on any domain, as long as
// the mesh type matches this application's configuration, but naturally one would use a mesh discretizing the
// unit-square domain [0,1]^d for the 'taygre' problem and the unit-circle domain for the 'ribovo' problem.
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
// --dirichlet [<meshparts...>]
// Specifies the names of the mesh-parts that have to be treated as Dirichlet boundary condition regions. If you only
// specify '--dirichlet' without supplying a list of mesh-part names, then all mesh-parts of the mesh file are treated
// as Dirichlet boundary regions.
// If you do not specify the '--dirichlet' parameter, then the entire boundary is treated as pure do-nothing boundary,
// which results in a homogeneous Neumann boundary region in the case of a simple Poisson equation.
//
// --slip [<meshparts...>]
// Specifies the names of the mesh-parts that have to be treated as Dirichlet boundary condition regions.
// The mesh-parts are parsed in the same way as the ones from the --dirichlet parameter.
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
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/distance_function.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/multigrid.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/blocked_basic.hpp>

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
  args.support("dirichlet");
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

  // Q1 finite element space defined on a standard first order transformation
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  // get our shape-dimension
  static constexpr int shape_dim = ShapeType::dimension;

  // parse our problem type
  String problem;
  if(args.parse("problem", problem) < 1)
  {
    comm.print(std::cerr, "ERROR: Mandatory option '--problem <problem>' is missing!");
    FEAT::Runtime::abort();
  }
  if(!problem.is_one_of("taygre ribovo"))
  {
    comm.print(std::cerr, "ERROR: Invalid problem type: '" + problem + "'");
    FEAT::Runtime::abort();
  }
  if constexpr(shape_dim != 2)
  {
    if(problem.is_one_of("taygre ribovo"))
    {
      comm.print(std::cerr, "ERROR: Problem types 'taygre' and 'ribovo' are only available in 2D");
      FEAT::Runtime::abort();
    }
  }

  // use deformation tensor?
  const bool deformation_tensor = (args.check("deform") >= 0);

  StopWatch watch_domain_setup;
  watch_domain_setup.start();

  // create our domain control
  typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
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
  comm.print("Selected Tensor.: " + String(deformation_tensor ? "deformation tensor" : "gradient tensor"));

  comm.print("Assembling system...");

  // define our desired data and index types
  typedef double DataType;
  typedef Index IndexType;

  // define a set of analytic solutions
  Analytic::Common::TaylorGreenVortexVelo2D<DataType> func_taygre;
  Analytic::Common::RigidBodyVortexVelo2D<DataType> func_ribovo;

  // define our system level
  typedef Control::BlockedCombinedSystemLevel<shape_dim, DataType, IndexType> SystemLevelType;

  // create system levels deque
  std::deque<std::shared_ptr<SystemLevelType>> system_levels;
  for (Index i(0); i < num_levels; ++i)
    system_levels.push_back(std::make_shared<SystemLevelType>());

  StopWatch watch_system_setup;
  watch_system_setup.start();

  // define our cubature formula
  const String cubature("auto-degree:" + stringify(2*SpaceType::local_degree+1));

  // assemble gate and compile domain assemblers
  for (Index i(0); i < num_levels; ++i)
  {
    domain.at(i)->domain_asm.compile_all_elements();
    system_levels.at(i)->assemble_gate(domain.at(i));
  }

  // assemble muxers and transfers
  for (Index i(0); (i < domain.size_physical()) && ((i+1) < domain.size_virtual()); ++i)
  {
    system_levels.at(i)->assemble_coarse_muxer(domain.at(i+1));
    if((i+1) < domain.size_physical())
      system_levels.at(i)->assemble_transfer(*system_levels.at(i+1), domain.at(i), domain.at(i+1), cubature);
    else
      system_levels.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);
  }

  // assemble Laplace matrices with either gradient or deformation tensor
  Assembly::Common::LaplaceOperatorBlocked<shape_dim> laplace_operator; // for gradient tensor
  Assembly::Common::DuDvOperatorBlocked<shape_dim> dudv_operator; // for deformation tensor

  for (Index i(0); i < num_levels; ++i)
  {
    // assemble matrix structure and then its actual contents
    Assembly::SymbolicAssembler::assemble_matrix_std1(system_levels.at(i)->matrix_sys.local(), domain.at(i)->space);
    system_levels.at(i)->matrix_sys.local().format();
      Assembly::assemble_bilinear_operator_matrix_1(domain.at(i)->domain_asm, system_levels.at(i)->matrix_sys.local(),
        laplace_operator, domain.at(i)->space, cubature);

    if(deformation_tensor)
      Assembly::assemble_bilinear_operator_matrix_1(domain.at(i)->domain_asm, system_levels.at(i)->matrix_sys.local(),
        dudv_operator, domain.at(i)->space, cubature);
    else
      Assembly::assemble_bilinear_operator_matrix_1(domain.at(i)->domain_asm, system_levels.at(i)->matrix_sys.local(),
        laplace_operator, domain.at(i)->space, cubature);
  }

  // do we need to assemble Dirichlet boundary conditions?
  const bool have_diri = (args.check("dirichlet") >= 0);
  const bool have_slip = (args.check("slip") >= 0);
  if(have_diri)
  {
    // get the Dirichlet mesh part names
    const std::deque<String> bnd_names = args.query("dirichlet")->second;
    if(!bnd_names.empty())
      comm.print("Assembling Dirichlet BCs on mesh-parts '" + stringify_join(bnd_names, "', '") + "'...");
    else
      comm.print("Assembling Dirichlet BCs on all available mesh-parts...");

    // assemble unit filter on all levels
    for (Index i(0); i < num_levels; ++i)
    {
      if constexpr(shape_dim == 2)
      {
        if(problem == "ribovo")
          system_levels.at(i)->assemble_unit_filter(*domain.at(i), domain.at(i)->space, "dirichlet", stringify_join(bnd_names, " "), func_ribovo);
        else if(problem == "taygre")
          system_levels.at(i)->assemble_unit_filter(*domain.at(i), domain.at(i)->space, "dirichlet", stringify_join(bnd_names, " "), func_taygre);
      }
    }
  }
  if(have_slip)
  {
    // get the Dirichlet mesh part names
    const std::deque<String> bnd_names = args.query("slip")->second;
    if(!bnd_names.empty())
      comm.print("Assembling Slip BCs on mesh-parts '" + stringify_join(bnd_names, "', '") + "'...");
    else
      comm.print("Assembling Slip BCs on all available mesh-parts...");

    // assemble unit filter on all levels
    for (Index i(0); i < num_levels; ++i)
    {
      // we need a separate slip filter for each boundary
      for(const String& name : bnd_names)
        system_levels.at(i)->assemble_slip_filter(*domain.at(i), domain.at(i)->space, name, name);
    }
  }
  // neither Dirichlet nor Slip BCs ==> do-nothing/Neumann BCs
  if(!have_diri && !have_slip)
  {
    // all Neumann boundary: we would need a mean filter, but this isn't implemented yet
    comm.print("WARNING: No Dirichlet or Slip boundary conditions given; resulting system will be singular!");
  }

  // get our assembled vector types
  typedef typename SystemLevelType::LocalSystemVector LocalSystemVector;
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
  if(problem == "ribovo")
  {
    // rhs = 0
  }
  else if(problem == "taygre")
  {
    Assembly::Common::LaplaceFunctional<decltype(func_taygre)> force_taygre(func_taygre);
    Assembly::assemble_linear_functional_vector(the_domain_level.domain_asm, vec_rhs.local(), force_taygre, the_domain_level.space, cubature);
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

  // PCG ( MG-VCycle ( S: Richardson ( Jacobi )  / C: Richardson ( Jacobi )  )  )
  auto multigrid_hierarchy = std::make_shared<
    Solver::MultiGridHierarchy<
    typename SystemLevelType::GlobalSystemMatrix,
    typename SystemLevelType::GlobalSystemFilter,
    typename SystemLevelType::GlobalSystemTransfer
      > >(domain.size_virtual());

  // loop over all multigrid levels from the finest to the coarsest
  for (Index i(0); i < num_levels; ++i)
  {
    const SystemLevelType& lvl = *system_levels.at(i);

    // create a damped Jacobi smoother
    auto jacobi = Solver::new_jacobi_precond(lvl.matrix_sys, lvl.filter_sys, 0.7);
    auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, 1.0, jacobi);
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
  auto solver = Solver::new_pcg(the_system_level.matrix_sys, the_system_level.filter_sys, mgv);

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

  if (args.check("no-err") < 0)
  {
    // use a cubature formula of higher order for error computation
    const String cubature_error("auto-degree:" + stringify(2*SpaceType::local_degree+2));

    // create errors result structure
    Assembly::DiscreteFunctionIntegral<LocalSystemVector, SpaceType>::Type errors;

    // compute errors against reference solution
    if constexpr(shape_dim == 2)
    {
      if(problem == "ribovo")
        errors = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_ribovo, vec_sol.local(), the_domain_level.space, cubature_error);
      else if(problem == "taygre")
        errors = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_taygre, vec_sol.local(), the_domain_level.space, cubature_error);
    }

    // synchronize errors over all processes
    errors.synchronize(comm);

    // print error information
    comm.print("\nError Analysis:\n" + errors.print_norms());
  }

  // ==================================================================================================================
  // Post-Processing: VTK writing phase
  // ==================================================================================================================

  if (args.check("vtk") >= 0)
  {
    // build VTK name
    String vtk_name = String("./poisson-simple-blocked");
    vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
    vtk_name += "-n" + stringify(comm.size());
    args.parse("vtk", vtk_name);

    comm.print("\nWriting VTK file to'" + vtk_name + ".[p]vtu'...");

    // Create a VTK exporter for our mesh
    Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

    // project solution and rhs vectors to the vertices
    typename SystemLevelType::LocalSystemVector vtx_sol, vtx_rhs;
    Assembly::DiscreteVertexProjector::project(vtx_sol, vec_sol.local(), the_domain_level.space);
    Assembly::DiscreteVertexProjector::project(vtx_rhs, vec_rhs.local(), the_domain_level.space);

    // write velocity
    exporter.add_vertex_vector("sol", vtx_sol);
    exporter.add_vertex_vector("rhs", vtx_rhs);

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
