// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// ====================================================================================================================
// Simple Scalar Poisson Example Application with fully-grown MPI-parallel Geometric Multigrid Solver
// --------------------------------------------------------------------------------------------------------------------
// This application is a simple scalar Poisson solver, which utilized a fully MPI-parallel geometric multigrid as a
// preconditioner for a parallel PCG solver. This application is configured to use a conforming Q1 discretization
// defined on a unstructured 2D quadrilateral mesh However, both the Shape-type of the underlying mesh as well as
// the finite-element space can be easily changed by changing the corresponding typedefs in the code below.
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
// * "one": solve -Laplace(u) = 1 (no analytic solution given)
// * "sin": solve -Laplace(u) = f with analytic solution
//          2D: u(x,y)   = sin(pi*x)*sin(pi*y)
//          3D: u(x,y,z) = sin(pi*x)*sin(pi*y)*sin(pi*z)
//          Note: this solution has homogeneous Dirichlet BCs on [0,1]^d
// * "cos": solve -Laplace(u) = f with analytic solution
//          2D: u(x,y)   = cos(pi*x)*cos(pi*y)
//          3D: u(x,y,z) = cos(pi*x)*cos(pi*y)*cos(pi*z)
//          Note: this solution has homogeneous Neumann BCs on [0,1]^d
// * "exp": solve -Laplace(u) = f with analytic solution
//          2D: u(x,y)   = (exp(-(2*x-1)^2) - exp(-1)) / (1 - exp(-1)) * (exp(-(2*y-1)^2) - exp(-1)) / (1 - exp(-1))
//          3D: u(x,y,z) = u(x,y) * (exp(-(2*z-1)^2) - exp(-1)) / (1 - exp(-1))
//          Note: this solution has homogeneous Dirichlet BCs on [0,1]^d
// * "sad": solve -Laplace(u) = 0 with analytic solution (2D only)
//          2D: u(x,y)   = x^2 - y^2
//
// --mesh <meshfile>
// Specifies the mesh file to read. In principle, all of the above problems can be solved on any domain, as long as
// the mesh type matches this application's configuration, but naturally one would use a mesh discretizing the
// unit-square domain [0,1]^d for the 'sin', 'cos' and 'exp' problems and the unit-circle domain for the 'sad' problem.
// The 'one' problem can be used on any domain, since it does not have an analytic reference solution.
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
// --vtk [<filename>]
// Specifies that the application should export the final solution to a (partitioned) VTU file. If no filename is
// given, then a default filename starting with 'poisson-simple-scalar' will be used.
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
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/multigrid.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>

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
  args.support("problem");
  args.support("no-err");
  args.support("dirichlet");
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
  if(!problem.is_one_of("one sin cos exp sad"))
  {
    comm.print(std::cerr, "ERROR: Invalid problem type: '" + problem + "'");
    FEAT::Runtime::abort();
  }
  if constexpr(shape_dim != 2)
  {
    if(problem == "sad")
    {
      comm.print(std::cerr, "ERROR: Problem type 'sad' is only available in 2D");
      FEAT::Runtime::abort();
    }
  }

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

  comm.print("Assembling system...");

  // define our desired data and index types
  typedef double DataType;
  typedef Index IndexType;

  // define a set of analytic solutions
  Analytic::Common::SineBubbleFunction<shape_dim> func_sin;
  Analytic::Common::CosineWaveFunction<shape_dim> func_cos;
  Analytic::Common::ExpBubbleFunction<shape_dim> func_exp;
  Analytic::Common::ConstantFunction<shape_dim, DataType> func_null(DataType(0));
  Analytic::Common::ConstantFunction<shape_dim, DataType> func_one(DataType(1));

  // create saddle function (2D only): x^2 - y^2
  auto func_sad = Analytic::create_lambda_function_scalar_2d(
    [](DataType x, DataType y) {return x*x - y*y;},
    [](DataType x, DataType  ) {return  DataType(2)*x;},
    [](DataType  , DataType y) {return -DataType(2)*y;});

  // define our system level
  typedef Control::ScalarCombinedSystemLevel<DataType, IndexType> SystemLevelType;

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

  // assemble scalar Laplace matrices
  Assembly::Common::LaplaceOperator laplace_operator;
  for (Index i(0); i < num_levels; ++i)
  {
    // assemble matrix structure and then its actual contents
    Assembly::SymbolicAssembler::assemble_matrix_std1(system_levels.at(i)->matrix_sys.local(), domain.at(i)->space);
    system_levels.at(i)->matrix_sys.local().format();
    Assembly::assemble_bilinear_operator_matrix_1(domain.at(i)->domain_asm, system_levels.at(i)->matrix_sys.local(),
      laplace_operator, domain.at(i)->space, cubature);
  }

  // do we need to assemble Dirichlet boundary conditions?
  if(args.check("dirichlet") >= 0)
  {
    // get the Dirichlet mesh part names and join them into a single string
    String dirichet_bnd_names = stringify_join(args.query("dirichlet")->second, " ");

    comm.print("Assembling Dirichlet BCs on mesh-parts " + dirichet_bnd_names);

    // assemble unit filter on all levels
    for (Index i(0); i < num_levels; ++i)
    {
      if(problem == "one")
      {
        system_levels.at(i)->assemble_unit_filter(*domain.at(i), domain.at(i)->space, "dirichlet", dirichet_bnd_names);
      }
      else if(problem == "sin")
      {
        system_levels.at(i)->assemble_unit_filter(*domain.at(i), domain.at(i)->space, "dirichlet", dirichet_bnd_names, func_sin);
      }
      else if(problem == "cos")
      {
        system_levels.at(i)->assemble_unit_filter(*domain.at(i), domain.at(i)->space, "dirichlet", dirichet_bnd_names, func_cos);
      }
      else if(problem == "exp")
      {
        system_levels.at(i)->assemble_unit_filter(*domain.at(i), domain.at(i)->space, "dirichlet", dirichet_bnd_names, func_exp);
      }
      else if(problem == "sad")
      {
        if constexpr (shape_dim == 2)
          system_levels.at(i)->assemble_unit_filter(*domain.at(i), domain.at(i)->space, "dirichlet", dirichet_bnd_names, func_sad);
      }
    }
  }
  else // no Dirichlet ==> all Neumann boundaries
  {
    // assemble mean filter on all levels instead
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_mean_filter(domain.at(i)->space, cubature);
    }
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
  if(problem == "one")
  {
    // rhs = 1
    Assembly::assemble_force_function_vector(the_domain_level.domain_asm, vec_rhs.local(), func_one, the_domain_level.space, cubature);
  }
  else if(problem == "sin")
  {
    Assembly::Common::LaplaceFunctional<decltype(func_sin)> force_sin(func_sin);
    Assembly::assemble_linear_functional_vector(the_domain_level.domain_asm, vec_rhs.local(), force_sin, the_domain_level.space, cubature);
  }
  else if(problem == "cos")
  {
    Assembly::Common::LaplaceFunctional<decltype(func_cos)> force_cos(func_cos);
    Assembly::assemble_linear_functional_vector(the_domain_level.domain_asm, vec_rhs.local(), force_cos, the_domain_level.space, cubature);
  }
  else if(problem == "exp")
  {
    Assembly::Common::LaplaceFunctional<decltype(func_exp)> force_exp(func_exp);
    Assembly::assemble_linear_functional_vector(the_domain_level.domain_asm, vec_rhs.local(), force_exp, the_domain_level.space, cubature);
  }
  else if(problem == "sad")
  {
    // rhs = 0
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
    if(problem == "one")
    {
      // no analytic result available, compare against null function
      errors = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_null, vec_sol.local(), the_domain_level.space, cubature_error);
    }
    else if(problem == "sin")
    {
      errors = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_sin, vec_sol.local(), the_domain_level.space, cubature_error);
    }
    else if(problem == "cos")
    {
      errors = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_cos, vec_sol.local(), the_domain_level.space, cubature_error);
    }
    else if(problem == "exp")
    {
      errors = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_exp, vec_sol.local(), the_domain_level.space, cubature_error);
    }
    else if(problem == "sad")
    {
      if constexpr(shape_dim == 2) // 2D only
        errors = Assembly::integrate_error_function<1>(the_domain_level.domain_asm, func_sad, vec_sol.local(), the_domain_level.space, cubature_error);
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
    String vtk_name = String("./poisson-simple-scalar");
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
    exporter.add_vertex_scalar("sol", vtx_sol.elements());
    exporter.add_vertex_scalar("rhs", vtx_rhs.elements());

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
