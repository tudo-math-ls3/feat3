// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Steady CCnD solver for DFG95 Flow-Around-A-Cylinder Benchmarks
// ------------------------------------------------------------------------------------------------
// This application implements a parallel steady CCND solver, which is pre-configured to solve
// the infamous unsteady "flow-around-a-cylinder" benchmark problem, which is defined in
//
//     M. Schaefer and S. Turek: Benchmark Computations of Laminar Flow Around a Cylinder
//
// The system is discretized using an isoparametric Q2/P1dc finite element discretization.
// The monolithic nonlinear Oseen systems are solved using an adaptive Newton-Multigrid solver
// with an additive matrix-based Vanka smoother ("AmaVanka") and using UMFPACK (if available)
// as a coarse grid solver. This application supports recursive partitioning.
//
//
// ------------------------------------
// Basic Setup and Mandatory Parameters
// ------------------------------------
// This application defines default values for most of its parameters, however, three parameters
// are mandatory and always have to be specified explicitly:
//
// --bench <1|7>
// Specifies which of the 2 steady benchmarks is to be solved:
//   --bench 1 corresponds to the test cases 2D-1 and 3D-1Z of flow-around-a-cylinder
//   --bench 7 corresponds to the 3D flow-around-a-sphere-inside-cylinder benchmark that is
//             analogous to the bench 1 simulation
//
// --mesh <meshfiles...>
// Specifies the input mesh file(s).
//
// --level <level-max> [levels...] <level-min>
// Specifies the mesh refinement levels in the syntax according to Control::PartiDomainControl.
//
//
// ----------------------------
// System Definition Parameters
// ----------------------------
// Some of the basic parameters of the underlying system can be configured by the following
// parameters.
//
// --nu <nu>
// Specifies the viscosity parameter nu for the diffusion operator. Defaults to 1E-3.
//
// --defo
// If specified, the deformation tensor formulation is used for the diffusion operator,
// otherwise the gradient tensor formulation is used. Defaults to gradient tensor.
//
// --v-max <vmax>
// Specifies the maximum velocity of the inflow boundary condition.
// Defaults to 1.5 in 2D and 2.25 in 3D.
//
// --upsam <ups>
// Specifies the stabilization parameter <ups> for the streamline diffusion stabilization.
// Defaults to 0, i.e. unstabilized.
//
// --stokes
// If specified, only the steady-state Stokes equations are solved and the nonlinear iteration
// for solving the Navier-Stokes equations is skipped entirely.
//
//
// -------------------------------
// Solver Configuration Parameters
// -------------------------------
// This section describes the parameters that control the non-linear Newton/Picard solver
// as well as its multigrid preconditioner and its smoother component.
//
// --picard
// If specified, the nonlinear system in each time step will be solved using a simple
// Picard iteration instead of the Newton iteration.
//
// --plot-mg-iter
// If specified, the convergence plot of the multigrid solver in each nonlinear solver iteration
// is printed.
//
// --min-nl-iter <N>
// Specifies the minimum number of nonlinear (Newton/Picard) solver iterations per time step.
// Defaults to 1.
//
// --max-nl-iter <N>
// Specifies the maximum number of nonlinear (Newton/Picard) solver iterations per time step.
// Defaults to 10.
//
// --min-mg-iter <N>
// Specifies the minimum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 1.
//
// --max-mg-iter <N>
// Specifies the maximum number of multigrid iterations per nonlinear solver iteration.
// Defaults to 5.
//
// --smooth-steps <N>
// Specifies the number of pre- and post-smoothing AmaVanka steps. Defaults to 8.
//
// --smooth-damp <omega>
// Specifies the damping parameter for the AmaVanka smoother. Defaults to 0.7.
//
// --no-umfpack
// If specified, the multigrid solver will use a BiCGStab-AmaVanka solver as the coarse grid solver
// instead of the UMFPACK direct solver. Note that UMFPACK is only used if it is included in the
// build id and if the coarse system is solved on a single process.
//
// --nl-tol-abs <tol>
// Specifies the absolute tolerance for the nonlinear solver. Defaults to 1E-8.
//
// --mg-tol-rel <tol>
// If given, specifies the relative tolerance for the multigrid solver.
// If not given, then the tolerance for the multigrid solver is chosen in an adaptive
// manner depending on the two previous nonlinear solver defects, which is the default case.
// The adaptive tolerance is chosen in each nonlinear iteration by analyzing the nonlinear
// defect improvement in the previous nonlinear (Newton/Picard) solver iteration in the
// following manner: Let def_{j} and def_{j-1} denote the two previous nonlinear defect norms,
// then the next nonlinear defect norm def_{j+1} should approximately fulfill the equation
//
//        (def_{j+1} / def_{j}) \approx (def_{j} / def_{j+1})^C
//
// where C \in {1,2} is the convergence speed of the nonlinear solver, i.e. C=2 for Newton
// and C=1 for Picard iteration. Multiplying the above equation by def_{j} gives us an
// estimate for the next nonlinear defect norm def_{j+1}:
//
//        def_{j+1} \approx def_{j} * (def_{j} / def_{j+1})^C
//
// To obtain an absolute tolerance for the next multigrid solver application, we simply
// multiply the above estimate by 0.1.
//
//
// -----------------------------------------------------
// Initial Solution Read-In and Final Solution Write-Out
// -----------------------------------------------------
//
// --save-sol <filename>
// Specifies that the application should write the final (partitioned) solution (and the
// partitioning) to a single binary output file. The output file can be loaded by the --load-sol
// option if the input mesh, the refinement level as well as the partitioning (and thus the
// process count) are identical.
//
// --save-joined-sol <filename>
// Specifies that the application should write the final solution into a joined binary output
// file by utilizing the base splitter. The output file can be loaded by the --load-joined-sol
// option if the input mesh and refinement level are identical, however, the process count
// and/or partitioning may differ. This feature should only be used with at most one parallel
// domain layer and moderate process counts.
//
// --load-sol <filename>
// Specifies that the application should read in the initial (partitioned) solution guess
// from a single binary output file, which was written by a --save-sol from a previous run.
// If specified, the solving of the Stokes system to obtain an initial guess is skipped.
//
// --load-joined-sol <filename> [<scale>]
// Specifies that the application should read in the initial joined solution guess from a
// single binary output file, which was written by a --save-joined-sol from a previous run.
// The second option argument <scale> is given, then the loaded solution is scaled by that factor,
// which can be used if the loaded solution was computed with a different inflow velocity.
// If specified, the solving of the Stokes system to obtain an initial guess is skipped.
//
// ------------------------
// Miscellaneous Parameters
// ------------------------
// This section describes miscellaneous parameters that do not fit into any other section and
// which do not deserve a custom section of their own.
//
// --vtk <filename> [<refined-filename>]
// Specifies that the application should write a VTK visualization output file. The second
// optional parameter specifies the filename for a VTK output file on a once refined mesh,
// if given. Note: it is possible to output only the refined VTKs by passing a whitespace
// string as the first filename, e.g.: --vtk " " myvtk
//
// --ext-stats
// If given, specifies that the application should output extensive statistics at the end of the
// program run, including detailed MPI timings.
//
// --test-mode
// If given, specifies that the application should run in test-mode rather than its normal mode.
// In test mode, the application only perform 3 time steps and writes some additional output which
// is parsed by the test system.
//
//
// \author Peter Zajac
//

// include common CCND header
#include "ccnd_common.hpp"
#include "ccnd_common_dfg95.hpp"

// open the namespace and define the DomainLevel and SystemLevel classes
namespace CCND
{
  /// our domain level
  typedef DomainLevelBase DomainLevel;

  /// our system level
  typedef SystemLevelBase SystemLevel;

} //  namespace CCND

// include steady application base header
#include "ccnd_steady_appbase.hpp"

namespace CCND
{
  /// our actual Application class
  class Application :
    public SteadyAppBase
  {
  public:
    /// our base class
    typedef SteadyAppBase BaseClass;

    /// which benchmark to solve? 1 or 7?
    int bench = 1;

    /// maximum inflow velocity: 0.3 for 2D bench 1, 0.45 for 3D bench 1; 1 for bench 7
    DataType v_max = DataType(1);

    /// DFG95 benchmark summary
    DFG95::BenchmarkAnalysis bench_analysis;

    explicit Application(const Dist::Comm& comm_, SimpleArgParser& args_) :
      BaseClass(comm_, args_)
    {
      args.support("bench");
      args.support("v-max");
      args.support("fbm");

      // this application always has homogeneous RHS
      homogeneous_rhs = true;
    }

    virtual void parse_args() override
    {
      args.parse("bench", bench);
      if((bench != 1) && (bench != 7))
      {
        comm.print("ERROR: invalid benchmark specified: " + stringify(bench));
        Runtime::abort();
      }
      if((bench == 7) && (dim != 3))
      {
        comm.print("ERROR: bench 7 is only valid for 3D");
        Runtime::abort();
      }

      // enable FBM if desired
      enable_fbm = (args.check("fbm") >= 0);

      // pre-define default parameters
      nu = DataType(0.001);
      if(bench == 1)
        v_max = DataType(dim) * DataType(0.15);

      args.parse("v-max", v_max);

      // set default filename
      default_filename = String("cc") + stringify(dim) + "d-steady-dfg95-bench" + stringify(bench);
      if(enable_fbm)
        default_filename += "-fbm";

      // parse remaining arguments
      BaseClass::parse_args();
    }

    virtual void print_problem() override
    {
      BaseClass::print_problem();
      comm.print("\nDFG95 Benchmark Parameters:");
      comm.print(String("Benchmark").pad_back(pad_len, pad_char) + ": bench " + stringify(bench));
      comm.print(String("V-Max").pad_back(pad_len, pad_char) + ": " + stringify(v_max));
    }

    virtual void create_domain() override
    {
      BaseClass::create_domain();

      // create FBM meshpart?
      if(enable_fbm)
      {
        // obstacle midpoint: 2D: (0.2,0.2)
        Tiny::Vector<DataType, dim> mp(DataType(0.2));
        if constexpr (dim == 3)
          mp[0] = DataType(0.5);

        Geometry::SphereHitTestFunction<DataType, dim> fbm_hit_func(mp, 0.05);

        // create meshparts on all levels
        for(std::size_t i(0); i < domain.size_physical(); ++i)
        {
          auto* mesh_node = domain.at(i)->get_mesh_node();
          Geometry::HitTestFactory<decltype(fbm_hit_func), MeshType> hit_factory(fbm_hit_func, *mesh_node->get_mesh());
          mesh_node->add_mesh_part("fbm", hit_factory.make_unique());

          // create FBM assembler
          domain.at(i)->create_fbm_assembler(domain.at(i).layer().comm(), "fbm");
        }
      }

      // add all isoparametric mesh-parts
      BaseClass::add_all_isoparam_parts();
    }

    virtual void create_benchmark_analysis()
    {
      // let the analysis get its mesh parts and charts
      bench_analysis.create(domain.get_atlas(), *domain.front()->get_mesh_node(), domain.front()->space_velo);
    }

    virtual void assemble_filters()
    {
      // our inflow BC function
      DFG95::SteadyInflowFunction<dim> inflow_func(v_max);

      Tiny::Vector<DataType, dim> pipe_origin, pipe_axis;
      pipe_origin = DataType(0.205);
      pipe_axis = DataType(0);
      pipe_axis[0] = DataType(1);
      Analytic::Common::PoiseuillePipeFlow<DataType, dim> pipe_inflow_func(pipe_origin, pipe_axis, DataType(0.205), v_max);

      // the names of the mesh parts on which to assemble
      std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

      for(std::size_t i(0); i < domain.size_physical(); ++i)
      {
        // get our local velocity filters
        auto& filter_v_noflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("noflow");
        auto& filter_v_inflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("inflow");

        // create unit-filter assembler
        Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow_1, unit_asm_inflow_2, unit_asm_noflow;

        // loop over all boundary parts except for the right one, which is outflow
        for(const auto& name : part_names)
        {
          // skip non-boundary mesh-parts
          if(!name.starts_with("bnd:"))
            continue;

          // try to fetch the corresponding mesh part node
          auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
          XASSERT(mesh_part_node != nullptr);

          auto* mesh_part = mesh_part_node->get_mesh();
          if (mesh_part != nullptr)
          {
            if(name == "bnd:l")
            {
              // inflow
              unit_asm_inflow_1.add_mesh_part(*mesh_part);
            }
            else if(name == "bnd:in")
            {
              // inflow
              unit_asm_inflow_2.add_mesh_part(*mesh_part);
            }
            else if((name != "bnd:r") && (name != "bnd:out"))
            {
              // outflow
              unit_asm_noflow.add_mesh_part(*mesh_part);
            }
          }
        }

        // assemble the filters
        unit_asm_inflow_1.assemble(filter_v_inflow, domain.at(i)->space_velo, inflow_func);
        unit_asm_inflow_2.assemble(filter_v_inflow, domain.at(i)->space_velo, pipe_inflow_func);
        unit_asm_noflow.assemble(filter_v_noflow, domain.at(i)->space_velo);

        // assemble FBM filters?
        if(enable_fbm)
        {
          auto& fbm_asm = *domain.at(i)->fbm_assembler;
          auto& filter_fbm_p = system.at(i)->get_local_pres_unit_filter();
          auto& filter_fbm_v = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("fbm");
          auto& filter_fbm_int_v = system.at(i)->filter_interface_fbm;

          // assemble velocity unit filter
          fbm_asm.assemble_inside_filter(filter_fbm_v, domain.at(i)->space_velo);
          fbm_asm.assemble_inside_filter(filter_fbm_p, domain.at(i)->space_pres);

          // assemble interface filter
          fbm_asm.assemble_interface_filter(filter_fbm_int_v, domain.at(i)->space_velo, system.at(i)->matrix_a, system.at(i)->velo_mass_matrix);

          // assemble mask vectors on finest level
          if(i == 0u)
          {
            auto& mask_v = system.at(i)->fbm_mask_velo;
            mask_v.reserve(domain.at(i)->space_velo.get_num_dofs());
            for(int d(0); d <= dim; ++d)
            {
              for(auto k : fbm_asm.get_fbm_mask_vector(d))
                mask_v.push_back(k);
            }
            system.at(i)->fbm_mask_pres = fbm_asm.get_fbm_mask_vector(dim);
          }
        }

        // compile system filter
        system.at(i)->compile_system_filter();
      }
    }

    virtual void create_vectors() override
    {
      BaseClass::create_vectors();

      // format the vectors
      vec_sol.format();
      vec_rhs.format();

      // apply filters
      system.front()->filter_sys.filter_sol(vec_sol);
      system.front()->filter_sys.filter_rhs(vec_rhs);
    }

    virtual void perform_solution_analysis()
    {
      watch_sol_analysis.start();

      // compute drag & lift by line integration
      if(bench == 1)
        bench_analysis.compute_body_forces_line(comm, vec_sol.local().at<0>(), vec_sol.local().at<1>(),
          domain.front()->space_velo, domain.front()->space_pres, nu, v_max);

      // compute drag & lift coefficients via volume integration from unsynchronized final defect
      bench_analysis.compute_body_forces_vol(comm, vec_def_unsynced.local().first(), nu, v_max, bench);

      // compute pressure difference
      bench_analysis.compute_pressure_difference(comm, vec_sol.local().at<1>(), domain.front()->space_pres);

      // compute fluxes
      bench_analysis.compute_fluxes(comm, vec_sol.local().first(), domain.front()->space_velo);

      // perform analysis of velocity field
      bench_analysis.compute_velocity_info(comm, vec_sol.local().first(), domain.front()->space_velo);

      // print analysis
      comm.print(bench_analysis.format(bench == 1));

      watch_sol_analysis.stop();
    }
  }; // class Application


  void main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));
    comm.print("Floating Point Type: " + String(fp_typename) + " precision");

    // create arg parser
    SimpleArgParser args(argc, argv);

    // create application
    Application app(comm, args);

    // check for unsupported arguments
    app.check_args();

    // parse arguments
    app.parse_args();

    // read mesh and create domain
    app.create_domain();

    // print problem parameters
    app.print_problem();

    // create benchmark analysis stuff
    app.create_benchmark_analysis();

    // create system
    app.create_system();

    // assemble filters
    app.assemble_filters();

    // compile system systems
    app.compile_system_matrices();

    // create vectors
    app.create_vectors();

    // collect system statistics
    app.collect_system_statistics();

    // create solver
    app.create_solver();

    // initialize solver
    app.init_solver_symbolic();

    // load initial solution or solve Stokes
    if(!app.load_initial_solution())
      app.solve_stokes();

    // solve Navier-Stokes if desired
    if(app.navier)
      app.solve_navier_stokes();
    else // compute unsynchronized defect for body forces computation
      app.compute_unsynced_defect();

    // collect statistics
    app.collect_solver_statistics();

    // release solver
    app.done_solver_symbolic();

    // write solutions if desired
    app.save_final_solution();

    // write VTK files if desired
    app.write_vtk();

    // perform solution analysis
    app.perform_solution_analysis();

    // collect final statistics
    app.collect_final_statistics();

    // print statistics
    app.print_statistics();
  }
} // namespace CCND

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  try
  {
    CCND::main(argc, argv);
  }
  catch(std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return 0;
}
