// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Unsteady CCnD solver for DFG95 Flow-Around-A-Cylinder Benchmarks
// ------------------------------------------------------------------------------------------------
// This application implements a parallel unsteady CCND solver, which is pre-configured to solve
// the infamous unsteady "flow-around-a-cylinder" benchmark problems, which are defined in
//
//     M. Schaefer and S. Turek: Benchmark Computations of Laminar Flow Around a Cylinder
//
// This application supports all parameters from the CCND::UnsteadyAppBase class, so please refer to
// the documentation of the ccnd_steady_appbase.hpp header file for an up-to-date list of all
// command line parameters supported by that base class.
//
// In addition, this application adds support for the following parameters:
//
// --bench <2|3>
// Specifies which of the 3 unsteady benchmarks is to be solved:
//   --bench 2 corresponds to the test cases 2D-2 and 3D-2Z, which is defined by a steady inflow
//             boundary condition function and is usually solved on the time interval [0, 30]
//   --bench 3 corresponds to the test cases 2D-3 and 3D-3Z, which is defined by an unsteady inflow
//             boundary condition function and is solved on the time interval [0, 8]
//
// --v-max <vmax>
// Specifies the maximum velocity of the inflow boundary condition.
// Defaults to 1.5 in 2D and 2.25 in 3D.
//
// --fbm
// Enables the fictitious boundary method (FBM) to enforce the boundary conditions for the circular
// or cylindrical obstacle.
//
// Happy vortex shedding! (^_^)
//
// \author Peter Zajac
//

// include common CCND headers
#include "ccnd_common.hpp"
#include "ccnd_common_dfg95.hpp"

// open the namespace and define the DomainLevel and SystemLevel classes
namespace CCND
{
  /// our domain level
  typedef CCND::DomainLevelBase DomainLevel;

  /// define our system level
  typedef CCND::SystemLevelBase SystemLevel;

} //  namespace CCND

  // include unsteady application base header
#include "ccnd_unsteady_appbase.hpp"

namespace CCND
{
  /// our actual Application class
  class Application :
    public UnsteadyAppBase
  {
  public:
    /// our base class
    typedef UnsteadyAppBase BaseClass;

    /// which benchmark to solve? 2 or 3?
    int bench = 2;

    /// maximum inflow velocity: 1.5 for 2D bench 2; 2.25 for 3D bench 2
    DataType v_max = DataType(dim) * DataType(0.75);

    /// DFG95 benchmark summary
    DFG95::BenchmarkAnalysis bench_analysis;

    explicit Application(const Dist::Comm& comm_, SimpleArgParser& args_) :
      BaseClass(comm_, args_)
    {
      args.support("bench");
      args.support("v-max");
      args.support("fbm");

      // this application has homogeneous RHS
      homogeneous_unsteady_rhs = true;
    }

    virtual void parse_args() override
    {
      args.parse("bench", bench);
      if((bench != 2) && (bench != 3))
      {
        comm.print("ERROR: invalid benchmark specified: " + stringify(bench));
        Runtime::abort();
      }

      // in case of bench 3, set the time interval to 8 and delta-t to 1/400
      if(bench == 3)
      {
        time_max = DataType(8);
        max_time_steps = 3200;
      }

      // enable FBM if desired
      enable_fbm = (args.check("fbm") >= 0);

      args.parse("v-max", v_max);

      // set default filename
      default_filename = String("cc") + stringify(dim) + "d-unsteady-dfg95-bench" + stringify(bench);
      if(enable_fbm)
        default_filename += "-fbm";

      // pre-define default parameters
      nu = DataType(0.001);

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

      // the names of the mesh parts on which to assemble
      std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

      for(std::size_t i(0); i < domain.size_physical(); ++i)
      {
        // get our local velocity filters
        auto& filter_v_noflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("noflow");
        auto& filter_v_inflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("inflow");

        // create unit-filter assembler
        Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow, unit_asm_noflow;

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
              unit_asm_inflow.add_mesh_part(*mesh_part);
            }
            else if((name != "bnd:r") && (name != "bnd:out"))
            {
              // outflow
              unit_asm_noflow.add_mesh_part(*mesh_part);
            }
          }
        }

        // assemble the filters
        unit_asm_inflow.assemble(filter_v_inflow, domain.at(i)->space_velo, inflow_func);
        unit_asm_noflow.assemble(filter_v_noflow, domain.at(i)->space_velo);

        // assemble FBM filters?
        if(enable_fbm)
          system.at(i)->assemble_fbm_filters(*domain.at(i)->fbm_assembler, domain.at(i)->space_velo, domain.at(i)->space_pres, i == 0u);

        // compile system filter
        system.at(i)->compile_system_filter();
      }
    }

    virtual void update_filters()
    {
      // reassemble inflow BCs for bench 3
      if(bench == 3)
      {
        DFG95::UnsteadyInflowFunction<dim> unsteady_inflow_func(v_max, cur_time);

        for(std::size_t i(0); i < domain.size_physical(); ++i)
        {
          // get our local velocity filter
          auto& filter_v_inflow = system.at(i)->get_local_velo_unit_filter_seq().find_or_add("inflow");

          // create unit-filter assembler
          Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow;

          // try to fetch the inflow mesh part node
          auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node("bnd:l");
          XASSERT(mesh_part_node != nullptr);

          auto* mesh_part = mesh_part_node->get_mesh();
          if (mesh_part != nullptr)
          {
            unit_asm_inflow.add_mesh_part(*mesh_part);
            unit_asm_inflow.assemble(filter_v_inflow, domain.at(i)->space_velo, unsteady_inflow_func);
          }

          // finally, compile the system filter
          system.at(i)->compile_system_filter();
        }
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
      // only perform the analysis in main steps
      if((main_step == 0u) || (cur_step  % main_step > 0u))
        return;

      watch_sol_analysis.start();

      // compute drag & lift by line integration
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

      // print the analysis
      comm.print(bench_analysis.format_compact(line_prefix() + " > "));

      watch_sol_analysis.stop();
      stats.times[Times::sol_analysis] = watch_sol_analysis.elapsed();
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
    app.create_system(true);

    // assemble filters
    app.assemble_filters();

    // compile systems
    app.compile_system_matrices();

    // create vectors
    app.create_vectors();

    // collect system statistics
    app.collect_system_statistics();

    // create solver
    app.create_solver();

    // initialize solver
    app.init_solver_symbolic();

    // restart from checkpoint?
    if(app.load_checkpoint())
    {
      // ok, application restarted from checkpoint
    }
    else if(app.load_initial_solution())
    {
      // ok, initial solution read in from file
    }
    else if(app.bench == 2)
    {
      // solve stokes to obtain initial guess
      app.solve_stokes();
    }

    comm.print("\nStarting Time-Loop...\n");

    // write initial VTK if desired
    app.write_vtk();

    // time step loop
    while(app.next_time_step())
    {
      // update filters for this time-step
      app.update_filters();

      // initialize solution for this time-step
      app.initialize_step_solution();

      // compute RHS for time step
      app.assemble_step_rhs();

      // solve Navier-Stokes for this step
      if(!app.solve_time_step())
      {
        break;
      }

      // perform solution analysis
      app.perform_solution_analysis();

      // analyze time derivative
      app.analyze_time_derivative();

      // write checkpoint
      app.save_checkpoint();

      // write VTK files if desired
      app.write_vtk();
    } // end of time step loop

    // collect statistics
    app.collect_solver_statistics();

    // release solver
    app.done_solver_symbolic();

    // write solutions if desired
    app.save_final_solution();

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
