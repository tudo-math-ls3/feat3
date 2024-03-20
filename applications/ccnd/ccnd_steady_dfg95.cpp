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
// This application supports all parameters from the CCND::SteadyAppBase class, so please refer to
// the documentation of the ccnd_steady_appbase.hpp header file for an up-to-date list of all
// command line parameters supported by that base class.
//
// In addition, this application adds support for the following parameters:
//
// --bench <1|7>
// Specifies which of the 2 steady benchmarks is to be solved:
//   --bench 1 corresponds to the test cases 2D-1 and 3D-1Z of flow-around-a-cylinder
//   --bench 7 corresponds to the 3D flow-around-a-sphere-inside-cylinder benchmark that is
//             analogous to the bench 1 simulation
//
// --v-max <vmax>
// Specifies the maximum velocity of the inflow boundary condition.
// Defaults to 0.3 in 2D and 0.45 in 3D.
//
// --fbm
// Enables the fictitious boundary method (FBM) to enforce the boundary conditions for the circular
// or cylindrical obstacle.
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
          system.at(i)->assemble_fbm_filters(*domain.at(i)->fbm_assembler, domain.at(i)->space_velo, domain.at(i)->space_pres, i == 0u);

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

      // compute pressure drop
      bench_analysis.compute_pressure_drop(comm, vec_sol.local().at<1>(), domain.front()->space_pres, bench);

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

    // load initial solution if desired
    if(app.load_initial_solution())
      app.newton_starts_with_picard = false; // skip Picard iteration in first Newton step

    // solve non-linear system
    app.solve_nonlinear_system();

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
