// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Steady CCnD solver for Y-Pipe Benchmark
// ------------------------------------------------------------------------------------------------
// This application implements a parallel steady CCND solver, which is pre-configured to solve
//
// --v-max <vmax>
// Specifies the maximum velocity of the inflow boundary condition. Defaults to 1.
//
// --fbm
// Enables the fictitious boundary method (FBM) to enforce the boundary conditions for the circular
// or cylindrical obstacle.
//
// \author Peter Zajac
//

#ifdef FEAT_CCND_APP_DIM
#if FEAT_CCND_APP_DIM != 2
#error FEAT_CCND_APP_DIM must be equal to 2 for this application
#endif
#else
#define FEAT_CCND_APP_DIM 2
#endif

// include common CCND header
#include "ccnd_common.hpp"

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

    /// maximum inflow velocity: 0.3 for 2D bench 1, 0.45 for 3D bench 1; 1 for bench 7
    DataType v_max = DataType(1);

    /// trace assemblers for inflow and outflows
    std::unique_ptr<Assembly::TraceAssembler<TrafoType>> flux_asm_in, flux_asm_out_b, flux_asm_out_t;


    explicit Application(const Dist::Comm& comm_, SimpleArgParser& args_) :
      BaseClass(comm_, args_)
    {
      args.support("v-max");
      args.support("fbm");

      // this application always has homogeneous RHS
      homogeneous_rhs = true;
    }

    virtual void parse_args() override
    {
      // enable FBM if desired
      enable_fbm = (args.check("fbm") >= 0);

      // v-max given?
      args.parse("v-max", v_max);

      // set default filename
      default_filename = String("cc") + stringify(dim) + "d-steady-ypipe";
      if(enable_fbm)
        default_filename += "-fbm";

      // parse remaining arguments
      BaseClass::parse_args();
    }

    virtual void print_problem() override
    {
      BaseClass::print_problem();
      comm.print("\nY-Pipe Benchmark Parameters:");
      comm.print(String("V-Max").pad_back(pad_len, pad_char) + ": " + stringify(v_max));
    }

    virtual void create_domain() override
    {
      BaseClass::create_domain();

      // add all isoparametric mesh-parts
      BaseClass::add_all_isoparam_parts();

      // create FBM meshpart?
      if(enable_fbm)
      {
        ChartBaseType* chart = domain.get_atlas().find_mesh_chart("ypipe:fbm");
        if(chart == nullptr)
        {
          comm.print("ERROR: chart 'ypipe:fbm' not found in mesh!");
          Runtime::abort();
        }

        // create meshparts on all levels
        for(std::size_t i(0); i < domain.size_physical(); ++i)
        {
          // create mesh-part from chart by using a hit-test factory
          Geometry::ChartHitTestFactory<MeshType> factory(domain.at(i)->get_mesh(), *chart);
          domain.at(i)->get_mesh_node()->add_mesh_part("fbm", factory.make_unique());

          // create FBM assembler
          domain.at(i)->create_fbm_assembler(domain.at(i).layer().comm(), "fbm");
        }
      }

      // create trace assemblers on finest level
      {
        flux_asm_in.reset(new Assembly::TraceAssembler<TrafoType>(domain.front()->trafo));
        flux_asm_out_t.reset(new Assembly::TraceAssembler<TrafoType>(domain.front()->trafo));
        flux_asm_out_b.reset(new Assembly::TraceAssembler<TrafoType>(domain.front()->trafo));

        const MeshNodeType& mesh_node = *domain.front()->get_mesh_node();
        const MeshPartType* mesh_part_bnd_l = mesh_node.find_mesh_part("bnd:l");
        const MeshPartType* mesh_part_bnd_in = mesh_node.find_mesh_part("bnd:in");
        const MeshPartType* mesh_part_bnd_out_t = mesh_node.find_mesh_part("bnd:out:t");
        const MeshPartType* mesh_part_bnd_out_b = mesh_node.find_mesh_part("bnd:out:b");

        if(mesh_part_bnd_l)
          flux_asm_in->add_mesh_part(*mesh_part_bnd_l);
        if(mesh_part_bnd_in)
          flux_asm_in->add_mesh_part(*mesh_part_bnd_in);
        if(mesh_part_bnd_out_t)
          flux_asm_out_t->add_mesh_part(*mesh_part_bnd_out_t);
        if(mesh_part_bnd_out_b)
          flux_asm_out_b->add_mesh_part(*mesh_part_bnd_out_b);

        flux_asm_in->compile();
        flux_asm_out_t->compile();
        flux_asm_out_b->compile();
      }
    }

    virtual void assemble_filters()
    {
      // our inflow BC function
      Analytic::Common::ParProfileVector<DataType> inflow_func(0.0, -0.45, 0.0, 0.45, v_max);

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
            if((name == "bnd:l") || (name == "bnd:in"))
            {
              // inflow
              unit_asm_inflow.add_mesh_part(*mesh_part);
            }
            else if((name != "bnd:r") && !name.starts_with("bnd:out"))
            {
              // not outflow
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

      Cubature::DynamicFactory cubature_factory("gauss-legendre:2");

      const LocalVeloVector& vec_sol_v = this->vec_sol.local().first();
      const SpaceVeloType& space_v = domain.front()->space_velo;

      DataType fx[3] =
      {
        flux_asm_in->assemble_discrete_integral(vec_sol_v, space_v, cubature_factory)[0],
        flux_asm_out_t->assemble_discrete_integral(vec_sol_v, space_v, cubature_factory)[0],
        flux_asm_out_b->assemble_discrete_integral(vec_sol_v, space_v, cubature_factory)[0]
      };

      comm.allreduce(fx, fx, 3u, Dist::op_sum);

      comm.print("Solution Analysis:");
      comm.print("Flux In...: " + stringify_fp_fix(fx[0], fp_num_digs));
      comm.print("Flux Out T: " + stringify_fp_fix(fx[1], fp_num_digs) + " [" + stringify_fp_fix(100.0*fx[1]/fx[0],3,7) + "%]");
      comm.print("Flux Out B: " + stringify_fp_fix(fx[2], fp_num_digs) + " [" + stringify_fp_fix(100.0*fx[2]/fx[0],3,7) + "%]");
      comm.print("Mass Loss.: " + stringify_fp_fix(100.0*(1.0 - (fx[1]+fx[2])/fx[0]), 5) + "%\n");

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

    // collect statistics
    app.collect_solver_statistics();

    // release solver
    app.done_solver_symbolic();

    // write solutions if desired
    app.save_final_solution();

    // write VTK files if desired
    app.write_vtk();

    // analyze solution
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
