// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Steady CCnD solver for 2D Taylor-Green Vortex benchmark
// ------------------------------------------------------------------------------------------------
//
// This application solves the steady-state 2D Taylor-Green benchmark and compares the discrete
// solution to the analytic solution by computing the H0/H1/H2 errors. The analytical solution
// for the velocity and pressure is given by
//
//  v(x,y) = [ sin(pi*x)*cos(pi*y) ; -cos(pi*x)*sin(pi*y) ]
//  p(x,y) = (cos(pi*x)^2 + cos(pi*y)^2 - 1) / 2
//
// The Taylor-Green benchmark is typically solved on the 2D unit square domain [0,1]x[0,1], but
// it can be practically defined on any other 2D domain, as long as one prescribes Dirichlet
// boundary conditions. On the unit square domain the velocity field also fulfills slip boundary
// conditions, so one may alternatively choose to enforce slip boundary conditions instead of the
// more restrictive Dirichlet boundary conditions by specifying the --slip parameter.
//
// This application supports all parameters from the CCND::SteadyAppBase class, so please refer to
// the documentation of the ccnd_steady_appbase.hpp header file for an up-to-date list of all
// command line parameters supported by that base class.
//
// In addition, this application adds support for the following parameters:
//
// --slip
// Enforce Slip boundary conditions instead of Dirichlet boundary conditions. This only yields the
// correct results on rectangular domains where the X/Y-coordinates of the boundary edges are
// integers, such as e.g. the unit square domain.
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

  /// our partidomaincontrol
  template<typename DomainLevel_>
  using PartitionControl = Control::Domain::PartiDomainControl<DomainLevel_>;

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

    /// slip filter BCs ?
    bool slip_bc = false;

    explicit Application(const Dist::Comm& comm_, SimpleArgParser& args_) :
      BaseClass(comm_, args_)
    {
      args.support("slip");

      // this application always has inhomogeneous RHS
      homogeneous_rhs = false;
    }

    virtual void parse_args() override
    {
      // use slip BCs ?
      slip_bc = (args.check("slip") >= 0);

      // set default filename
      default_filename = String("cc") + stringify(dim) + "d-steady-taygre";

      // parse remaining arguments
      BaseClass::parse_args();
    }

    virtual void print_problem() override
    {
      BaseClass::print_problem();
      comm.print("\nTaylor-Green Parameters:");
      comm.print(String("Boundary Conditions").pad_back(pad_len, pad_char) + ": " + (slip_bc ? "Slip" : "Dirichlet"));
    }

    virtual void assemble_filters()
    {
      // the names of the mesh parts on which to assemble
      std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

      Analytic::Common::TaylorGreenVortexVelo2D<DataType> taylor_green_velo;

      for(std::size_t i(0); i < domain.size_physical(); ++i)
      {
        // loop over all boundary parts except for the right one, which is outflow
        for(const auto& name : part_names)
        {
          // skip non-boundary mesh-parts
          if(!name.starts_with("bnd:"))
            continue;

          // try to fetch the corresponding mesh part node
          const auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
          XASSERT(mesh_part_node != nullptr);
          const MeshPartType* mesh_part = mesh_part_node->get_mesh();

          if(slip_bc)
          {
            Assembly::SlipFilterAssembler<TrafoType> slip_asm(domain.at(i)->trafo);
            if(mesh_part)
              slip_asm.add_mesh_part(*mesh_part);
            slip_asm.assemble(system.at(i)->get_local_velo_slip_filter_seq().find_or_add(name), domain.at(i)->space_velo);
          }
          else
          {
            Assembly::UnitFilterAssembler<MeshType> unit_asm;
            if(mesh_part)
              unit_asm.add_mesh_part(*mesh_part);
            unit_asm.assemble(system.at(i)->get_local_velo_unit_filter_seq().find_or_add(name), domain.at(i)->space_velo, taylor_green_velo);
          }
        }

        // assemble pressure mean filter
        system.at(i)->assemble_pressure_mean_filter(domain.at(i)->space_pres, false);

        // synchronize slip filters
        if(slip_bc)
          system.at(i)->sync_velocity_slip_filters();

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

      // assemble Taylor-Green RHS: the RHS 'f' of the Taylor-Green vortex for the steady-state
      // Navier-Stokes equations is equal to 2*nu*pi^2*v, where 'v' is the velocity solution
      Analytic::Common::TaylorGreenVortexVelo2D<DataType> taylor_green_velo;
      Assembly::assemble_force_function_vector(domain.front()->domain_asm, vec_rhs.local().first(), taylor_green_velo,
        domain.front()->space_velo, "gauss-legendre:3", DataType(2) * this->nu * Math::sqr(Math::pi<DataType>()));

      // if we only solve Stokes rather than Navier-Stokes, then we have to subtract the convective
      // term from the RHS, which corresponds (by pure coincidence, of course) to -grad(P):
      if(!navier)
      {
        Analytic::Common::TaylorGreenVortexPres2D<DataType> taylor_green_pres;
        Analytic::Gradient<decltype(taylor_green_pres)> grad_p(taylor_green_pres);
        Assembly::assemble_force_function_vector(domain.front()->domain_asm, vec_rhs.local().first(), grad_p,
          domain.front()->space_velo, "gauss-legendre:3", DataType(1));
      }

      // synchronize
      vec_rhs.sync_0();

      // apply filters
      system.front()->filter_sys.filter_sol(vec_sol);
      system.front()->filter_sys.filter_rhs(vec_rhs);
    }

    virtual void perform_solution_analysis()
    {
      watch_sol_analysis.start();

      Analytic::Common::TaylorGreenVortexVelo2D<DataType> taylor_green_velo;
      Analytic::Common::TaylorGreenVortexPres2D<DataType> taylor_green_pres;

      // compute velocity error
      auto error_v = Assembly::integrate_error_function<2>(domain.front()->domain_asm, taylor_green_velo,
        vec_sol.local().at<0>(), domain.front()->space_velo, "gauss-legendre:5");
      auto error_p = Assembly::integrate_error_function<1>(domain.front()->domain_asm, taylor_green_pres,
        vec_sol.local().at<1>(), domain.front()->space_pres, "gauss-legendre:4");

      error_v.synchronize(comm);
      error_p.synchronize(comm);

      comm.print("Velocity Errors:\n" + error_v.print_norms(15, pad_len) + "\n");
      comm.print("Pressure Errors:\n" + error_p.print_norms(15, pad_len) + "\n");

      // write all relevant errors in a single line for easy parsing
      comm.print(String("Errors V:H0/H1/H2 P:H0/H1: ")
        + stringify_fp_sci(Math::sqrt(error_v.norm_h0_sqr)) + " "
        + stringify_fp_sci(Math::sqrt(error_v.norm_h1_sqr)) + " "
        + stringify_fp_sci(Math::sqrt(error_v.norm_h2_sqr)) + " "
        + stringify_fp_sci(Math::sqrt(error_p.norm_h0_sqr)) + " "
        + stringify_fp_sci(Math::sqrt(error_p.norm_h1_sqr)) + "\n");

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
    std::cerr << "ERROR: " << e.what() << "\n";
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception\n";
    FEAT::Runtime::abort();
  }
  return 0;
}
