// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
//
// ------------------------------------------------------------------------------------------------
// Steady CCnD solver for 2D/3D Driven Cavity benchmark
// ------------------------------------------------------------------------------------------------
// This application implements a parallel steady CCND solver, which is pre-configured to solve
// the driven cavity benchmark on the unit-square or unit-cube domain, where the X-velocity
// at the top boundary edge/face (i.e. where y=1) is set to 1 in the unregularized case or to
// sin(pi*x)^gamma (2D) or (sin(pi*x)*sin(pi*z)^gamma (3D) and all other velocity components are
// set to 0.
//
// This application supports all parameters from the CCND::SteadyAppBase class, so please refer to
// the documentation of the ccnd_steady_appbase.hpp header file for an up-to-date list of all
// command line parameters supported by that base class.
//
// In addition, this application adds support for the following parameters:
// --gamma <gamma>
// Specifies the regularization parameter gamma for the boundary condition definition for the
// X-velocity component at the top edge of the domain. Defaults to 0.
// - If gamma = 0 ==> v_x := 1 at top edge, i.e. where y=1
// - If gamma > 0 in 2D ==> v_x := sin(pi*x)^gamma at top edge
// - If gamma > 0 in 3D ==> v_x := (sin(pi*x)*sin(pi*z))^gamma at top edge
//
// \author Peter Zajac
//

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
  class DrivenCavityBoundaryFunction :
    public Analytic::Function
  {
  public:
    typedef Analytic::Image::Vector<dim> ImageType;
    static constexpr int domain_dim = dim;
    static constexpr bool can_value = true;
    static constexpr bool can_grad = false;
    static constexpr bool can_hess = false;

    /// regularization parameter
    DataType _gamma;

    explicit DrivenCavityBoundaryFunction(DataType gamma) :
      _gamma(gamma)
    {
    }

    template<typename EvalTraits_>
    class Evaluator :
      public Analytic::Function::Evaluator<EvalTraits_>
    {
    public:
      typedef typename EvalTraits_::DataType DataType;
      typedef typename EvalTraits_::PointType PointType;
      typedef typename EvalTraits_::ValueType ValueType;

    private:
      const DataType _gamma;
      const DataType _pi;
      const bool _regularized;

    public:
      explicit Evaluator(const DrivenCavityBoundaryFunction& function) :
        _gamma(function._gamma),
        _pi(Math::pi<DataType>()),
        _regularized(_gamma > 1E-8)
      {
      }

      ValueType value(const PointType& point)
      {
        ValueType val;
        val.format();
        if(!_regularized)
          val[0] = DataType(1);
        else if constexpr(dim == 3)
          val[0] = Math::pow(Math::sin(_pi * point[0]) * Math::sin(_pi * point[2]), _gamma);
        else
          val[0] = Math::pow(Math::sin(_pi * point[0]), _gamma);
        return val;
      }
    };
  }; // class DrivenCavityBoundaryFunction

  /// our actual Application class
  class Application :
    public SteadyAppBase
  {
  public:
    /// our base class
    typedef SteadyAppBase BaseClass;

    /// boundary condition regularization parameter gamma
    DataType bc_gamma = DataType(0);

    explicit Application(const Dist::Comm& comm_, SimpleArgParser& args_) :
      BaseClass(comm_, args_)
    {
      args.support("gamma");

      // this application always has homogeneous RHS
      homogeneous_rhs = true;
    }

    virtual void parse_args() override
    {
      args.parse("gamma", bc_gamma);
      if(bc_gamma < DataType(0))
      {
        comm.print("ERROR: invalid bc regularization parameter gamma: " + stringify(bc_gamma));
        Runtime::abort();
      }

      // set default filename
      default_filename = String("cc") + stringify(dim) + "d-steady-dricav";

      // parse remaining arguments
      BaseClass::parse_args();
    }

    virtual void print_problem() override
    {
      BaseClass::print_problem();
      comm.print("\nDriven Cavity Benchmark Parameters:");
      comm.print(String("Regularization Gamma").pad_back(pad_len, pad_char) + ": " + stringify(bc_gamma));
    }

    virtual void assemble_filters()
    {
      // the names of the mesh parts on which to assemble
      std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

      /*auto top_func = Analytic::create_lambda_function_vector_2d(
        [](DataType x, DataType ) {return 2.0*x*(1.0-x);},
        [](DataType, DataType) {return 0.0;}
      );*/

      DrivenCavityBoundaryFunction dcbc_function(this->bc_gamma);

      for(std::size_t i(0); i < domain.size_physical(); ++i)
      {
        Assembly::UnitFilterAssembler<MeshType> unit_asm_top, unit_asm_noflow;

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
          if(mesh_part)
          {
            if(name == "bnd:t")
              unit_asm_top.add_mesh_part(*mesh_part);
            else
              unit_asm_noflow.add_mesh_part(*mesh_part);
          }
        }

        unit_asm_top.assemble(system.at(i)->get_local_velo_unit_filter_seq().find_or_add("top"), domain.at(i)->space_velo, dcbc_function);
        unit_asm_noflow.assemble(system.at(i)->get_local_velo_unit_filter_seq().find_or_add("noflow"), domain.at(i)->space_velo);

        // assemble pressure mean filter
        system.at(i)->assemble_pressure_mean_filter(domain.at(i)->space_pres, false);

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
      Assembly::DiscreteFunctionIntegralJob<LocalVeloVector, SpaceVeloType, 1> velo_job(vec_sol.local().template at<0>(), domain.front()->space_velo, "gauss-legendre:3");
      Assembly::DiscreteFunctionIntegralJob<LocalPresVector, SpacePresType, 0> pres_job(vec_sol.local().template at<1>(), domain.front()->space_pres, "gauss-legendre:2");
      domain.front()->domain_asm.assemble(velo_job);
      domain.front()->domain_asm.assemble(pres_job);

      velo_job.result().synchronize(comm);
      pres_job.result().synchronize(comm);

      comm.print("Velocity Field Analysis:");
      comm.print(velo_job.result().print_norms(15, pad_len));
      comm.print(velo_job.result().print_field_info(15, pad_len));
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
