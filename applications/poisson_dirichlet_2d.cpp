#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/solver/basic_vcycle.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/precon_wrapper.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>

#include <control/domain/partitioner_domain_control.hpp>
#include <control/scalar_basic.hpp>

namespace PoissonDirichlet2D
{
  using namespace FEAST;

  template<typename Space_>
  class PoissonDirichletAssemblerLevel :
    public Control::ScalarBasicAssemblerLevel<Space_>
  {
  public:
    typedef Control::ScalarBasicAssemblerLevel<Space_> BaseClass;
    typedef Space_ SpaceType;
    typedef typename SpaceType::MeshType MeshType;

  public:
    explicit PoissonDirichletAssemblerLevel(typename BaseClass::DomainLevelType& dom_lvl) :
      BaseClass(dom_lvl)
    {
    }

    template<typename SystemLevel_>
    void assemble_system_matrix(SystemLevel_& sys_level)
    {
      // get the global matrix
      typename SystemLevel_::GlobalSystemMatrix& mat_glob = sys_level.matrix_sys;

      // get the local matrix
      typename SystemLevel_::LocalSystemMatrix& mat_loc = *mat_glob;

      // assemble matrix structure?
      if (mat_loc.empty())
      {
        Assembly::SymbolicMatrixAssembler<>::assemble1(mat_loc, this->space);
      }

      // assemble velocity laplace matrix
      {
        mat_loc.format();
        Assembly::Common::LaplaceOperator laplace_op;
        Assembly::BilinearOperatorAssembler::assemble_matrix1(mat_loc, laplace_op, this->space, this->cubature);
      }
    }

    template<typename SystemLevel_>
    void assemble_system_filter(SystemLevel_& sys_level)
    {
      // get our global system filter
      typename SystemLevel_::GlobalSystemFilter& fil_glob = sys_level.filter_sys;

      // get our local system filter
      typename SystemLevel_::LocalSystemFilter& fil_loc = *fil_glob;

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> unit_asm;

      std::deque<String> part_names = this->domain_level.get_mesh_node()->get_mesh_part_names();
      for(const auto& name : part_names)
      {
        if(name.starts_with('_'))
          continue;

        auto* mesh_part_node = this->domain_level.get_mesh_node()->find_mesh_part_node(name);

        // found it?
        if (mesh_part_node == nullptr)
          throw InternalError("Mesh Part Node 'boundary' not found!");

        // let's see if we have that mesh part
        // if it is nullptr, then our patch is not adjacent to that boundary part
        auto* mesh_part = mesh_part_node->get_mesh();
        if (mesh_part != nullptr)
        {

        // add to boundary assembler
        unit_asm.add_mesh_part(*mesh_part);
        }
      }

      // finally, assemble the filter
      unit_asm.assemble(fil_loc, this->space);
    }


    template<typename SystemLevel_, typename SolFunc_>
    typename SystemLevel_::GlobalSystemVector assemble_rhs_vector(SystemLevel_& sys_level, const SolFunc_& sol_func)
    {
      // create new vector
      typename SystemLevel_::GlobalSystemVector vec_rhs = sys_level.matrix_sys.create_vector_r();
      vec_rhs.format();

      // get the local vector
      typename SystemLevel_::LocalSystemVector& vec_f = *vec_rhs;

      // assemble the force
      Assembly::Common::LaplaceFunctional<SolFunc_> force_func(sol_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_f, force_func, this->space, this->cubature);

      // sync the vector
      vec_rhs.sync_0();

      // and filter it
      sys_level.filter_sys.filter_rhs(vec_rhs);
      return vec_rhs;
    }

    template<typename SystemLevel_>
    typename SystemLevel_::GlobalSystemVector assemble_sol_vector(SystemLevel_& sys_level)
    {
      typename SystemLevel_::GlobalSystemVector vec_sol = sys_level.matrix_sys.create_vector_r();
      vec_sol.format();
      sys_level.filter_sys.filter_sol(vec_sol);
      return vec_sol;
    }

    template<typename SystemLevel_, typename SolFunc_>
    void analyse_sol_vector(bool plot, SystemLevel_& sys_level, const typename SystemLevel_::GlobalSystemVector& vec_sol, const SolFunc_& sol_func)
    {
      typedef typename SystemLevel_::DataType DataType;

      // Compute and print the H0-/H1-errors
      Assembly::ScalarErrorInfo<DataType> errors = Assembly::ScalarErrorComputer<1>::compute(
        *vec_sol, sol_func, this->space, this->cubature);

      // synhronise all local errors
      errors.norm_h0 = sys_level.gate_sys.norm2(errors.norm_h0);
      errors.norm_h1 = sys_level.gate_sys.norm2(errors.norm_h1);

      // print errors
      if (plot)
      {
        std::cout << errors << std::endl;
      }
    }
  };

  template<typename MeshType_, typename TargetMatrixSolve_>
  void run(const int rank, const int nprocs, SimpleArgParser& args, Control::Domain::DomainControl<MeshType_>& domain)
  {
    // define our mesh type
    typedef MeshType_ MeshType;

    // define our arch types
    typedef typename Mem::Main MemType;
    typedef typename TargetMatrixSolve_::DataType DataType;
    typedef typename TargetMatrixSolve_::IndexType IndexType;

    // choose our desired analytical solution
    Analytic::Common::ExpBubbleFunction<2> sol_func;

    // define our domain type
    typedef Control::Domain::DomainControl<MeshType_> DomainControlType;

    // define our system level
    typedef Control::ScalarUnitFilterSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    // define our transfer level
    typedef Control::ScalarBasicTransferLevel<SystemLevelType> TransferLevelType;

    // define our trafo and FE spaces
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange1::Element<TrafoType> SpaceType;

    // define our assembler level
    typedef typename DomainControlType::LevelType DomainLevelType;
    typedef PoissonDirichletAssemblerLevel<SpaceType> AssemblerLevelType;

    // get our domain level and layer
    typedef typename DomainControlType::LayerType DomainLayerType;
    const DomainLayerType& layer = *domain.get_layers().back();
    const std::deque<DomainLevelType*>& domain_levels = domain.get_levels();

    std::deque<SystemLevelType*> system_levels;
    std::deque<AssemblerLevelType*> asm_levels;
    std::deque<TransferLevelType*> transfer_levels;

    const Index num_levels = Index(domain_levels.size());

    //Lin-Solve phase related typedefs
    //Main-CSR or CUDA-ELL
    typedef typename TargetMatrixSolve_::MemType MemTypeSolve;
    typedef Control::ScalarUnitFilterSystemLevel<MemTypeSolve, DataType, IndexType, TargetMatrixSolve_> SystemLevelTypeSolve;
    typedef Control::ScalarBasicTransferLevel<SystemLevelTypeSolve> TransferLevelTypeSolve;

    std::deque<SystemLevelTypeSolve*> system_levels_solve;
    std::deque<TransferLevelTypeSolve*> transfer_levels_solve;

    // create stokes and system levels
    for (Index i(0); i < num_levels; ++i)
    {
      asm_levels.push_back(new AssemblerLevelType(*domain_levels.at(i)));
      system_levels.push_back(new SystemLevelType());
      if (i > 0)
      {
        transfer_levels.push_back(new TransferLevelType(*system_levels.at(i - 1), *system_levels.at(i)));
      }
    }

    /* ***************************************************************************************** */

    if (rank == 0)
    {
      std::cout << "Creating gates..." << std::endl;
    }

    for (Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_gates(layer, *system_levels.at(i));
    }

    /* ***************************************************************************************** */

    if (rank == 0)
    {
      std::cout << "Assembling system matrices..." << std::endl;
    }

    for (Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_system_matrix(*system_levels.at(i));
    }

    /* ***************************************************************************************** */

    if (rank == 0)
    {
      std::cout << "Assembling system filters..." << std::endl;
    }

    for (Index i(0); i < num_levels; ++i)
    {
      asm_levels.at(i)->assemble_system_filter(*system_levels.at(i));
    }

    /* ***************************************************************************************** */

    if (rank == 0)
    {
      std::cout << "Assembling transfer matrices..." << std::endl;
    }

    for (Index i(0); (i + 1) < num_levels; ++i)
    {
      asm_levels.at(i + 1)->assemble_system_transfer(*transfer_levels.at(i), *asm_levels.at(i));
    }

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain_levels.back();
    SystemLevelType& the_system_level = *system_levels.back();
    AssemblerLevelType& the_asm_level = *asm_levels.back();

    // create our RHS and SOL vectors
    GlobalSystemVector vec_rhs = the_asm_level.assemble_rhs_vector(the_system_level, sol_func);
    GlobalSystemVector vec_sol = the_asm_level.assemble_sol_vector(the_system_level);

    ////////////////// solver type conversion ////////////////////////

    // get our global solver system types
    typedef typename SystemLevelTypeSolve::GlobalSystemVector GlobalSystemVectorSolve;
    typedef typename SystemLevelTypeSolve::GlobalSystemMatrix GlobalSystemMatrixSolve;
    typedef typename SystemLevelTypeSolve::GlobalSystemFilter GlobalSystemFilterSolve;

    if (rank == 0)
    {
      std::cout << "Converting assembled linear system from " + SystemLevelType::LocalScalarMatrix::name() <<
        ", Mem:" << MemType::name() << " to " << SystemLevelTypeSolve::LocalScalarMatrix::name() << ", Mem:" <<
        MemTypeSolve::name() << "..." << std::endl;
    }

    //convert system and transfer levels
    for (Index i(0); i < num_levels; ++i)
    {
      //system levels must be converted first, because transfer levels use their converted gates
      system_levels_solve.push_back(new SystemLevelTypeSolve());
      system_levels_solve.back()->convert(*system_levels.at(i));
      if (i > 0)
      {
        transfer_levels_solve.push_back(new TransferLevelTypeSolve());
        transfer_levels_solve.back()->convert(*system_levels_solve.at(i-1), *system_levels_solve.at(i), *transfer_levels.at(i-1));
      }
    }

    // get our global solve matrix and filter
    GlobalSystemMatrixSolve& matrix_solve = (*system_levels_solve.back()).matrix_sys;
    GlobalSystemFilterSolve& filter_solve = (*system_levels_solve.back()).filter_sys;

    //convert rhs and sol vectors
    GlobalSystemVectorSolve vec_rhs_solve;
    vec_rhs_solve.convert(&system_levels_solve.back()->gate_sys, vec_rhs);
    GlobalSystemVectorSolve vec_sol_solve;
    vec_sol_solve.convert(&system_levels_solve.back()->gate_sys, vec_sol);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a multigrid cycle
    auto mgv = std::make_shared<
      Solver::BasicVCycle<
      GlobalSystemMatrixSolve,
      GlobalSystemFilterSolve,
      typename TransferLevelTypeSolve::GlobalSystemTransferMatrix
      > >();

    // scaling factor
    DataType omega = DataType(0.2);

    // create coarse grid solver
    auto coarse_solver = Solver::new_scale_precond(system_levels_solve.front()->filter_sys, omega);
    mgv->set_coarse_level(system_levels_solve.front()->matrix_sys, system_levels_solve.front()->filter_sys, coarse_solver);

    // push levels into MGV
    auto jt = transfer_levels_solve.begin();
    for (auto it = ++system_levels_solve.begin(); it != system_levels_solve.end(); ++it, ++jt)
    {
      auto smoother = Solver::new_scale_precond((*it)->filter_sys, omega);
      mgv->push_level((*it)->matrix_sys, (*it)->filter_sys, (*jt)->prol_sys, (*jt)->rest_sys, smoother, smoother);
    }

    // create our solver
    //auto solver = Control::SolverFactory::create_scalar_solver(system_levels_solver, transfer_levels_solver, Runtime::global_property(), "linsolver");
    auto solver = Solver::new_pcg(matrix_solve, filter_solve, mgv);
    //auto solver = Solver::new_bicgstab(matrix_solve, filter_solve, mgv);
    //auto solver = Solver::new_fgmres(matrix_solve, filter_solve, 8, 0.0, mgv);
    //auto solver = Solver::new_richardson(matrix_solve, filter_solve, 1.0, mgv);

    // enable plotting
    solver->set_plot(rank == 0);

    // set tolerance
    solver->set_tol_rel(1E-8);
    solver->set_max_iter(1000);

    // initialise
    solver->init();

    Statistics::reset_flops();
    Statistics::reset_times();
    Statistics::reset_solver_statistics();

    TimeStamp at;

    // solve
    Solver::solve(*solver, vec_sol_solve, vec_rhs_solve, matrix_solve, filter_solve);
    TimeStamp bt;

    std::size_t la_size(0);
    std::for_each(system_levels.begin(), system_levels.end(), [&] (SystemLevelType * n) { la_size += n->bytes(); });
    std::for_each(transfer_levels.begin(), transfer_levels.end(), [&] (TransferLevelType * n) { la_size += n->bytes(); });
    std::size_t mpi_size(0);
    std::for_each(system_levels.begin(), system_levels.end(), [&] (SystemLevelType * n) { mpi_size += n->gate_sys.bytes(); });
    if (rank == 0 && args.check("statistics") >= 0)
    {
      std::cout<<std::endl<<solver->get_formated_solver_tree().trim()<<std::endl;
      String flops = Statistics::get_formated_flops(bt.elapsed(at), nprocs);
      std::cout<<"\nComplete solver TOE: "<<bt.elapsed(at)<<std::endl;
      std::cout<<flops<<std::endl<<std::endl;
      std::cout<<Statistics::get_formated_times(bt.elapsed(at))<<std::endl<<std::endl;
      std::cout<<"Domain size: " << double(domain.bytes())  / (1024. * 1024.)  << " MByte" << std::endl;
      std::cout<<"LA size: " << double(la_size) / (1024. * 1024.) << " MByte" << std::endl << std::endl;
      std::cout<<"MPI size: " << double(mpi_size) / (1024. * 1024.) << " MByte" << std::endl << std::endl;
      if (args.check("statistics") > 0) // provided parameter full or whatever
      {
        std::cout<<Statistics::get_formated_solvers();
      }
    }
    /// \todo add mpi related allocation size
    if (args.check("statistics") > 0) // provided parameter full or whatever
      Statistics::write_out_solver_statistics(rank, la_size, domain.bytes(), mpi_size);

    // release solver
    solver->done();

    // download solution
    vec_sol.convert(&system_levels.back()->gate_sys, vec_sol_solve);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("no-err") < 0)
    {
      the_asm_level.analyse_sol_vector(rank == 0, the_system_level, vec_sol, sol_func);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./poisson-dirichlet-2d");
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(nprocs);

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

      // project velocity and pressure
      typename SystemLevelType::LocalSystemVector vtx_sol, vtx_rhs;
      Assembly::DiscreteVertexProjector::project(vtx_sol, (*vec_sol), the_asm_level.space);
      Assembly::DiscreteVertexProjector::project(vtx_rhs, (*vec_rhs), the_asm_level.space);

      // write velocity
      exporter.add_scalar_vertex("sol", vtx_sol.elements());
      exporter.add_scalar_vertex("rhs", vtx_rhs.elements());

      // finally, write the VTK file
      exporter.write(vtk_name, rank, nprocs);
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // clean up
    while (!transfer_levels_solve.empty())
    {
      delete transfer_levels_solve.back();
      transfer_levels_solve.pop_back();
    }
    while (!system_levels_solve.empty())
    {
      delete system_levels_solve.back();
      system_levels_solve.pop_back();
    }

    while (!transfer_levels.empty())
    {
      delete transfer_levels.back();
      transfer_levels.pop_back();
    }
    while (!system_levels.empty())
    {
      delete system_levels.back();
      system_levels.pop_back();
    }
    while (!asm_levels.empty())
    {
      delete asm_levels.back();
      asm_levels.pop_back();
    }
  }

  int main(int argc, char* argv [])
  {
    int rank = 0;
    int nprocs = 0;

    // initialise
    FEAST::Runtime::initialise(argc, argv, rank, nprocs);
#ifndef SERIAL
    if (rank == 0)
    {
      std::cout << "NUM-PROCS: " << nprocs << std::endl;
    }
#endif

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    args.support("level");
    args.support("no-err");
    args.support("vtk");
    args.support("statistics");
    args.support("mem");
    args.support("part_min_elems");
    args.support("meshfile");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      if(rank == 0)
      {
        for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
          std::cerr << "ERROR: unknown option '--" << (*it).second << "'" << std::endl;
      }
      // abort
      FEAST::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;

    int lvl_max = 3;
    int lvl_min = 0;
    args.parse("level", lvl_max, lvl_min);

    FEAST::String mem_string = "main";
    args.parse("mem", mem_string);

#ifndef DEBUG
    try
#endif
    {
      TimeStamp stamp1;

      // let's create our domain
      if (rank == 0)
      {
        std::cout << "Preparing domain..." << std::endl;
      }
      int min_elems_partitioner(1);
      args.parse("part_min_elems", min_elems_partitioner);

      // fetch the mandatory mesh filename
      String meshfile;
      if(args.parse("meshfile", meshfile) <= 0)
      {
        if(rank == 0)
          std::cerr << "ERROR: Mandatory option --meshfile is missing!" << std::endl;
        FEAST::Runtime::abort();
      }
#ifdef FEAST_HAVE_PARMETIS
      Control::Domain::PartitionerDomainControl<Foundation::PExecutorParmetis<Foundation::ParmetisModePartKway>, MeshType> domain(lvl_max, lvl_min, Index(min_elems_partitioner), meshfile);
#elif !defined(SERIAL)
      Control::Domain::PartitionerDomainControl<Foundation::PExecutorFallback<double, Index>, MeshType> domain(lvl_max, lvl_min, Index(min_elems_partitioner), meshfile);
#else
      Control::Domain::PartitionerDomainControl<Foundation::PExecutorNONE<double, Index>, MeshType> domain(lvl_max, lvl_min, Index(min_elems_partitioner), meshfile);
#endif

      // plot our levels
      if (rank == 0)
      {
        std::cout << "LVL-MIN: " << domain.get_levels().front()->get_level_index() << " [" << lvl_min << "]" << std::endl;
        std::cout << "LVL-MAX: " << domain.get_levels().back()->get_level_index() << " [" << lvl_max << "]" << std::endl;
      }

      // run our application
      if (mem_string == "main")
      {
        run<MeshType, LAFEM::SparseMatrixCSR<Mem::Main, double, Index> >(rank, nprocs, args, domain);
      }
#ifdef FEAST_BACKENDS_CUDA
      else if(mem_string == "cuda")
      {
        run<MeshType, LAFEM::SparseMatrixELL<Mem::CUDA, double, Index> >(rank, nprocs, args, domain);
      }
#endif
      else
      {
        throw InternalError("Memory type " + mem_string + " not known!");
      }

      TimeStamp stamp2;

      // get times
      long long time1 = stamp2.elapsed_micros(stamp1);

      // accumulate times over all processes
      long long time2 = time1 * (long long) nprocs;

      // print time
      if (rank == 0)
      {
        std::cout << "Run-Time: "
          << TimeStamp::format_micros(time1, TimeFormat::m_s_m) << " ["
          << TimeStamp::format_micros(time2, TimeFormat::m_s_m) << "]" << std::endl;
      }
    }
#ifndef DEBUG
    catch (const std::exception& exc)
    {
      std::cerr << "ERROR: unhandled exception: " << exc.what() << std::endl;
      FEAST::Runtime::abort();
    }
    catch (...)
    {
      std::cerr << "ERROR: unknown exception" << std::endl;
      FEAST::Runtime::abort();
    }
#endif // DEBUG

    // okay
    return FEAST::Runtime::finalise();
  }
} // namespace PoissonDirichlet2D

int main(int argc, char* argv [])
{
  return PoissonDirichlet2D::main(argc, argv);
}
