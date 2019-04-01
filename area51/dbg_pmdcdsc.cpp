// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/global/symmetric_lumped_schur_matrix.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/scale_precond.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>

#include <kernel/global/pmdcdsc_matrix.hpp>
#include <kernel/solver/hypre.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/domain/unit_cube_domain_control.hpp>
#include <control/scalar_mixed.hpp>
#include <control/statistics.hpp>

namespace PoissonMixed
{
  using namespace FEAT;

  template<typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    int dbg_rank = -1;
    args.parse("debug", dbg_rank);

    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our arch types
    typedef Mem::Main MemType;
    typedef double DataType;
    typedef Index IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh and shape types
    typedef typename DomainLevelType::SpacePresType SpacePresType;
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainControlType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // This has homogeneous Neumann BCs
    //Analytic::Common::CosineWaveFunction<dim> sol_func;
    //Analytic::Common::SineBubbleFunction<dim> sol_func;
    Analytic::Common::ExpBubbleFunction<dim> sol_func;

    // define our system level
    typedef Control::ScalarMixedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    Cubature::DynamicFactory cubature("auto-degree:7");

    // get domain level
    const DomainLevelType& the_domain_level = *domain.front();

    // create system level
    SystemLevelType the_system_level;

#ifdef FEAT_COMPILER_MICROSOFT
    if(comm.rank() == dbg_rank) __debugbreak();
#endif

    // assemble gates
    the_system_level.assemble_gates(domain.front());

    // assemble matrices
    the_system_level.assemble_velo_struct(the_domain_level.space_velo);
    the_system_level.assemble_grad_div_matrices(the_domain_level.space_velo, the_domain_level.space_pres, cubature);

    // assemble velocity mass matrix
    {
      the_system_level.matrix_a.local().format();
      Assembly::Common::IdentityOperatorBlocked<dim> identity_op;
      Assembly::BilinearOperatorAssembler::assemble_block_matrix1(
        the_system_level.matrix_a.local(), identity_op, the_domain_level.space_velo, cubature);
    }

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalPresVector GlobalSystemVector;

    auto& matrix_a = the_system_level.matrix_a;
    auto& matrix_b = the_system_level.matrix_b;
    auto& matrix_d = the_system_level.matrix_d;

    // lump velocity matrix
    typename SystemLevelType::GlobalVeloVector lumped_matrix_a(matrix_a.create_vector_l());
    matrix_a.local().lump_rows(lumped_matrix_a.local());
    lumped_matrix_a.sync_0();

    auto inv_lumped_a = lumped_matrix_a.clone();
    inv_lumped_a.component_invert(lumped_matrix_a);

    // create new vector
    GlobalSystemVector vec_rhs = the_system_level.matrix_b.create_vector_r();
    vec_rhs.format();

    {
      // get the local vector
      typename SystemLevelType::LocalPresVector& vec_f = vec_rhs.local();

      // assemble the force
      Assembly::Common::LaplaceFunctional<decltype(sol_func)> force_func(sol_func);
      Assembly::LinearFunctionalAssembler::assemble_vector(vec_f, force_func, the_domain_level.space_pres, cubature);

      // sync the vector
      vec_rhs.sync_0();
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    typedef Global::SymmetricLumpedSchurMatrix
    <
      typename SystemLevelType::GlobalVeloVector,
      typename SystemLevelType::GlobalMatrixBlockB,
      typename SystemLevelType::GlobalMatrixBlockD,
      typename SystemLevelType::GlobalVeloFilter
    > GlobalLumpedSystemMatrix;

    // schur matrix
    GlobalLumpedSystemMatrix matrix_gls(lumped_matrix_a, matrix_b, matrix_d, the_system_level.filter_velo);

    /* ***************************************************************************************** */

    typedef Global::PMDCDSCMatrix<
      typename SystemLevelType::GlobalMatrixBlockB,
      typename SystemLevelType::GlobalMatrixBlockD
    > GlobalPMDCDSCMatrix;

    GlobalPMDCDSCMatrix matrix_pmdcdsc(inv_lumped_a, matrix_b, matrix_d);

    comm.print("\nInitialising PMDCDSCMatrix...");
    TimeStamp ts_1;
    matrix_pmdcdsc.init_symbolic();
    TimeStamp ts_2;
    matrix_pmdcdsc.init_numeric();
    TimeStamp ts_3;

    comm.print("Symbolic Init Time: " + ts_2.elapsed_string(ts_1, TimeFormat::s_m));
    comm.print("Numeric  Init Time: " + ts_3.elapsed_string(ts_2, TimeFormat::s_m));

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if(args.check("bench") >= 0)
    {
      comm.print("\nPerforming Mat-Vec-Mult Benchmark...");

      GlobalSystemVector vec_tmp1 = vec_rhs.clone();
      GlobalSystemVector vec_tmp2 = vec_rhs.clone();

      int num_mv = 100;
      args.parse("bench", num_mv);

      // perform mat-vec benchmark
      TimeStamp ts_mv1;
      for(int i(0); i < num_mv; ++i)
      {
        matrix_gls.apply(vec_tmp2, vec_tmp1);
      }
      double time_mv_1 = ts_mv1.elapsed_now() / double(num_mv);
      comm.print("Old Time per MV: " + stringify_fp_fix(time_mv_1, 5, 9));

      TimeStamp ts_mv2;
      for(int i(0); i < num_mv; ++i)
      {
        matrix_pmdcdsc.apply(vec_tmp2, vec_tmp1);
      }
      double time_mv_2 =  ts_mv2.elapsed_now() / double(num_mv);
      comm.print("New Time per MV: " + stringify_fp_fix(time_mv_2, 5, 9));

      comm.print("Speedup........: " + stringify_fp_fix(time_mv_1/time_mv_2, 3, 7));
    }


    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */


    double time_1 = 0.0, time_2 = 0.0;

    auto vec_sol_1 = vec_rhs.clone();
    vec_sol_1.format();

    if(args.check("pcg") >= 0)
    {
      comm.print("\nSolving System using SymmetricLumpedSchurMatrix...");
      //auto precond = Solver::new_jacobi_precond(matrix_gls, the_system_level.filter_pres);
      auto solver = Solver::new_pcg(matrix_gls, the_system_level.filter_pres/*, precond*/);

      // enable plotting
      solver->set_plot_mode(Solver::PlotMode::all);

      // set tolerance
      solver->set_tol_rel(1E-5);
      solver->set_max_iter(1000);

      // initialise
      solver->init();

      TimeStamp ts;
      solver->correct(vec_sol_1, vec_rhs);
      time_1 = ts.elapsed_now();

      // release solver
      solver->done();
    }

    auto vec_sol_2 = vec_rhs.clone();
    vec_sol_2.format();

    if(args.check("pcg") >= 0)
    {
      comm.print("\nSolving System using PMDCDSCMatrix...");
      //auto precond = Solver::new_jacobi_precond(matrix_pmdcdsc, the_system_level.filter_pres);
      auto solver = Solver::new_pcg(matrix_pmdcdsc, the_system_level.filter_pres/*, precond*/);

      // enable plotting
      solver->set_plot_mode(Solver::PlotMode::all);

      // set tolerance
      solver->set_tol_rel(1E-5);
      solver->set_max_iter(1000);

      // initialise
      solver->init();

      TimeStamp ts;
      solver->correct(vec_sol_2, vec_rhs);
      time_2 = ts.elapsed_now();

      // release solver
      solver->done();
    }


    comm.print("");
    comm.print("Old   PCG Solve Time: " + stringify_fp_fix(time_1, 3, 9));
    comm.print("New   PCG Solve Time: " + stringify_fp_fix(time_2, 3, 9));
    comm.print("Old vs New Speedup..: " + stringify_fp_fix(time_2 > 1E-12*time_1 ? time_1/time_2 : 0.0, 3, 9));

    // print timings
    comm.print("");
    comm.print(matrix_pmdcdsc.format_timings());

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    //if (args.check("no-err") < 0)
    {
      // Compute and print the H0-/H1-errors
      static constexpr int max_der = SpacePresType::local_degree > 0 ? 1 : 0;
      Assembly::ScalarErrorInfo<DataType> errors_1 = Assembly::ScalarErrorComputer<max_der>::
        compute(vec_sol_1.local(), sol_func, the_domain_level.space_pres, cubature);
      Assembly::ScalarErrorInfo<DataType> errors_2 = Assembly::ScalarErrorComputer<max_der>::
        compute(vec_sol_2.local(), sol_func, the_domain_level.space_pres, cubature);

      // synchronise all local errors
      errors_1.synchronise(comm);
      errors_2.synchronise(comm);

      // print errors
      comm.print("");
      comm.print(errors_1.format_string());
      comm.print(errors_2.format_string());
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    DataType omega = 0.1;

    if(args.check("rich") >= 0)
    {
      comm.print("\nApplying Richardson to SymmetricLumpedSchurMatrix...");
      auto solver = Solver::new_richardson(matrix_gls, the_system_level.filter_pres, omega);

      // enable plotting
      solver->set_plot_mode(Solver::PlotMode::all);

      // set tolerance
      solver->set_tol_rel(1E-5);
      solver->set_max_iter(10);

      // initialise
      solver->init();

      auto vec_sol = vec_rhs.clone();
      vec_sol.format();

      solver->correct(vec_sol, vec_rhs);

      // release solver
      solver->done();
    }

    if(args.check("rich") >= 0)
    {
      comm.print("\nApplying Richardson to PMDCDSCMatrix...");
      auto solver = Solver::new_richardson(matrix_pmdcdsc, the_system_level.filter_pres, omega);

      // enable plotting
      solver->set_plot_mode(Solver::PlotMode::all);

      // set tolerance
      solver->set_tol_rel(1E-5);
      solver->set_max_iter(10);

      // initialise
      solver->init();

      auto vec_sol = vec_rhs.clone();
      vec_sol.format();

      solver->correct(vec_sol, vec_rhs);

      // release solver
      solver->done();
    }

    if(args.check("rich") >= 0)
    {
      comm.print("\nApplying Richardson to PMDCDSCMatrix in ADP Manner...");

      Index glob_dof_offset(0), glob_dof_count(0);

      auto matrix_s = matrix_pmdcdsc.asm_adp_symbolic(glob_dof_offset, glob_dof_count);
      matrix_pmdcdsc.asm_adp_numeric(matrix_s);

      // validate sorting of column indices
      {
        const IndexType* row_ptr_s = matrix_s.row_ptr();
        const IndexType* col_idx_s = matrix_s.col_ind();
        for(Index i = 0; i < matrix_s.rows(); ++i)
        {
          for(IndexType j(row_ptr_s[i]); j + 1 < row_ptr_s[i + 1]; ++j)
          {
            if(col_idx_s[j + 1] <= col_idx_s[j])
            {
              std::cout << "WARNING: ADP Schur matrix row " <<  i << " is unsorted!" << std::endl;
            }
          }
        }
      }

      std::vector<int> all_glob_dof_offset((std::size_t)comm.size());
      std::vector<int> all_glob_dof_counts((std::size_t)comm.size());
      int my_nod = int(matrix_s.rows());
      int my_off = int(glob_dof_offset);

      // allgather offsets + counts
      comm.allgather(&my_off, 1, all_glob_dof_offset.data(), 1);
      comm.allgather(&my_nod, 1, all_glob_dof_counts.data(), 1);

      // allocate a global vector
      auto vec_x = matrix_s.create_vector_r();
      XASSERT(vec_x.size() == glob_dof_count);

      auto vsol_1 = vec_rhs.clone();
      auto vsol_2 = vec_rhs.clone();
      auto vdef_1 = vec_rhs.clone();
      auto vdef_2 = vec_rhs.clone();
      auto verr = vec_rhs.clone();
      vsol_1.format();
      vsol_2.format();

      //typename SystemLevelType::LocalPresVector vtx_sol_1, vtx_sol_2, vtx_def_1, vtx_def_2;
      //typename SystemLevelType::LocalPresVector cell_sol_1, cell_sol_2, cell_def_1, cell_def_2;


      for(int iter(0); iter <= 10; ++iter)
      {
        matrix_pmdcdsc.apply(vdef_1, vsol_1);
        //matrix_pmdcdsc.apply(vdef_1, vsol_1, vec_rhs, -DataType(1));
        //matrix_pmdcdsc.get_local_schur_matrix().apply(vdef_1.local(), vsol_1.local(), vec_rhs.local(), -1.0);

        comm.allgatherv(vsol_2.local().elements(), (std::size_t)my_nod, vec_x.elements(),
          all_glob_dof_counts.data(), all_glob_dof_offset.data());
        //matrix_s.apply(vdef_2.local(), vec_x, vec_rhs.local(), -DataType(1));
        matrix_s.apply(vdef_2.local(), vec_x);

        //Assembly::DiscreteVertexProjector::project(vtx_def_1, vdef_1.local(), the_domain_level.space_pres);
        //Assembly::DiscreteVertexProjector::project(vtx_def_2, vdef_2.local(), the_domain_level.space_pres);
        //Assembly::DiscreteCellProjector::project(cell_def_1, vdef_1.local(), the_domain_level.space_pres);
        //Assembly::DiscreteCellProjector::project(cell_def_2, vdef_2.local(), the_domain_level.space_pres);

        vdef_1.axpy(vdef_1, vec_rhs, -1.0);
        vdef_2.axpy(vdef_2, vec_rhs, -1.0);

        // compute defect norm
        DataType norm_1 = vdef_1.norm2();
        DataType norm_2 = vdef_2.norm2();

        vsol_1.axpy(vdef_1, vsol_1, omega);
        vsol_2.axpy(vdef_2, vsol_2, omega);

        // compute difference of solutions
        verr.axpy(vsol_1, vsol_2, -1.0);
        DataType norm_e = verr.norm2();

        comm.print("ADP-Richardson: " + stringify(iter).pad_front(3) + ": " + stringify_fp_sci(norm_1)
          + " | " + stringify_fp_sci(norm_2)+ " : " + stringify_fp_sci(norm_e));

        /**f(true)
        {
          // build VTK name
          String vtk_name = String("./pmdcdsc");
          //vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
          vtk_name += "-n" + stringify(comm.size());
          vtk_name += "." + stringify(iter).pad_front(2, '0');

          // Create a VTK exporter for our mesh
          Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

          // project velocity and pressure
          Assembly::DiscreteVertexProjector::project(vtx_sol_1, vsol_1.local(), the_domain_level.space_pres);
          Assembly::DiscreteVertexProjector::project(vtx_sol_2, vsol_2.local(), the_domain_level.space_pres);
          //Assembly::DiscreteVertexProjector::project(vtx_sol_3, vsol_3.local(), the_domain_level.space_pres);
          //Assembly::DiscreteVertexProjector::project(vtx_rhs, vec_rhs.local(), the_domain_level.space_pres);
          Assembly::DiscreteCellProjector::project(cell_sol_1, vsol_1.local(), the_domain_level.space_pres);
          Assembly::DiscreteCellProjector::project(cell_sol_2, vsol_2.local(), the_domain_level.space_pres);

          // write velocity
          exporter.add_vertex_scalar("vsol_1", vtx_sol_1.elements());
          exporter.add_vertex_scalar("vsol_2", vtx_sol_2.elements());
          exporter.add_vertex_scalar("vdef_1", vtx_def_1.elements());
          exporter.add_vertex_scalar("vdef_2", vtx_def_2.elements());
          //exporter.add_vertex_scalar("sol_3", vtx_sol_3.elements());
          //exporter.add_vertex_scalar("rhs", vtx_rhs.elements());
          exporter.add_cell_scalar("csol_1", cell_sol_1.elements());
          exporter.add_cell_scalar("csol_2", cell_sol_2.elements());
          exporter.add_cell_scalar("cdef_1", cell_def_1.elements());
          exporter.add_cell_scalar("cdef_2", cell_def_2.elements());

          // finally, write the VTK file
          //comm.print("Writing "+vtk_name);
          exporter.write(vtk_name, comm);
        }*/
      }
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    if (args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name = String("./poisson-mixed");
      vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
      vtk_name += "-n" + stringify(comm.size());

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

      // project velocity and pressure
      typename SystemLevelType::LocalPresVector vtx_sol_1, vtx_sol_2, vtx_rhs;
      Assembly::DiscreteVertexProjector::project(vtx_sol_1, vec_sol_1.local(), the_domain_level.space_pres);
      Assembly::DiscreteVertexProjector::project(vtx_sol_2, vec_sol_2.local(), the_domain_level.space_pres);
      Assembly::DiscreteVertexProjector::project(vtx_rhs, vec_rhs.local(), the_domain_level.space_pres);

      // write velocity
      exporter.add_vertex_scalar("sol_1", vtx_sol_1.elements());
      exporter.add_vertex_scalar("sol_2", vtx_sol_2.elements());
      exporter.add_vertex_scalar("rhs", vtx_rhs.elements());

      // finally, write the VTK file
      comm.print("Writing "+vtk_name);
      exporter.write(vtk_name, comm);
    }
  }

  template<int dim>
  void run_dim(const Dist::Comm& comm, SimpleArgParser& args)
  {
    typedef Shape::Hypercube<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    typedef Space::Lagrange2::Element<TrafoType> AuxSpace;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1> > SpaceType;

    //typedef Space::CroRavRanTur::Element<TrafoType> AuxSpace;
    //typedef Space::Discontinuous::ElementP0<TrafoType> SpaceType;

    // let's create our domain
    comm.print("Dimension: " + stringify(dim));

    // parse levels
    int level = 3;
    args.parse("level", level);

    // let's create our domain
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, AuxSpace, SpaceType> DomainLevelType;
    Control::Domain::UnitCubeDomainControl<DomainLevelType> domain(comm, level, level);

    // plot our levels
    comm.print("Level: " + stringify(domain.max_level_index()) + " [" + stringify(level) + "]");

    // run our application
    run(args, domain);
  }

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // Add the DomainControl's supported arguments
    // check command line arguments
    args.support("dim");
    args.support("debug");
    args.support("bench");
    args.support("pcg");
    args.support("level");
    args.support("vtk");
    args.support("amg");
    args.support("rich");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
      {
        comm.print(std::cerr, "\nERROR: unknown option '--" + (*it).second + "'");
      }

      // abort
      FEAT::Runtime::abort();
    }

    // create a time-stamp
    TimeStamp time_stamp;

    // define our mesh type
    int dim = 2;
    args.parse("dim", dim);
    switch(dim)
    {
    // case 1: run_dim<1>(comm, args); break;
    case 2: run_dim<2>(comm, args); break;
    case 3: run_dim<3>(comm, args); break;
    default:
      comm.print(std::cerr, "\nERROR: invalid dimension");
      Runtime::abort();
    }

    // print elapsed runtime
    comm.print("\nTotal Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }
} // namespace PoissonMixed

int main(int argc, char* argv [])
{
  // initialise
  FEAT::Runtime::initialise(argc, argv);

  PoissonMixed::main(argc, argv);
  // okay
  return FEAT::Runtime::finalise();
}
