// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once
#ifndef APPLICATIONS_CCND_COMMON_HPP
#define APPLICATIONS_CCND_COMMON_HPP 1

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/assembly/burgers_assembly_job.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/assembly/stokes_fbm_assembler.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/direct_stokes_solver.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/checkpoint_control.hpp>
#include <control/statistics.hpp>

#include <vector>

#ifndef FEAT_CCND_APP_DIM
#error FEAT_CCND_APP_DIM must be defined to either 2 or 3
#endif

namespace CCND
{
  using namespace FEAT;

  // define application dimension
  static constexpr int dim = FEAT_CCND_APP_DIM;
  static_assert((dim == 2) || (dim == 3), "invalid dimension");

  // output padding length
  static constexpr std::size_t pad_len = 30u;

  // output padding character
  static constexpr char pad_char = '.';

  // helper function to parse arguments
  template<typename T_>
  T_ parse(SimpleArgParser& args, const String& name, T_ default_value)
  {
    args.parse(name, default_value);
    return default_value;
  }

  // our one and only index type
  typedef Index IndexType;

  // depending on whether FEAT_CCND_APP_QUADMATH is defined,
  // we use quad precision or double precision
#ifdef FEAT_CCND_APP_QUADMATH
#  define Q_(x) (x##Q)
  typedef __float128 DataType;
  static constexpr int fp_num_digs = 35;
  static const char* fp_typename = "quadruple";
#else
#  define Q_(x) x
  typedef double DataType;
  static constexpr int fp_num_digs = 15;
  static const char* fp_typename = "double";
#endif

  // our matrix types
  typedef LAFEM::SparseMatrixBCSR<CCND::DataType, CCND::IndexType, dim, dim> LocalMatrixBlockA;
  typedef LAFEM::SparseMatrixBCSR<CCND::DataType, CCND::IndexType, dim, 1> LocalMatrixBlockB;
  typedef LAFEM::SparseMatrixBCSR<CCND::DataType, CCND::IndexType, 1, dim> LocalMatrixBlockD;
  typedef LAFEM::SparseMatrixCSR<CCND::DataType, CCND::IndexType> LocalScalarMatrix;

  // our vector types
  typedef LAFEM::DenseVectorBlocked<CCND::DataType, CCND::IndexType, dim> LocalVeloVector;
  typedef LAFEM::DenseVector<CCND::DataType, CCND::IndexType> LocalPresVector;

  // define our mesh type and other geometry related classes
  typedef Shape::Hypercube<dim> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType, dim, DataType> MeshType;
  typedef Geometry::MeshPart<MeshType> MeshPartType;
  typedef Geometry::RootMeshNode<MeshType> MeshNodeType;
  typedef Geometry::Atlas::ChartBase<MeshType> ChartBaseType;
  typedef Geometry::MeshAtlas<MeshType> MeshAtlasType;

  // define our trafo type: standard or isoparametric
#ifdef FEAT_CCND_APP_ISOPARAM
  static constexpr bool isoparam = true;
  typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
#else
  static constexpr bool isoparam = false;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
#endif

  // define FE space types
  typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
  typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class DomainLevelBase :
    public Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType>
  {
  public:
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> BaseClass;

    /// the FBM assembler for this level
    std::unique_ptr<Assembly::StokesFBMAssembler<MeshType>> fbm_assembler;

    // inherit constructor
    using BaseClass::BaseClass;

    void create_fbm_assembler(const Dist::Comm& comm, const String& fbm_meshpart_name)
    {
      // get out mesh node
      const MeshNodeType& mesh_node = *this->get_mesh_node();

      // allocate assembler
      fbm_assembler.reset(new Assembly::StokesFBMAssembler<MeshType>(*mesh_node.get_mesh()));

      // find meshpart node
      const auto* part_node = mesh_node.find_mesh_part_node(fbm_meshpart_name);
      XASSERTM(part_node != nullptr, "FBM meshpart node not found");

      // get the meshpart if it exists
      const MeshPartType* fbm_meshpart = part_node->get_mesh();
      if(fbm_meshpart)
        fbm_assembler->add_fbm_meshpart(*fbm_meshpart);

      // synchronize over comm
      fbm_assembler->sync(mesh_node, comm);

      // compile the assembler
      fbm_assembler->compile();
    }
  }; // class DomainLevel

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * \brief Navier-Stokes System Level base class
   */
  class SystemLevelBase :
    public Control::StokesBlockedCombinedSystemLevel<dim, DataType, IndexType, LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockD, LocalScalarMatrix>
  {
  public:
    /// out base class
    typedef Control::StokesBlockedCombinedSystemLevel<dim, DataType, IndexType, LocalMatrixBlockA, LocalMatrixBlockB, LocalMatrixBlockD, LocalScalarMatrix> BaseClass;

    /// the filtered local system matrix for Vanka
    LocalSystemMatrix local_matrix_sys;

    /// the velocity mass matrix
    GlobalMatrixBlockA velo_mass_matrix;

    /// the local interface filter
    LocalVeloUnitFilter filter_interface_fbm;

    /// FBM mask vectors of velocity and pressure
    std::vector<int> fbm_mask_velo, fbm_mask_pres;

    void assemble_velocity_laplace_matrix(Assembly::DomainAssembler<TrafoType>& dom_asm,
      const SpaceVeloType& space_velo, const DataType nu, bool defo, String cubature = "")
    {
      auto& loc_a = this->matrix_a.local();
      loc_a.format();
      Assembly::Common::LaplaceOperatorBlocked<dim> lapl_op;
      Assembly::Common::DuDvOperatorBlocked<dim> dudv_op;
      if(cubature.empty())
        cubature = "auto-degree:" + stringify(2*SpaceVeloType::local_degree+2);
      if(defo)
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_a, dudv_op, space_velo, cubature, nu);
      else
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_a, lapl_op, space_velo, cubature, nu);
    }

    void assemble_velocity_mass_matrix(Assembly::DomainAssembler<TrafoType>& dom_asm, const SpaceVeloType& space_velo, String cubature = "")
    {
      this->velo_mass_matrix = this->matrix_a.clone(LAFEM::CloneMode::Weak);
      auto& loc_m = this->velo_mass_matrix.local();
      loc_m.format();

      if(cubature.empty())
        cubature = "auto-degree:" + stringify(2*SpaceVeloType::local_degree+2);
      Assembly::Common::IdentityOperatorBlocked<dim> id_op;
      Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_m, id_op, space_velo, cubature);
    }

    void compile_local_matrix()
    {
      // convert local matrices
      this->local_matrix_sys.block_a() = this->matrix_a.convert_to_1();
      this->local_matrix_sys.block_b() = this->matrix_b.local().clone(LAFEM::CloneMode::Weak);
      this->local_matrix_sys.block_d() = this->matrix_d.local().clone(LAFEM::CloneMode::Weak);

      // apply velocity unit filters to A and B
      for(const auto& filter : this->get_local_velo_unit_filter_seq())
      {
        filter.second.filter_mat(this->local_matrix_sys.block_a());
        filter.second.filter_offdiag_row_mat(this->local_matrix_sys.block_b());
      }
    }
  }; // class SystemLevelBase

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  struct Counts
  {
    static constexpr std::size_t velo_dofs = 0u;        // number of velocity dofs nodes
    static constexpr std::size_t pres_dofs = 1u;        // number of pressure dofs
    static constexpr std::size_t total_dofs = 2u;       // total dofs = dim*velo_dofs + pres_dofs
    static constexpr std::size_t nonlin_iter = 3u;      // number of nonlinear iterations (newton/picard)
    static constexpr std::size_t linsol_iter = 4u;      // number of linear solver iterations
    static constexpr std::size_t nnze_a = 5u;           // number of non-zero blocks in matrix-block A
    static constexpr std::size_t nnze_b = 6u;           // number of non-zero blocks in matrix-block B
    static constexpr std::size_t nnze_total = 7u;       // total number of non-zero entries in matrix
    static constexpr std::size_t vanka_data = 8u;       // total number of data entries for Vanka
    static constexpr std::size_t elements = 9u;         // total number of mesh elements
    static constexpr std::size_t count = 10u;
  };

  struct Times
  {
    static constexpr std::size_t total_run = 0u;        // total runtime
    static constexpr std::size_t nonlin_total = 1u;     // nonlinear solver time
    static constexpr std::size_t nonlin_asm_def = 2u;   // nonlinear defect assembly time
    static constexpr std::size_t nonlin_asm_mat = 3u;   // nonlinear matrix assembly time
    static constexpr std::size_t linsol_init = 4u;      // linear solver init time
    static constexpr std::size_t linsol_apply = 5u;     // linear solver apply time
    static constexpr std::size_t mg_defect = 6u;        // multigrid defect compute time
    static constexpr std::size_t mg_smooth = 7u;        // multigrid smooth apply time
    static constexpr std::size_t mg_coarse = 8u;        // multigrid coarse grid solve time
    static constexpr std::size_t mg_transfer = 9u;      // multigrid grid transfer time
    static constexpr std::size_t stokes_solve = 10u;    // stokes solve time
    static constexpr std::size_t vanka_init_sym = 11u;  // vanka symbolic factorization time
    static constexpr std::size_t vanka_init_num = 12u;  // vanka numeric factorization time
    static constexpr std::size_t vanka_apply = 13u;     // vanka apply time
    static constexpr std::size_t vtk_write = 14u;       // VTK output time
    static constexpr std::size_t sol_analysis = 15u;    // solution analysis time
    static constexpr std::size_t sol_init_read = 16u;   // initial solution read time
    static constexpr std::size_t sol_final_write = 17u; // final solution write time
    static constexpr std::size_t checkpoint_save = 18u; // checkpoint save time
    static constexpr std::size_t checkpoint_load = 19u; // checkpoint load time
    static constexpr std::size_t create_domain = 20u;   // create domain time
    static constexpr std::size_t create_system = 21u;   // create system time
    static constexpr std::size_t create_solver = 22u;   // create solver time
    static constexpr std::size_t count = 23u;
  };

  struct Bytes
  {
    static constexpr std::size_t peak_p = 0u;           // peak physical memory usage
    static constexpr std::size_t peak_v = 1u;           // peak virtual memory usage
    static constexpr std::size_t mesh = 2u;             // mesh node memory usage
    static constexpr std::size_t gate = 3u;             // gate memory usage
    static constexpr std::size_t muxer = 4u;            // muxer memory usage
    static constexpr std::size_t matrix = 5u;           // system matrices memory usage
    static constexpr std::size_t matrix_struct = 6u;    // matrix structures memory usage
    static constexpr std::size_t matrix_values = 7u;    // matrix values memory usage
    static constexpr std::size_t transfer = 8u;         // transfers operators memory usage
    static constexpr std::size_t vanka = 9u;            // vanka memory usage
    static constexpr std::size_t count = 10u;
  };

  // class for collection benchmark statistics
  class BenchmarkStats
  {
  public:
    // timings (in seconds)
    double times[Times::count];     // of this process
    double times_min[Times::count]; // minimum over all processes
    double times_max[Times::count]; // maximum over all processes

    // counts
    unsigned long long counts[Counts::count];     // of this process
    unsigned long long counts_sum[Counts::count]; // summed up over all processes
    //unsigned long long counts_min[Counts::count]; // minimum over all processes
    //unsigned long long counts_max[Counts::count]; // maximum over all processes

    // sizes (in bytes)
    unsigned long long bytes[Bytes::count];     // of this process
    unsigned long long bytes_sum[Bytes::count]; // summed up over all processes
    unsigned long long bytes_min[Bytes::count]; // minimum over all processes
    unsigned long long bytes_max[Bytes::count]; // maximum over all processes

    // --------------------------------------------------------------

    // formatted output padding length
    static constexpr std::size_t padlen = std::size_t(30);

    BenchmarkStats()
    {
      for(std::size_t i(0); i < Counts::count; ++i)
        counts[i] = std::size_t(0);
      for(std::size_t i(0); i < Times::count; ++i)
        times[i] = std::size_t(0);
      for(std::size_t i(0); i < Bytes::count; ++i)
        bytes[i] = std::size_t(0);
    }

    String format() const
    {
      String s;
      const char pc = '.';

      s = "Basic Statistics:\n";
      s += String("Mesh Elements").pad_back(padlen, pc) + ": " + stringify(counts_sum[Counts::elements]) + "\n";
      s += String("Velocity Dof Nodes").pad_back(padlen, pc) + ": " + stringify(counts[Counts::velo_dofs]) + "\n";
      s += String("Pressure Dofs").pad_back(padlen, pc) + ": " + stringify(counts[Counts::pres_dofs]) + "\n";
      s += String("Total Dofs").pad_back(padlen, pc) + ": " + stringify(counts[Counts::total_dofs]) + "\n";
      s += String("Nonzero Blocks A").pad_back(padlen, pc) + ": " + stringify(counts_sum[Counts::nnze_a]) + "\n";
      s += String("Nonzero Blocks B").pad_back(padlen, pc) + ": " + stringify(counts_sum[Counts::nnze_b]) + "\n";
      s += String("Total Nonzero Entries").pad_back(padlen, pc) + ": " + stringify(counts_sum[Counts::nnze_total]) + "\n";
      s += String("Vanka Nonzero Entries").pad_back(padlen, pc) + ": " + stringify(counts_sum[Counts::vanka_data]) + "\n";

      // solver statistics
      s += "\nSolver Statistics:\n";
      s += String("Multigrid Iterations").pad_back(padlen, pc) + ": " + stringify(counts[Counts::linsol_iter]) + "\n";
      s += String("Nonlinear Iterations").pad_back(padlen, pc) + ": " + stringify(counts[Counts::nonlin_iter]) + "\n";

      // append timing info
      s += String("\nTiming Statistics:\n");
      s += String("Total Runtime").pad_back(padlen, pc) + ": " + stringify_fp_fix(times[Times::total_run], 3, 10) + " sec\n";
      s += format_subtime("Domain Create", times[Times::create_domain]);
      s += format_subtime("System Create", times[Times::create_system]);
      s += format_subtime("Solver Create", times[Times::create_solver]);
      s += format_subtime("Initial Stokes Solver", times[Times::stokes_solve]);
      s += format_subtime("Navier-Stokes Solver Total", times[Times::nonlin_total]);
      s += format_subtime("Nonlinear Defect Assembly", times[Times::nonlin_asm_def]);
      s += format_subtime("Nonlinear Matrix Assembly", times[Times::nonlin_asm_mat]);
      s += format_subtime("Nonlinear Precond Initialize", times[Times::linsol_init]);
      s += format_subtime("Nonlinear Precond Apply", times[Times::linsol_apply]);
      s += format_subtime("Multigrid Defect Compute", times[Times::mg_defect]);
      s += format_subtime("Multigrid Smooth Apply", times[Times::mg_smooth]);
      s += format_subtime("Multigrid Coarse Solve", times[Times::mg_coarse]);
      s += format_subtime_mm("Multigrid Transfer Apply", times[Times::mg_transfer], times_min[Times::mg_transfer], times_max[Times::mg_transfer]);
      s += format_subtime_mm("Vanka Symbolic Initialize", times[Times::vanka_init_sym], times_min[Times::vanka_init_sym], times_max[Times::vanka_init_sym]);
      s += format_subtime_mm("Vanka Numeric Initialize", times[Times::vanka_init_num], times_min[Times::vanka_init_num], times_max[Times::vanka_init_num]);
      s += format_subtime_mm("Vanka Local Apply", times[Times::vanka_apply], times_min[Times::vanka_apply], times_max[Times::vanka_apply]);
      s += format_subtime("Solution Analysis", times[Times::sol_analysis]);
      s += format_subtime("VTK Write", times[Times::vtk_write]);
      s += format_subtime("Initial Solution Read", times[Times::sol_init_read]);
      s += format_subtime("Final Solution Write", times[Times::sol_final_write]);
      s += format_subtime("Checkpoint Save", times[Times::checkpoint_save]);
      s += format_subtime("Checkpoint Load", times[Times::checkpoint_load]);

      // append memory info
      s += String("\nMemory Statistics:\n");
      s += String("Peak Physical Memory").pad_back(padlen,pc) + ": "
        + stringify_fp_fix(double(bytes_sum[Bytes::peak_p]) / (1024.0*1024.0), 3, 12) + " MB             { "
        + stringify_fp_fix(double(bytes_min[Bytes::peak_p]) / (1024.0*1024.0), 3, 12) + " MB / "
        + stringify_fp_fix(double(bytes_max[Bytes::peak_p]) / (1024.0*1024.0), 3, 12) + " MB }\n";
      s += String("Peak Virtual Memory").pad_back(padlen,pc) + ": "
        + stringify_fp_fix(double(bytes_sum[Bytes::peak_v]) / (1024.0*1024.0), 3, 12) + " MB             { "
        + stringify_fp_fix(double(bytes_min[Bytes::peak_v]) / (1024.0*1024.0), 3, 12) + " MB / "
        + stringify_fp_fix(double(bytes_max[Bytes::peak_v]) / (1024.0*1024.0), 3, 12) + " MB }\n";
      s += format_submemory_mm("Mesh-Node Size", bytes_sum[Bytes::mesh], bytes_min[Bytes::mesh], bytes_max[Bytes::mesh]);
      s += format_submemory_mm("Gate Size", bytes_sum[Bytes::gate], bytes_min[Bytes::gate], bytes_max[Bytes::gate]);
      s += format_submemory_mm("Muxer Size", bytes_sum[Bytes::muxer], bytes_min[Bytes::muxer], bytes_max[Bytes::muxer]);
      s += format_submemory_mm("Matrix Total Size", bytes_sum[Bytes::matrix], bytes_min[Bytes::matrix], bytes_max[Bytes::matrix]);
      s += format_submemory_mm("Matrix Struct Size", bytes_sum[Bytes::matrix_struct], bytes_min[Bytes::matrix_struct], bytes_max[Bytes::matrix_struct]);
      s += format_submemory_mm("Matrix Values Size", bytes_sum[Bytes::matrix_values], bytes_min[Bytes::matrix_values], bytes_max[Bytes::matrix_values]);
      s += format_submemory_mm("Transfer Size", bytes_sum[Bytes::transfer], bytes_min[Bytes::transfer], bytes_max[Bytes::transfer]);
      s += format_submemory_mm("Vanka Size", bytes_sum[Bytes::vanka], bytes_min[Bytes::vanka], bytes_max[Bytes::vanka]);

      return s;
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(times, times_min, Times::count, Dist::op_min);
      comm.allreduce(times, times_max, Times::count, Dist::op_max);
      comm.allreduce(counts, counts_sum, Counts::count, Dist::op_sum);
      //comm.allreduce(counts, counts_min, Counts::count, Dist::op_min);
      //comm.allreduce(counts, counts_max, Counts::count, Dist::op_max);
      comm.allreduce(bytes, bytes_sum, Bytes::count, Dist::op_sum);
      comm.allreduce(bytes, bytes_min, Bytes::count, Dist::op_min);
      comm.allreduce(bytes, bytes_max, Bytes::count, Dist::op_max);
    }

    String format_submemory(String s, unsigned long long m) const
    {
      return format_submemory(s, m, bytes_sum[Bytes::peak_p]);
    }

    String format_submemory(String s, unsigned long long m, unsigned long long total) const
    {
      double mb = double(m) / (1024.0*1024.0);
      return s.pad_back(padlen, '.') + ": " + stringify_fp_fix(mb, 3, 12) + " MB [" + stringify_fp_fix(100.0*double(m)/double(total), 3, 7) + "% ]\n";
    }

    String format_submemory_mm(String s, unsigned long long m, unsigned long long mmin, unsigned long long mmax) const
    {
      return format_submemory_mm(s, m, bytes_sum[Bytes::peak_p], mmin, mmax);
    }

    String format_submemory_mm(String s, unsigned long long m, unsigned long long total, unsigned long long mmin, unsigned long long mmax) const
    {
      double mb = double(m) / (1024.0*1024.0);
      double mmi = double(mmin) / (1024.0*1024.0);
      double mma = double(mmax) / (1024.0*1024.0);
      return s.pad_back(padlen, '.') + ": " + stringify_fp_fix(mb, 3, 12) + " MB [" + stringify_fp_fix(100.0*double(m)/double(total), 3, 7)
        + "% ] { " + stringify_fp_fix(mmi, 3, 12) + " MB / " + stringify_fp_fix(mma, 3, 12) + " MB }\n";
    }

    String format_subtime(String s, double t) const
    {
      return format_subtime(s, t, times[Times::total_run]);
    }

    String format_subtime(String s, double t, double total) const
    {
      return s.pad_back(padlen, '.') + ": " + stringify_fp_fix(t, 3, 10) + " sec [" + stringify_fp_fix(100.0*t/total, 3, 7) + "% ]\n";
    }

    String format_subtime_mm(String s, double t, double tmin, double tmax) const
    {
      return format_subtime_mm(s, t, times[Times::total_run], tmin, tmax);
    }

    String format_subtime_mm(String s, double t, double total, double tmin, double tmax) const
    {
      return s.pad_back(padlen, '.') + ": " + stringify_fp_fix(t, 3, 10) + " sec [" + stringify_fp_fix(100.0*t/total, 3, 7)
        + "% ] { " + stringify_fp_fix(tmin, 3, 10) + " / " + stringify_fp_fix(tmax, 3, 10) + " }\n";
    }
  }; // struct BenchmarkStats<...>
} // namespace CCND

#endif // APPLICATIONS_CCND_COMMON_HPP
