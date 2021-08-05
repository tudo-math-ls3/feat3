#pragma once
#ifndef AREA51_CCND_FIBER_COMMON_HPP
#define AREA51_CCND_FIBER_COMMON_HPP 1

#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <area51/ccnd_fiber/orientation_moment_burgers_assembler.hpp>
#include <kernel/assembly/burgers_assembler.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/vanka.hpp>
#include <kernel/solver/bicgstab.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/umfpack.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/analytic/common.hpp>
// FEAT-Cubature includes
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory

#include<kernel/assembly/fe_interpolator.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/checkpoint_control.hpp>
#include <control/statistics.hpp>

#include <vector>
#include <area51/ccnd_fiber/ccnd_steady_function.hpp>


namespace CCND_FIBER
{
  using namespace FEAT;

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
    static constexpr std::size_t checkpoint = 15u;      // checkpoint write time
    static constexpr std::size_t sol_analysis = 16u;    // solution analysis time
    static constexpr std::size_t tensor_load = 17u;     //time to load tensors
    static constexpr std::size_t basic_matrix_asm = 18u;//time to assemble the basic matrices for stokes
    static constexpr std::size_t count = 19u;
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
    static constexpr std::size_t tensor = 10u;          // tensor memory usage
    static constexpr std::size_t count = 11u;
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
      s += format_subtime("Stokes Solver", times[Times::stokes_solve]);
      s += format_subtime("Nonlinear Solver Total", times[Times::nonlin_total]);
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
      s += format_subtime_mm("Tensor loading", times[Times::tensor_load], times_min[Times::tensor_load], times_max[Times::tensor_load]);
      s += format_subtime_mm("Basic Stokes Matrix Assembly", times[Times::basic_matrix_asm], times_min[Times::basic_matrix_asm], times_max[Times::basic_matrix_asm]);
      s += format_subtime("VTK Write", times[Times::vtk_write]);
      s += format_subtime("Checkpointing", times[Times::checkpoint]);
      s += format_subtime("Solution Analysis", times[Times::sol_analysis]);

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
      s += format_submemory_mm("Tensor Size", bytes_sum[Bytes::tensor], bytes_min[Bytes::tensor], bytes_max[Bytes::tensor]);

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
  }; // struct BenchmarkSummary<...>


  // helper function to parse arguments
  template<typename T_>
  T_ parse(SimpleArgParser& args, const String& name, T_ default_value)
  {
    args.parse(name, default_value);
    return default_value;
  }


  typedef double DataType;
  // our one and only memory type
  typedef Mem::Main MemType;
  // our one and only index type
  typedef Index IndexType;



/**
   * \brief modified Navier-Stokes System Level class
   */
  template<
    int dim_,
    typename MemType_ = Mem::Main,
    typename DataType_ = DataType,
    typename IndexType_ = Index,
    typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
    typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
    typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
    typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>>
  class ModNavierStokesBlockedSystemLevel :
    public Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
  {
  public:
    typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;

    // the noflow velocity filter, used in bench3 simulation
    typename BaseClass::LocalVeloFilter local_velo_filter_noflow;

    // the filtered local system matrix for Vanka
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    // the local velocity mass matrix
    MatrixBlockA_ local_velo_mass_matrix;

    //the standard BurgersAssembler
    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType alpha = DataType(0))
    {
      Assembly::BurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = DataType(0);
      burgers_mat.theta = alpha;
      burgers_mat.deformation = true;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, space_velo, cubature);
    }

    // ????!
//     //the standard BurgersAssembler
//     template<typename SpaceV_, typename Cubature_>
//     void assemble_velocity_chaning_defo_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ defo_threshold, const DataType alpha = DataType(0))
//     {
//       Assembly::ChangeGradToDefoBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
//       burgers_mat.nu = nu;
//       burgers_mat.beta = DataType(0);
//       burgers_mat.theta = alpha;
//
//       auto& loc_a = this->matrix_a.local();
//       loc_a.format();
//
//       // create dummy convection vector
//       auto vec_c = loc_a.create_vector_l();
//       vec_c.format();
//
//       burgers_mat.assemble_matrix(loc_a, vec_c, space_velo, cubature, defo_threshold);
//     }



    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const Tiny::Matrix<DataType_, dim_, dim_>& orientation, const DataType alpha = DataType(0.))
    {
      Assembly::ModBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, orientation, space_velo, cubature);
    }

//     template<typename SpaceV_, typename Cubature_>
//     void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType alpha = DataType(0.))
//     {
//       Assembly::ModBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
//       burgers_mat.nu = nu;
//       burgers_mat.beta = 0.0;
//       burgers_mat.theta = alpha;
//       burgers_mat.N_s = N_s;
//
//       //create placeholder orientation matrix:
//       Tiny::Matrix<DataType_, dim_, dim_> orientation;
//       orientation.format();
//       for(int i = 0; i < dim_; ++i)
//       {
//         orientation(i, i) = DataType_(1);
//       }
//
//       auto& loc_a = this->matrix_a.local();
//       loc_a.format();
//
//       // create dummy convection vector
//       auto vec_c = loc_a.create_vector_l();
//       vec_c.format();
//
//       burgers_mat.assemble_matrix(loc_a, vec_c, orientation, space_velo, cubature);
//     }

    //for Orientation matrix
    template<typename SpaceV_, typename Cubature_, typename Function_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const Function_& func, const DataType alpha = DataType(0.))
    {
      Assembly::ModBurgersSDAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, func, space_velo, cubature);
    }

    //for Full Tensor formulation matrix
    template<typename SpaceV_, typename Cubature_, typename Function2_, typename Function4_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType_ N_p, const Function2_& func2, const Function4_& func4, const DataType alpha = DataType(0.))
    {
      Assembly::FullTensorBurgersSDAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;
      burgers_mat.N_p = N_p;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, func2, func4, space_velo, cubature);
    }

    //for discrete Velo FE Full Tensor formulation matrix
    template<typename SpaceV_, typename Cubature_, typename VeloFunction_, typename Function2_, typename Function4_>
    void assemble_discrete_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType_ N_p, const VeloFunction_& velo_func, const IndexType_ level_diff, const Function2_& func2, const Function4_& func4, const DataType alpha = DataType(0.))
    {
      Assembly::FullDiscreteTensorBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;
      burgers_mat.N_p = N_p;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, velo_func, level_diff, func2, func4, space_velo, cubature);
    }

    //for discrete Velo FE Full Tensor formulation matrix
    template<typename SpaceV_, typename Cubature_, typename Function2_, typename Function4_>
    void assemble_truncated_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType_ N_p, const LAFEM::DenseVectorBlocked<MemType_, DataType_, IndexType_, dim_>& velo_vec, const Function2_& func2, const Function4_& func4, const DataType_ alpha = DataType(0.))
    {
      Assembly::FullFETensorBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;
      burgers_mat.N_p = N_p;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, velo_vec, func2, func4, space_velo, cubature);
    }

    //for fiber orientation formulation matrix
    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_laplace_matrix_fiber(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType_ N_p, const LAFEM::DenseVectorBlocked<MemType_, DataType_, IndexType_, dim_*(dim_+1)/2>& tensor2, const LAFEM::DenseVectorBlocked<MemType_, DataType_, IndexType_, dim_*(dim_+1)*(dim_+2)*(dim_+3)/24>& tensor4, const DataType alpha = DataType(0.))
    {
      Assembly::FullFiberOrientationTensorBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;
      burgers_mat.N_p = N_p;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, tensor2, tensor4, space_velo, cubature);
    }

    //for Viscosity
    template<typename SpaceV_, typename Cubature_, typename Function_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const Function_& func, const DataType alpha = DataType(0.))
    {
      Assembly::ModBurgersSDViscosityAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, func, space_velo, cubature);
    }

    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_mass_matrix(const SpaceV_& space_velo, const Cubature_& cubature)
    {
      local_velo_mass_matrix = this->matrix_a.local().clone(LAFEM::CloneMode::Weak);
      local_velo_mass_matrix.format();

      Assembly::Common::IdentityOperatorBlocked<dim_> id_op;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(local_velo_mass_matrix, id_op, space_velo, cubature);
    }

    void compile_local_matrix()
    {
      // convert local matrices
      this->local_matrix_sys.block_a() = this->matrix_a.convert_to_1();
      this->local_matrix_sys.block_b() = this->matrix_b.local().clone(LAFEM::CloneMode::Weak);
      this->local_matrix_sys.block_d() = this->matrix_d.local().clone(LAFEM::CloneMode::Weak);

      // apply filter to A and B
      this->filter_velo.local().filter_mat(this->local_matrix_sys.block_a());
      this->filter_velo.local().filter_offdiag_row_mat(this->local_matrix_sys.block_b());
    }
  }; // class ModNavierStokesBlockedSystemLevel


/**
   * \brief modified Navier-Stokes System Level class
   */
  template<
    int dim_,
    typename MemType_ = Mem::Main,
    typename DataType_ = DataType,
    typename IndexType_ = Index,
    typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
    typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
    typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
    typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>>
  class ModNavierStokesBlockedMeanPresSystemLevel :
    public Control::StokesBlockedUnitVeloMeanPresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
  {
  public:
    typedef Control::StokesBlockedUnitVeloMeanPresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;

    // the noflow velocity filter, used in bench3 simulation
    typename BaseClass::LocalVeloFilter local_velo_filter_noflow;
//     // the mean pres filter
//     typename BaseClass::LocalPresFilter local_pres_filter_mean;

    // the filtered local system matrix for Vanka
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    // the local velocity mass matrix
    MatrixBlockA_ local_velo_mass_matrix;

    //the standard BurgersAssembler
    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType alpha = DataType(0))
    {
      Assembly::BurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = DataType(0);
      burgers_mat.theta = alpha;
      burgers_mat.deformation = true;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, space_velo, cubature);
    }

    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const Tiny::Matrix<DataType_, dim_, dim_>& orientation, const DataType alpha = DataType(0.))
    {
      Assembly::ModBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, orientation, space_velo, cubature);
    }

//     template<typename SpaceV_, typename Cubature_>
//     void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType alpha = DataType(0.))
//     {
//       Assembly::ModBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
//       burgers_mat.nu = nu;
//       burgers_mat.beta = 0.0;
//       burgers_mat.theta = alpha;
//       burgers_mat.N_s = N_s;
//
//       //create placeholder orientation matrix:
//       Tiny::Matrix<DataType_, dim_, dim_> orientation;
//       orientation.format();
//       for(int i = 0; i < dim_; ++i)
//       {
//         orientation(i, i) = DataType_(1);
//       }
//
//       auto& loc_a = this->matrix_a.local();
//       loc_a.format();
//
//       // create dummy convection vector
//       auto vec_c = loc_a.create_vector_l();
//       vec_c.format();
//
//       burgers_mat.assemble_matrix(loc_a, vec_c, orientation, space_velo, cubature);
//     }

    //for Orientation matrix
    template<typename SpaceV_, typename Cubature_, typename Function_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const Function_& func, const DataType alpha = DataType(0.))
    {
      Assembly::ModBurgersSDAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, func, space_velo, cubature);
    }

    //for Full Tensor formulation matrix
    template<typename SpaceV_, typename Cubature_, typename Function2_, typename Function4_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType_ N_p, const Function2_& func2, const Function4_& func4, const DataType alpha = DataType(0.))
    {
      Assembly::FullTensorBurgersSDAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;
      burgers_mat.N_p = N_p;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, func2, func4, space_velo, cubature);
    }

    //for discrete Velo FE Full Tensor formulation matrix
    template<typename SpaceV_, typename Cubature_, typename VeloFunction_, typename Function2_, typename Function4_>
    void assemble_discrete_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType_ N_p, const VeloFunction_& velo_func, const IndexType_ level_diff, const Function2_& func2, const Function4_& func4, const DataType alpha = DataType(0.))
    {
      Assembly::FullDiscreteTensorBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;
      burgers_mat.N_p = N_p;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, velo_func, level_diff, func2, func4, space_velo, cubature);
    }

    //for discrete Velo FE Full Tensor formulation matrix
    template<typename SpaceV_, typename Cubature_, typename Function2_, typename Function4_>
    void assemble_truncated_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType_ N_p, const LAFEM::DenseVectorBlocked<MemType_, DataType_, IndexType_, dim_>& velo_vec, const Function2_& func2, const Function4_& func4, const DataType_ alpha = DataType(0.))
    {
      Assembly::FullFETensorBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;
      burgers_mat.N_p = N_p;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, velo_vec, func2, func4, space_velo, cubature);
    }

    //for fiber orientation formulation matrix
    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_laplace_matrix_fiber(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const DataType_ N_s, const DataType_ N_p, const LAFEM::DenseVectorBlocked<MemType_, DataType_, IndexType_, dim_*(dim_+1)/2>& tensor2, const LAFEM::DenseVectorBlocked<MemType_, DataType_, IndexType_, dim_*(dim_+1)*(dim_+2)*(dim_+3)/24>& tensor4, const DataType alpha = DataType(0.))
    {
      Assembly::FullFiberOrientationTensorBurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;
      burgers_mat.N_s = N_s;
      burgers_mat.N_p = N_p;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, tensor2, tensor4, space_velo, cubature);
    }

    //for Viscosity
    template<typename SpaceV_, typename Cubature_, typename Function_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, const Function_& func, const DataType alpha = DataType(0.))
    {
      Assembly::ModBurgersSDViscosityAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = alpha;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, func, space_velo, cubature);
    }

    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_mass_matrix(const SpaceV_& space_velo, const Cubature_& cubature)
    {
      local_velo_mass_matrix = this->matrix_a.local().clone(LAFEM::CloneMode::Weak);
      local_velo_mass_matrix.format();

      Assembly::Common::IdentityOperatorBlocked<dim_> id_op;
      Assembly::BilinearOperatorAssembler::assemble_matrix1(local_velo_mass_matrix, id_op, space_velo, cubature);
    }

    void compile_local_matrix()
    {
      // convert local matrices
      this->local_matrix_sys.block_a() = this->matrix_a.convert_to_1();
      this->local_matrix_sys.block_b() = this->matrix_b.local().clone(LAFEM::CloneMode::Weak);
      this->local_matrix_sys.block_d() = this->matrix_d.local().clone(LAFEM::CloneMode::Weak);

      // apply filter to A and B
      this->filter_velo.local().filter_mat(this->local_matrix_sys.block_a());
      this->filter_velo.local().filter_offdiag_row_mat(this->local_matrix_sys.block_b());
    }
  }; // class ModNavierStokesBlockedMeanPresSystemLevel
}
#endif
