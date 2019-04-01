// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// The nD Steady Navier-Stokes CC-Q2/P1dc Toy-Code Benchmark Solver (TM)
// ---------------------------------------------------------------------
// This application implements a "simple" parallel steady Navier-Stokes solver
// using the monolithic "CC" approach with Q2/P1dc space to solve the infamous
// steady-state "flow around a cylinder" benchmark.
//
// !!!!!!!!!!!!!!!!!!!!!!!
// !!!  W A R N I N G  !!!
// !!!!!!!!!!!!!!!!!!!!!!!
// This application is a "toy code" benchmark solver, i.e. it is meant as a playground
// for the HPC guys to tweak their parallel Oseen solvers for more interesting scenarios
// than Stokes on the unit-square. Do *NOT* use this application as a base for serious
// FEM analysis work !!!
// </Warning>
//
// \author Peter Zajac
//
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
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

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>
#include <control/statistics.hpp>

#include <vector>

namespace NvSCCNDQ2P1dc
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
    static constexpr std::size_t vanka_init_sym = 11u;  // vanka symbolic factorisation time
    static constexpr std::size_t vanka_init_num = 12u;  // vanka numeric factorisation time
    static constexpr std::size_t vanka_apply = 13u;     // vanka apply time
    static constexpr std::size_t count = 14u;
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

  template<typename DT_, int dim_>
  struct BenchmarkSummary
  {
    // post-processing data

    // drag/lift forces
    DT_ drag_coeff, lift_coeff, drag_err, lift_err;

    // pressure values
    DT_ pres_diff, pres_err;

    // flow through upper/lower region
    DT_ flux_upper, flux_lower;

    // velocity field information
    Assembly::VelocityInfo<DT_, dim_> velo_info;

    // --------------------------------------------------------------

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

    BenchmarkSummary() :
      drag_coeff(0.0), lift_coeff(0.0), drag_err(0.0), lift_err(0.0),
      pres_diff(0.0), pres_err(0.0),
      flux_upper(0.0), flux_lower(0.0)
    {
      for(std::size_t i(0); i < Counts::count; ++i)
        counts[i] = std::size_t(0);
      for(std::size_t i(0); i < Times::count; ++i)
        times[i] = std::size_t(0);
      for(std::size_t i(0); i < Bytes::count; ++i)
        bytes[i] = std::size_t(0);
    }

    String format(int prec = 16) const
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

      // append coefficients and velocity info
      s += "\nSolution Analysis:\n";
      s += String("Drag Coefficient").pad_back(padlen, pc) + ": " + stringify_fp_fix(drag_coeff, prec)
        + "   [ Error: " + stringify_fp_sci(drag_err, prec) + " ]\n";
      s += String("Lift Coefficient").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff, prec)
        + "   [ Error: " + stringify_fp_sci(lift_err, prec) + " ]\n";
      s += String("Pressure Difference").pad_back(padlen, pc) + ": " + stringify_fp_fix(pres_diff, prec)
        + "   [ Error: " + stringify_fp_sci(pres_err, prec) + " ]\n";
      s += String("Upper Flux").pad_back(padlen, pc) + ": " + stringify_fp_fix(flux_upper, prec) + "\n";
      s += String("Lower Flux").pad_back(padlen, pc) + ": " + stringify_fp_fix(flux_lower, prec) + "\n";
      s += velo_info.format_string(prec, padlen, pc);

      // append timing info
      s += String("\nTiming Statistics:\n");
      s += String("Total Runtime").pad_back(padlen, pc) + ": " + stringify_fp_fix(times[Times::total_run], 3, 10) + " sec\n";
      s += format_subtime("Stokes Solver", times[Times::stokes_solve]);
      s += format_subtime("Nonlinear Solver Total", times[Times::nonlin_total]);
      s += format_subtime("Nonlinear Defect Assembly", times[Times::nonlin_asm_def]);
      s += format_subtime("Nonlinear Matrix Assembly", times[Times::nonlin_asm_mat]);
      s += format_subtime("Nonlinear Precond Initialise", times[Times::linsol_init]);
      s += format_subtime("Nonlinear Precond Apply", times[Times::linsol_apply]);
      s += format_subtime("Multigrid Defect Compute", times[Times::mg_defect]);
      s += format_subtime("Multigrid Smooth Apply", times[Times::mg_smooth]);
      s += format_subtime("Multigrid Coarse Solve", times[Times::mg_coarse]);
      s += format_subtime_mm("Multigrid Transfer Apply", times[Times::mg_transfer], times_min[Times::mg_transfer], times_max[Times::mg_transfer]);
      s += format_subtime_mm("Vanka Symbolic Initialise", times[Times::vanka_init_sym], times_min[Times::vanka_init_sym], times_max[Times::vanka_init_sym]);
      s += format_subtime_mm("Vanka Numeric Initialise", times[Times::vanka_init_num], times_min[Times::vanka_init_num], times_max[Times::vanka_init_num]);
      s += format_subtime_mm("Vanka Local Apply", times[Times::vanka_apply], times_min[Times::vanka_apply], times_max[Times::vanka_apply]);

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
  }; // struct BenchmarkSummary<...>

  template<int dim_>
  class InflowFunction;

  template<>
  class InflowFunction<2> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 2;
    typedef Analytic::Image::Vector<2> ImageType;
    static constexpr bool can_value = true;

  protected:
    Real _vmax;

  public:
    explicit InflowFunction(Real vmax) :
      _vmax(vmax)
    {
    }

    template<typename Traits_>
    class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
    {
    protected:
      typedef typename Traits_::DataType DataType;
      typedef typename Traits_::PointType PointType;
      typedef typename Traits_::ValueType ValueType;

      const DataType _vmax, _d, _den;

    public:
      explicit Evaluator(const InflowFunction& function) :
        _vmax(function._vmax),
        _d(DataType(0.41)),
        _den(_d*_d)
      {
      }

      ValueType value(const PointType& point)
      {
        const DataType y = point[1];

        ValueType val;
        val[0] = (_vmax * DataType(4) * y * (_d - y)) / _den;
        val[1] = DataType(0);
        return val;
      }
    };
  }; // class InflowFunction<2>

  template<>
  class InflowFunction<3> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 3;
    typedef Analytic::Image::Vector<3> ImageType;
    static constexpr bool can_value = true;

  protected:
    Real _vmax;

  public:
    explicit InflowFunction(Real vmax) :
      _vmax(vmax)
    {
    }

    template<typename Traits_>
    class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
    {
    protected:
      typedef typename Traits_::DataType DataType;
      typedef typename Traits_::PointType PointType;
      typedef typename Traits_::ValueType ValueType;

      const DataType _vmax, _d, _den;

    public:
      explicit Evaluator(const InflowFunction& function) :
        _vmax(function._vmax),
        _d(DataType(0.41)),
        _den(_d*_d*_d*_d)
      {
      }

      ValueType value(const PointType& point)
      {
        const DataType y = point[1];
        const DataType z = point[2];

        ValueType val;
        val[0] = (_vmax * DataType(16) * y * (_d - y) * z * (_d - z)) / _den;
        val[1] = val[2] = DataType(0);
        return val;
      }
    };
  }; // class InflowFunction

  template<typename DataType_>
  class BenchBodyForceAccumulator
  {
  private:
    const bool _defo;
    const DataType_ _nu;
    const DataType_ _v_max;

  public:
    DataType_ drag;
    DataType_ lift;

    explicit BenchBodyForceAccumulator(bool defo, DataType_ nu, DataType_ v_max) :
      _defo(defo), _nu(nu), _v_max(v_max),
      drag(DataType_(0)), lift(DataType_(0))
    {
    }

    /// 2D variant
    template<typename T_>
    void operator()(
      const T_ omega,
      const Tiny::Vector<T_, 2, 2>& /*pt*/,
      const Tiny::Matrix<T_, 2, 1, 2, 1>& jac,
      const Tiny::Vector<T_, 2, 2>& /*val_v*/,
      const Tiny::Matrix<T_, 2, 2, 2, 2>& grad_v,
      const T_ val_p)
    {
      // compute normal and tangential
      const T_ n2 = T_(1) / Math::sqrt(jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0));
      const T_ tx = jac(0,0) * n2;
      const T_ ty = jac(1,0) * n2;
      const T_ nx = -ty;
      const T_ ny =  tx;

      /// \todo adjust this to support the deformation tensor!!!

      Tiny::Matrix<T_, 2, 2, 2, 2> nt;
      nt(0,0) = tx * nx;
      nt(0,1) = tx * ny;
      nt(1,0) = ty * nx;
      nt(1,1) = ty * ny;

      const T_ dut = Tiny::dot(nt, grad_v);
      const T_ dpf1 = _nu;
      const T_ dpf2 = (2.0 / (0.1*Math::sqr(_v_max*(2.0/3.0)))); // = 2 / (rho * U^2 * D)

      drag += DataType_(omega * dpf2 * ( dpf1 * dut * ny - val_p * nx));
      lift += DataType_(omega * dpf2 * (-dpf1 * dut * nx - val_p * ny));
    }

    /// 3D variant
    template<typename T_>
    void operator()(
      const T_ omega,
      const Tiny::Vector<T_, 3, 3>& /*pt*/,
      const Tiny::Matrix<T_, 3, 2, 3, 2>& jac,
      const Tiny::Vector<T_, 3, 3>& /*val_v*/,
      const Tiny::Matrix<T_, 3, 3, 3, 3>& grad_v,
      const T_ val_p)
    {
      // compute normal and tangential
      const T_ n2 = T_(1) / Math::sqrt(jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0));
      const T_ tx = jac(0,0) * n2;
      const T_ ty = jac(1,0) * n2;
      const T_ nx = -ty;
      const T_ ny =  tx;

      /// \todo adjust this to support the deformation tensor!!!

      Tiny::Matrix<T_, 3, 3, 3, 3> nt;
      nt.format();
      nt(0,0) = tx * nx;
      nt(0,1) = tx * ny;
      nt(1,0) = ty * nx;
      nt(1,1) = ty * ny;

      const T_ dut = Tiny::dot(nt, grad_v);
      const T_ dpf1 = _nu;
      const T_ dpf2 = (2.0 / (0.1*Math::sqr(_v_max*(4.0/9.0))* 0.41)); // = 2 / (rho * U^2 * D * H)

      drag += DataType_(omega * dpf2 * ( dpf1 * dut * ny - val_p * nx));
      lift += DataType_(omega * dpf2 * (-dpf1 * dut * nx - val_p * ny));
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(&drag, &drag, std::size_t(1), Dist::op_sum);
      comm.allreduce(&lift, &lift, std::size_t(1), Dist::op_sum);
    }
  }; // class BenchBodyForces<...,2>

  template<typename DataType_>
  class XFluxAccumulator
  {
  public:
    DataType_ flux;

    XFluxAccumulator() :
      flux(DataType_(0))
    {
    }

    template<typename T_, int d_, int d2_>
    void operator()(
      const T_ omega,
      const Tiny::Vector<T_, d_, d_>& /*pt*/,
      const Tiny::Matrix<T_, d_, d2_, d_, d2_>& /*jac*/,
      const Tiny::Vector<T_, d_, d_>& val_v,
      const Tiny::Matrix<T_, d_, d_, d_, d_>& /*grad_v*/,
      const T_ /*val_p*/)
    {
      flux += DataType_(omega * val_v[0]);
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(&flux, &flux, std::size_t(1), Dist::op_sum);
    }
  }; // class XFluxAccumulator

  /**
   * \brief Navier-Stokes System Level class
   */
  template<
    int dim_,
    typename MemType_ = Mem::Main,
    typename DataType_ = Real,
    typename IndexType_ = Index,
    typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, dim_>,
    typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, dim_, 1>,
    typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<MemType_, DataType_, IndexType_, 1, dim_>,
    typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>>
  class NavierStokesBlockedSystemLevel :
    public Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
  {
  public:
    typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, MemType_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;

    // the filtered local system matrix for Vanka
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    template<typename SpaceV_, typename Cubature_>
    void assemble_velocity_laplace_matrix(const SpaceV_& space_velo, const Cubature_& cubature, const DataType_ nu, bool defo)
    {
      Assembly::BurgersAssembler<DataType_, IndexType_, dim_> burgers_mat;
      burgers_mat.deformation = defo;
      burgers_mat.nu = nu;
      burgers_mat.beta = 0.0;
      burgers_mat.theta = 0.0;

      auto& loc_a = this->matrix_a.local();
      loc_a.format();

      // create dummy convection vector
      auto vec_c = loc_a.create_vector_l();
      vec_c.format();

      burgers_mat.assemble_matrix(loc_a, vec_c, space_velo, cubature);
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
  }; // class NavierStokesBlockedSystemLevel


  template<typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our arch types
    typedef Mem::Main MemType;
    typedef double DataType;
    typedef Index IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename DomainLevelType::SpaceVeloType SpaceVeloType;
    typedef typename DomainLevelType::SpacePresType SpacePresType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;

    // define our system level
    typedef NavierStokesBlockedSystemLevel<dim, MemType, DataType, IndexType> SystemLevelType;

    BenchmarkSummary<DataType, dim> summary;
    StopWatch watch_total_run;
    watch_total_run.start();

    /* ****************************************************************************************** */

    DataType nu = 1e-3;
    DataType v_max = DataType(dim) * DataType(0.15); // 2D: 0.3, 3D: 0.45
    DataType upsam = 0.0;       // streamline diffusion parameter
    Index max_nl_steps = 10;    // max. nonlinear solver iterations
    Index max_mg_steps = 25;    // max. multigrid iterations
    Index smooth_steps = 8;
    DataType smooth_damp = 0.7;
    DataType linsol_tol = 1E-5; // rel. tolerance for linear solver
    DataType nonlin_tol = 1E-8; // rel. tolerance for nonlinear solver
    Index gmres_k = 50;         // k of GMRES(k) coarse grid solver

    args.parse("nu", nu);
    args.parse("v-max", v_max);
    args.parse("upsam", upsam);
    args.parse("max-nl-steps", max_nl_steps);
    args.parse("max-mg-steps", max_mg_steps);
    args.parse("smooth-steps", smooth_steps);
    args.parse("smooth-damp", smooth_damp);
    args.parse("linsol-tol", linsol_tol);
    args.parse("nonlin-tol", nonlin_tol);
    args.parse("gmresk", gmres_k);

    const bool navier = (args.check("stokes") < 0);
    const bool newton = (args.check("picard") < 0);
    const bool defo = (args.check("defo") >= 0);
    const bool old_vanka = (args.check("old-vanka") >= 0);
    const bool testmode = (args.check("test-mode") >= 0);

#ifdef FEAT_HAVE_UMFPACK
    const bool umf_cgs = (domain.back_layer().comm().size() == 1);
#else
    const bool umf_cgs = false;
#endif

    {
      static constexpr std::size_t pl = 20u;
      static constexpr char pc = '.';
      comm.print("\nProblem Parameters:");
      comm.print(String("Nu").pad_back(pl, pc) + ": " + stringify(nu));
      comm.print(String("V-Max").pad_back(pl, pc) + ": " + stringify(v_max));
      comm.print(String("Upsam").pad_back(pl, pc) + ": " + stringify(upsam));
      comm.print(String("System").pad_back(pl, pc) + ": " + (navier ? "Navier-Stokes" : "Stokes"));
      comm.print(String("Tensor").pad_back(pl, pc) + ": " + (defo ? "Deformation" : "Gradient"));
      comm.print(String("NL-Solver").pad_back(pl, pc) + ": " + (newton ? "Newton" : "Picard"));
      comm.print(String("Smoothing Steps").pad_back(pl, pc) + ": " + stringify(smooth_steps));
      comm.print(String("Vanka Damping").pad_back(pl, pc) + ": " + stringify(smooth_damp));
      comm.print(String("LinSol RelTol").pad_back(pl, pc) + ": " + stringify_fp_sci(linsol_tol));
      comm.print(String("NonLin RelTol").pad_back(pl, pc) + ": " + stringify_fp_sci(nonlin_tol));
      comm.print(String("Max LinSol Iter").pad_back(pl, pc) + ": " + stringify(max_mg_steps));
      comm.print(String("Max NonLin Iter").pad_back(pl, pc) + ": " + stringify(max_nl_steps));
      comm.print(String("Vanka Type").pad_back(pl, pc) + ": " + (old_vanka ? "Old Vanka version" : "AmaVanka version"));
      if(umf_cgs)
        comm.print(String("Coarse Solver").pad_back(pl, pc) + ": UMFPACK");
      else
        comm.print(String("Coarse Solver").pad_back(pl, pc) + ": GMRES(" + stringify(gmres_k) + ")");
    }

    std::deque<Index> iters;
    if(testmode)
    {
      iters.push_back(Index(1));
      iters.push_back(Index(1));
    }
    else if(args.check("iters") > 0)
    {
      auto it = args.query("iters");
      for(const auto& s : it->second)
      {
        Index i;
        if(!s.parse(i))
        {
          comm.print("ERROR: Failed to parse iter count");
          Runtime::abort();
        }
        iters.push_back(i);
      }
    }

    const bool biters = !iters.empty();

    /* ****************************************************************************************** */

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = Index(domain.size_physical());

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("gauss-legendre:3");

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_gates(domain.at(i));

      if((i+1) < domain.size_virtual())
      {
        system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
        system_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature);
      }
    }

    if(navier)
    {
      // assemble velocity truncation operators -- we need those for the assembly of the
      // non-linear burgers operators on the coarser levels
      for (Index i(0); i < num_levels; ++i)
      {
        if(i+1 < num_levels)
          system_levels.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature, system_levels.at(i+1).get());
        else if(i+1 < domain.size_virtual())
          system_levels.at(i)->assemble_velocity_truncation(domain.at(i), domain.at(i+1), cubature);
      }
    }

    {
      auto tv = system_levels.front()->gate_velo._freqs.clone(LAFEM::CloneMode::Deep);
      tv.format(1.0);
      Index velo_dofs = Index(system_levels.front()->gate_velo.dot(tv, tv));
      Index locp_dofs = system_levels.front()->gate_pres._freqs.size();
      Index pres_dofs = Index(system_levels.front()->gate_pres.sum(DataType(locp_dofs)));
      summary.counts[Counts::velo_dofs] = velo_dofs/dim;
      summary.counts[Counts::pres_dofs] = pres_dofs;
      summary.counts[Counts::total_dofs] = velo_dofs+pres_dofs;
      summary.counts[Counts::elements] = domain.front()->get_mesh().get_num_elements();
    }

    /* ***************************************************************************************** */

    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->assemble_velo_struct(domain.at(i)->space_velo);
      system_levels.at(i)->assemble_pres_struct(domain.at(i)->space_pres);
      system_levels.at(i)->assemble_velocity_laplace_matrix(domain.at(i)->space_velo, cubature, nu, defo);
      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
      system_levels.at(i)->compile_system_matrix();
    }

    /* ***************************************************************************************** */

    // our inflow BC function
    InflowFunction<dim> inflow_func(v_max);

    // the names of the mesh parts on which to assemble
    std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

    for(Index i(0); i < num_levels; ++i)
    {
      // get our local velocity filter
      auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

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
          else if(name != "bnd:r")
          {
            // outflow
            unit_asm_noflow.add_mesh_part(*mesh_part);
          }
        }
      }

      // assemble the filters
      unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, inflow_func);
      unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo);

      // finally, compile the system filter
      system_levels.at(i)->compile_system_filter();
    }

    // finally, compile the local type-1 matrices
    for(Index i(0); i < num_levels; ++i)
    {
      system_levels.at(i)->compile_local_matrix();
    }

    Statistics::toe_assembly = stamp_ass.elapsed_now();

    // accumulate sizes
    {
      for(Index i(0); i < num_levels; ++i)
      {
        summary.bytes[Bytes::mesh] += domain.at(i)->get_mesh_node_ptr()->bytes();
        summary.bytes[Bytes::gate] += system_levels.at(i)->gate_sys.bytes();
        summary.bytes[Bytes::muxer] += system_levels.at(i)->coarse_muxer_sys.bytes();
        summary.bytes[Bytes::matrix] += system_levels.at(i)->matrix_sys.local().bytes();
        summary.bytes[Bytes::transfer] += system_levels.at(i)->transfer_sys.bytes();
        const auto& loc_a = system_levels.at(i)->matrix_sys.local().block_a();
        const auto& loc_b = system_levels.at(i)->matrix_sys.local().block_b();
        const auto& loc_d = system_levels.at(i)->matrix_sys.local().block_d();
        summary.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_a.used_elements() + loc_a.rows() + Index(1));
        summary.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_b.used_elements() + loc_b.rows() + Index(1));
        summary.bytes[Bytes::matrix_struct] += sizeof(IndexType) * std::size_t(loc_d.used_elements() + loc_d.rows() + Index(1));
        summary.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_a.template used_elements<LAFEM::Perspective::pod>());
        summary.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_b.template used_elements<LAFEM::Perspective::pod>());
        summary.bytes[Bytes::matrix_values] += sizeof(DataType) * std::size_t(loc_d.template used_elements<LAFEM::Perspective::pod>());
      }
    }

    /* ***************************************************************************************** */

    // get our global solver system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalSystemTransfer GlobalSystemTransfer;
    //typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    //typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // get our global solve matrix and filter
    GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
    GlobalSystemFilter& filter = the_system_level.filter_sys;

    // create new vectors
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_def = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_cor = the_system_level.matrix_sys.create_vector_r();

    // format the vectors
    vec_sol.format();
    vec_def.format();

    // and filter them
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_def);

    {
      // count non-zeros in a and b
      summary.counts[Counts::nnze_a] = the_system_level.matrix_sys.local().block_a().used_elements();
      summary.counts[Counts::nnze_b] = the_system_level.matrix_sys.local().block_b().used_elements();
      summary.counts[Counts::nnze_total] = the_system_level.matrix_sys.local().template used_elements<LAFEM::Perspective::pod>();
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a multigrid solver
    auto multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>>(domain.size_virtual());

    // array of Vanka pointers - this is only required to collect the memory usage
    // statistics of the Vankas, as we need the memory usage after factorisation
    std::deque<
      std::shared_ptr<
        Solver::Vanka<
          typename SystemLevelType::LocalSystemMatrix,
          typename SystemLevelType::LocalSystemFilter>>> old_vankas;
    std::deque<
      std::shared_ptr<
        Solver::AmaVanka<
          typename SystemLevelType::LocalSystemMatrix,
          typename SystemLevelType::LocalSystemFilter>>> ama_vankas;


    // push levels into multigrid
    for(std::size_t i(0); i < system_levels.size(); ++i)
    {
      SystemLevelType& lvl = *system_levels.at(i);

      // a pointer for our Schwarz-Vanka
      std::shared_ptr<Solver::SolverBase<GlobalSystemVector>> schwarz;

      if(old_vanka)
      {
        auto vanka = Solver::new_vanka(lvl.local_matrix_sys, lvl.filter_sys.local(), Solver::VankaType::block_full_add);
        old_vankas.push_back(vanka);
        schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
      }
      else
      {
        auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
        ama_vankas.push_back(vanka);
        schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);
      }

      if((i+1) < domain.size_virtual())
      {
        // create richardson smoother
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);
        smoother->set_min_iter(smooth_steps);
        smoother->set_max_iter(smooth_steps);
        smoother->skip_defect_calc(true); // skip defect calculation
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
#ifdef FEAT_HAVE_UMFPACK
      else if(umf_cgs)
      {
        // create UMFPACK coarse grid solver
        auto umfpack = Solver::new_generic_umfpack(lvl.local_matrix_sys);
        auto cgsolver = Solver::new_schwarz_precond(umfpack, lvl.filter_sys);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
#endif //  FEAT_HAVE_UMFPACK
      else
      {
        // create coarse grid solver
        auto cgsolver = Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, gmres_k, 0.0, schwarz);
        cgsolver->set_max_iter(1000);
        cgsolver->set_tol_rel(1e-3);
        //cgsolver->set_plot_mode(Solver::PlotMode::summary);

        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
    }

    // create our multigrid solver
    auto multigrid = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

    // create our solver
    auto solver = Solver::new_richardson(matrix, filter, 1.0, multigrid);
    //auto solver = Solver::new_fgmres(matrix, filter, 20, 0.0, multigrid);
    //auto solver = Solver::new_bicgstab(matrix, filter, multigrid);

    solver->set_plot_name("Multigrid");

    solver->set_max_iter(max_mg_steps);
    solver->set_tol_rel(linsol_tol);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a few watches
    StopWatch watch_stokes_solve;
    StopWatch watch_nonlin_loop;
    StopWatch watch_nonlin_def_asm;
    StopWatch watch_nonlin_mat_asm;
    StopWatch watch_nonlin_solver_init;
    StopWatch watch_nonlin_solver_apply;

    comm.print("\nSolving Stokes system...");

    // initialise solver
    watch_nonlin_solver_init.start();
    multigrid_hierarchy->init();
    solver->init();
    watch_nonlin_solver_init.stop();

    // accumulate vanka sizes
    if(old_vanka)
    {
      summary.counts[Counts::vanka_data] = Index(old_vankas.front()->data_size());
      summary.bytes[Bytes::vanka] = 0ull;
      for(auto& v : old_vankas)
        summary.bytes[Bytes::vanka] += v->bytes();
    }
    else
    {
      summary.counts[Counts::vanka_data] = Index(ama_vankas.front()->data_size());
      summary.bytes[Bytes::vanka] = 0ull;
      for(auto& v : ama_vankas)
        summary.bytes[Bytes::vanka] += v->bytes();
    }

    // solve Stokes system
    solver->set_plot_mode(Solver::PlotMode::iter);
    if(biters)
    {
      solver->set_min_iter(iters.front());
      solver->set_max_iter(iters.front());
      iters.pop_front();
    }
    watch_stokes_solve.start();
    Solver::Status stokes_status = Solver::solve(*solver, vec_sol, vec_def, matrix, filter);
    watch_stokes_solve.stop();

    solver->set_plot_mode(Solver::PlotMode::none);

    summary.counts[Counts::linsol_iter] = solver->get_num_iter();

    // release solvers
    solver->done_numeric();
    multigrid_hierarchy->done_numeric();

    FEAT::Statistics::compress_solver_expressions();

    if(!Solver::status_success(stokes_status))
    {
      comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
      if(testmode)
        comm.print("Test-Mode: FAILED");
      return;
    }

    if(navier)
    {
      comm.print("\nSolving Navier-Stokes system...");

      // setup burgers assembler for matrix
      Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_mat;
      burgers_mat.deformation = defo;
      burgers_mat.nu = nu;
      burgers_mat.beta = DataType(1);
      burgers_mat.frechet_beta = DataType(newton ? 1 : 0);
      burgers_mat.sd_delta = upsam;
      burgers_mat.sd_nu = nu;

      // setup burgers assembler for defect vector
      Assembly::BurgersAssembler<DataType, IndexType, dim> burgers_def;
      burgers_def.deformation = defo;
      burgers_def.nu = nu;
      burgers_def.beta = DataType(1);

      // assemble non-linear defect
      DataType def_nl_init = DataType(0);

      watch_nonlin_loop.start();

      // nonlinear loop
      for(Index nl_step(0); nl_step <= max_nl_steps; ++nl_step)
      {
        summary.counts[Counts::nonlin_iter] = nl_step;

        // assemble nonlinear defect vector
        watch_nonlin_def_asm.start();
        vec_def.format();
        // assemble burgers operator defect
        burgers_def.assemble_vector(vec_def.local().template at<0>(), vec_sol.local().template at<0>(),
          vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature, -1.0);
        // compute remainder of defect vector
        the_system_level.matrix_sys.local().block_b().apply(
          vec_def.local().template at<0>(), vec_sol.local().template at<1>(), vec_def.local().template at<0>(), -1.0);
        the_system_level.matrix_sys.local().block_d().apply(
          vec_def.local().template at<1>(), vec_sol.local().template at<0>(), vec_def.local().template at<1>(), -1.0);
        // sync and filter
        vec_def.sync_0();
        filter.filter_def(vec_def);
        watch_nonlin_def_asm.stop();

        // compute defect norm
        DataType def_nl = vec_def.norm2();

        if(nl_step == Index(0))
          def_nl_init = def_nl;

        String line = (newton ? "Newton: " : "Picard: ");
        line += stringify(nl_step).pad_front(2) + ": ";
        line += stringify_fp_sci(def_nl, 6) + " / " + stringify_fp_sci(def_nl/def_nl_init);

        if(def_nl > def_nl_init * 1E+3)
        {
          comm.print(line);
          comm.print("\nERROR: NONLINEAR SOLVER DIVERGED !!!\n");
          if(testmode)
            comm.print("Test-Mode: FAILED");
          return;
        }
        else if(biters)
        {
          if(iters.empty())
          {
            comm.print(line);
            break;
          }
        }
        else if((nl_step >= max_nl_steps) || (def_nl < def_nl_init * nonlin_tol))
        {
          comm.print(line);
          break;
        }

        // assemble burgers matrices on all levels
        watch_nonlin_mat_asm.start();
        {
          // get a clone of the global velocity vector
          typename SystemLevelType::GlobalVeloVector vec_conv(
            &the_system_level.gate_velo, vec_sol.local().template at<0>().clone());

          // set velocity norm for streamline diffusion (if enabled)
          if(Math::abs(upsam) > DataType(0))
          {
            // set norm by convection vector; we can use this for all levels
            burgers_mat.set_sd_v_norm(vec_conv);
          }

          // loop over all system levels
          for(std::size_t i(0); i < system_levels.size(); ++i)
          {
            // assemble our system matrix
            auto& loc_mat_a = system_levels.at(i)->matrix_sys.local().block_a();
            loc_mat_a.format();
            burgers_mat.assemble_matrix(loc_mat_a, vec_conv.local(), domain.at(i)->space_velo, cubature);
            system_levels.at(i)->compile_local_matrix();

            // restrict our convection vector
            if((i+1) >= domain.size_virtual())
              break;

            // does this process have another system level?
            if((i+1) < system_levels.size())
            {
              // create a coarse mesh velocity vector
              auto vec_crs = system_levels.at(i+1)->matrix_a.create_vector_l();

              // truncate fine mesh velocity vector
              system_levels.at(i)->transfer_velo.trunc(vec_conv, vec_crs);

              // the coarse vector is our next convection vector
              vec_conv = std::move(vec_crs);
            }
            else
            {
              // this process is a child, so send truncation to parent
              system_levels.at(i)->transfer_velo.trunc_send(vec_conv);
            }
          }
        }
        watch_nonlin_mat_asm.stop();

        // initialise linear solver
        watch_nonlin_solver_init.start();
        multigrid_hierarchy->init_numeric();
        solver->init_numeric();
        watch_nonlin_solver_init.stop();
        if(biters)
        {
          solver->set_min_iter(iters.front());
          solver->set_max_iter(iters.front());
          iters.pop_front();
        }

        // solve linear system
        watch_nonlin_solver_apply.start();
        Solver::Status status = solver->apply(vec_cor, vec_def);
        watch_nonlin_solver_apply.stop();

        summary.counts[Counts::linsol_iter] += solver->get_num_iter();

        line += String(" | ") + stringify(solver->get_num_iter()).pad_front(3) + ": "
          + stringify_fp_sci(solver->get_def_final() / solver->get_def_initial(), 6);
        comm.print(line);

        // release linear solver
        solver->done_numeric();
        multigrid_hierarchy->done_numeric();

        if(!Solver::status_success(status))
        {
          comm.print("\nERROR: LINEAR SOLVER BREAKDOWN\n");
          if(testmode)
            comm.print("Test-Mode: FAILED");
          return;
        }

        // update solution
        vec_sol.axpy(vec_cor, vec_sol, 1.0);

        FEAT::Statistics::compress_solver_expressions();
        // next non-linear iteration
      }

      watch_nonlin_loop.stop();

      // end of Navier-Stokes solve
    }

    // release solver
    solver->done_symbolic();
    multigrid_hierarchy->done_symbolic();

    // save timings
    summary.times[Times::stokes_solve] = watch_stokes_solve.elapsed();
    summary.times[Times::nonlin_total] = watch_nonlin_loop.elapsed();
    summary.times[Times::nonlin_asm_def] = watch_nonlin_def_asm.elapsed();
    summary.times[Times::nonlin_asm_mat] = watch_nonlin_mat_asm.elapsed();
    summary.times[Times::linsol_init] = watch_nonlin_solver_init.elapsed();
    summary.times[Times::linsol_apply] = watch_nonlin_solver_apply.elapsed();

    // get multigrid timings
    summary.times[Times::mg_defect] = multigrid_hierarchy->get_time_defect();
    summary.times[Times::mg_smooth] = multigrid_hierarchy->get_time_smooth();
    summary.times[Times::mg_coarse] = multigrid_hierarchy->get_time_coarse();
    summary.times[Times::mg_transfer] = multigrid_hierarchy->get_time_transfer();

    // accumulate vanka timings
    if(old_vanka)
    {
      for(auto& v : old_vankas)
      {
        summary.times[Times::vanka_init_sym] += v->time_init_symbolic();
        summary.times[Times::vanka_init_num] += v->time_init_numeric();
        summary.times[Times::vanka_apply] += v->time_apply();
      }
    }
    else
    {
      for(auto& v : ama_vankas)
      {
        summary.times[Times::vanka_init_sym] += v->time_init_symbolic();
        summary.times[Times::vanka_init_num] += v->time_init_numeric();
        summary.times[Times::vanka_apply] += v->time_apply();
      }
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // compute drag & lift coefficients
    {
      Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> body_asm(the_domain_level.trafo);

      auto* mesh_part = the_domain_level.get_mesh_node()->find_mesh_part("bnd:c");
      if(mesh_part != nullptr)
      {
        body_asm.add_mesh_part(*mesh_part);
      }

      body_asm.compile_facets();

      BenchBodyForceAccumulator<DataType> body_force_accum(defo, nu, v_max);
      body_asm.assemble_flow_accum(
        body_force_accum,
        vec_sol.local().template at<0>(),
        vec_sol.local().template at<1>(),
        the_domain_level.space_velo,
        the_domain_level.space_pres,
        cubature);
      body_force_accum.sync(comm);

      // reference values
      const DataType ref_drag = (dim == 2 ? 5.57953523384 : 6.18533);
      const DataType ref_lift = (dim == 2 ? 0.010618948146 : 0.009401);

      summary.drag_coeff = body_force_accum.drag;
      summary.lift_coeff = body_force_accum.lift;
      summary.drag_err = Math::abs((body_force_accum.drag - ref_drag) / ref_drag);
      summary.lift_err = Math::abs((body_force_accum.lift - ref_lift) / ref_lift);
    }

    // compute pressure difference
    {
      typedef Trafo::InverseMapping<typename SpacePresType::TrafoType, DataType> InvMappingType;
      InvMappingType inv_mapping(the_domain_level.trafo);

      // reference pressure points
      typename InvMappingType::ImagePointType v_a, v_e;
      if(dim == 2)
      {
        v_a[0] = 0.15;
        v_e[0] = 0.25;
        v_a[1] = v_e[1] = 0.2;
      }
      else
      {
        v_a[0] = 0.45;
        v_e[0] = 0.55;
        v_a[1] = v_e[1] = 0.2;
        v_a[2] = v_e[2] = 0.205;
      }

      // unmap points
      auto iv_a = inv_mapping.unmap_point(v_a, true);
      auto iv_e = inv_mapping.unmap_point(v_e, true);

      // evaluate pressure
      auto pval_a = Assembly::DiscreteEvaluator::eval_fe_function(
        iv_a, vec_sol.local().template at<1>(), the_domain_level.space_pres);
      auto pval_e = Assembly::DiscreteEvaluator::eval_fe_function(
        iv_e, vec_sol.local().template at<1>(), the_domain_level.space_pres);

      // compute pressure mean
      const auto p_a = pval_a.mean_value_dist(comm);
      const auto p_e = pval_e.mean_value_dist(comm);
      const auto d_p = p_a - p_e;

      // reference pressure values
      const DataType ref_d_p = (dim == 2 ? 0.11752016697 : 0.170826996);

      // compute error to reference values
      summary.pres_diff = d_p;
      summary.pres_err = Math::abs((d_p - ref_d_p) / ref_d_p);
    }

    // compute flow through upper region
    {
      Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> flux_asm(the_domain_level.trafo);

      auto* mesh_part = the_domain_level.get_mesh_node()->find_mesh_part("inner:u");
      if(mesh_part != nullptr)
      {
        flux_asm.add_mesh_part(*mesh_part);
      }

      flux_asm.compile_facets(false);

      XFluxAccumulator<DataType> flux_accum;
      flux_asm.assemble_flow_accum(
        flux_accum,
        vec_sol.local().template at<0>(),
        vec_sol.local().template at<1>(),
        the_domain_level.space_velo,
        the_domain_level.space_pres,
        cubature);
      flux_accum.sync(comm);

      summary.flux_upper = flux_accum.flux / DataType(2);
    }

    // compute flow through lower region
    {
      Assembly::TraceAssembler<typename SpaceVeloType::TrafoType> flux_asm(the_domain_level.trafo);

      auto* mesh_part = the_domain_level.get_mesh_node()->find_mesh_part("inner:l");
      if(mesh_part != nullptr)
      {
        flux_asm.add_mesh_part(*mesh_part);
      }

      flux_asm.compile_facets(false);

      XFluxAccumulator<DataType> flux_accum;
      flux_asm.assemble_flow_accum(
        flux_accum,
        vec_sol.local().template at<0>(),
        vec_sol.local().template at<1>(),
        the_domain_level.space_velo,
        the_domain_level.space_pres,
        cubature);
      flux_accum.sync(comm);

      summary.flux_lower = flux_accum.flux / DataType(2);
    }

    {
      Assembly::VelocityInfo<DataType, dim> vi = Assembly::VelocityAnalyser::compute(
        vec_sol.local().template at<0>(), the_domain_level.space_velo, cubature);
      vi.synchronise(comm);

      summary.velo_info = vi;
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    //*
    if(args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name, vtk_name2;
      int npars = args.parse("vtk", vtk_name, vtk_name2);
      if(npars < 1)
      {
        vtk_name = String("nvs-cc") + stringify(dim) + "d-q2p1";
        vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
        vtk_name += "-n" + stringify(comm.size());
      }

      {
        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

        // project velocity and pressure
        exporter.add_vertex_vector("v", vec_sol.local().template at<0>());

        // project pressure
        Cubature::DynamicFactory cub("gauss-legendre:2");
        LAFEM::DenseVector<Mem::Main, double, Index> vtx_p;
        Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), the_domain_level.space_pres, cub);

        // write pressure
        exporter.add_cell_scalar("pressure", vtx_p.elements());

        // finally, write the VTK file
        exporter.write(vtk_name, comm);
      }

      // write refined VTK?
      if(npars > 1)
      {
        // refine mesh
        Geometry::StandardRefinery<MeshType> refinery(the_domain_level.get_mesh());
        MeshType mesh(refinery);

        // Create a VTK exporter for our mesh
        Geometry::ExportVTK<MeshType> exporter(mesh);

        // project velocity and pressure
        exporter.add_vertex_vector("v", vec_sol.local().template at<0>());

        // finally, write the VTK file
        exporter.write(vtk_name2, comm);
      }
    }

    watch_total_run.stop();
    summary.times[Times::total_run] = watch_total_run.elapsed();

    {
      MemoryUsage mi;
      summary.bytes[Bytes::peak_p] = mi.get_peak_physical();
      summary.bytes[Bytes::peak_v] = mi.get_peak_virtual();
    }

    summary.sync(comm);

    comm.print(String("\n") + String(80, '=') + "\n");
    comm.print(summary.format(20));

    // print multigrid timings
    if(comm.rank() == 0)
    {
      comm.print("Multigrid Timings:");
      comm.print("              Defect /   Smoother /   Transfer /     Coarse");
      comm.print("Overall : " +
          stringify_fp_fix(multigrid_hierarchy->get_time_defect(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_smooth(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_transfer(), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_coarse(), 3, 10));
      for(int i(0); i < int(multigrid_hierarchy->size_physical()); ++i)
      {
        comm.print("Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(multigrid_hierarchy->get_time_defect(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_smooth(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_transfer(i), 3, 10) + " / " +
          stringify_fp_fix(multigrid_hierarchy->get_time_coarse(i), 3, 10));
      }
    }

    comm.print("\n");
    comm.print(FEAT::Statistics::get_formatted_flops(summary.times[Times::stokes_solve] + summary.times[Times::linsol_apply], comm.size()));
    comm.print(FEAT::Statistics::get_formatted_times(summary.times[Times::stokes_solve] + summary.times[Times::linsol_apply]));
    comm.print(FEAT::Statistics::get_formatted_solver_internals("default"));
    comm.print("\n");
    comm.print(FEAT::Statistics::get_formatted_solver_tree("default").trim());

    if(testmode)
      comm.print("\nTest-Mode: PASSED");
  }

  template<int dim_>
  void run_dim(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
    // define our mesh type
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;

    // create a time-stamp
    TimeStamp time_stamp;

    // create our domain control
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(mesh_reader);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now(TimeFormat::s_m));
  }

  void main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("level");
    args.support("vtk");
    args.support("mesh");
    args.support("nu");
    args.support("upsam");
    args.support("max-nl-steps");
    args.support("max-mg-steps");
    args.support("smooth-steps");
    args.support("smooth-damp");
    args.support("linsol-tol");
    args.support("nonlin-tol");
    args.support("picard");
    args.support("v-max");
    args.support("defo");
    args.support("stokes");
    args.support("iters");
    args.support("gmresk");
    args.support("test-mode");
    args.support("old-vanka");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      // abort
      FEAT::Runtime::abort();
    }

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
      FEAT::Runtime::abort();
    }
    if(args.check("level") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
      FEAT::Runtime::abort();
    }

    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run 2D or 3D ?
    if(mesh_type == "conformal:hypercube:2:2")
      run_dim<2>(args, comm, mesh_reader);
    else if(mesh_type == "conformal:hypercube:3:3")
      run_dim<3>(args, comm, mesh_reader);
    else
    {
      comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
      FEAT::Runtime::abort();
    }
  }
} // namespace NvSCCNDQ2P1dc

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialise(argc, argv);
  try
  {
    NvSCCNDQ2P1dc::main(argc, argv);
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
  return FEAT::Runtime::finalise();
}
