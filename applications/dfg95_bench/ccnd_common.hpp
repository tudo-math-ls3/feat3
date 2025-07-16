// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/assembly/burgers_assembly_job.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/assembly/trace_assembler_basic_jobs.hpp>
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
#include <control/checkpoint_control.hpp>
#include <control/statistics.hpp>

#include <vector>

namespace DFG95
{
  using namespace FEAT;

  // helper function to parse arguments
  template<typename T_>
  T_ parse(SimpleArgParser& args, const String& name, T_ default_value)
  {
    args.parse(name, default_value);
    return default_value;
  }

  // our one and only index type
  typedef Index IndexType;

  // depending on whether FEAT_CCND_USE_QUADMATH is defined,
  // we use use quad precision or double precision
#ifdef FEAT_CCND_USE_QUADMATH
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

  // template for the steady-state inflow function used in bench1 and bench2
  template<int dim_>
  class SteadyInflowFunction;

  // 2D steady inflow function used in bench1 and bench2
  template<>
  class SteadyInflowFunction<2> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 2;
    typedef Analytic::Image::Vector<2> ImageType;
    static constexpr bool can_value = true;

  protected:
    DataType _vmax;

  public:
    explicit SteadyInflowFunction(DataType vmax) :
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
      explicit Evaluator(const SteadyInflowFunction& function) :
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
  }; // class SteadyInflowFunction<2>

  // 3D steady inflow function used in bench1 and bench2
  template<>
  class SteadyInflowFunction<3> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 3;
    typedef Analytic::Image::Vector<3> ImageType;
    static constexpr bool can_value = true;

  protected:
    DataType _vmax;

  public:
    explicit SteadyInflowFunction(DataType vmax) :
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
      explicit Evaluator(const SteadyInflowFunction& function) :
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
  }; // class SteadyInflowFunction

  // template for the unsteady inflow function used in bench3
  template<int dim_>
  class UnsteadyInflowFunction;

  // 2D steady inflow function used in bench3
  template<>
  class UnsteadyInflowFunction<2> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 2;
    typedef Analytic::Image::Vector<2> ImageType;
    static constexpr bool can_value = true;

  protected:
    DataType _vmax;
    DataType _time;

  public:
    explicit UnsteadyInflowFunction(DataType vmax, DataType t) :
      _vmax(vmax),
      _time(t)
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

      const DataType _vmax, _t, _pi, _d, _scale;

    public:
      explicit Evaluator(const UnsteadyInflowFunction& function) :
        _vmax(function._vmax),
        _t(function._time),
        _pi(Math::pi<DataType>()),
        _d(DataType(0.41)),
        _scale(Math::sin(_pi*_t/DataType(8))/(_d*_d))
      {
      }

      ValueType value(const PointType& point)
      {
        const DataType y = point[1];

        ValueType val;
        val[0] = (_vmax * DataType(4) * y * (_d - y)) * _scale;
        val[1] = DataType(0);
        return val;
      }
    };
  }; // class UnsteadyInflowFunction<2>

  // 3D steady inflow function used in bench3
  template<>
  class UnsteadyInflowFunction<3> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 3;
    typedef Analytic::Image::Vector<3> ImageType;
    static constexpr bool can_value = true;

  protected:
    DataType _vmax;
    DataType _time;

  public:
    explicit UnsteadyInflowFunction(DataType vmax, DataType t) :
      _vmax(vmax),
      _time(t)
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

      const DataType _vmax, _t, _pi, _d, _scale;

    public:
      explicit Evaluator(const UnsteadyInflowFunction& function) :
        _vmax(function._vmax),
        _t(function._time),
        _pi(Math::pi<DataType>()),
        _d(DataType(0.41)),
        _scale(Math::sin(_pi*_t/DataType(8))/(_d*_d*_d*_d))
      {
      }

      ValueType value(const PointType& point)
      {
        const DataType y = point[1];
        const DataType z = point[2];

        ValueType val;
        val[0] = (_vmax * DataType(16) * y * (_d - y) * z * (_d - z)) * _scale;
        val[1] = val[2] = DataType(0);
        return val;
      }
    };
  }; // class UnsteadyInflowFunction

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
    static constexpr std::size_t count = 17u;
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

  // class for collecting benchmark postprocessing summary data
  template<typename DT_, int dim_>
  class BenchmarkSummary
  {
  public:
    // post-processing data

    // drag/lift forces by line integration
    DT_ drag_coeff_line, lift_coeff_line, drag_err_line, lift_err_line;

    // drag/lift forces by volume integration
    DT_ drag_coeff_vol, lift_coeff_vol, side_coeff_vol, drag_err_vol, lift_err_vol, side_err_vol;

    // pressure values
    DT_ pres_diff, pres_err;

    // flow through upper/lower region
    DT_ flux_upper, flux_lower;

    // velocity field information
    Assembly::VelocityInfo<DT_, dim_> velo_info;

    static constexpr std::size_t padlen = std::size_t(30);

    BenchmarkSummary() :
      drag_coeff_line(0.0), lift_coeff_line(0.0), drag_err_line(0.0), lift_err_line(0.0),
      drag_coeff_vol(0.0), lift_coeff_vol(0.0), side_coeff_vol(0.0), drag_err_vol(0.0), lift_err_vol(0.0), side_err_vol(0.0),
      pres_diff(0.0), pres_err(0.0),
      flux_upper(0.0), flux_lower(0.0)
    {
    }

    String format(int prec = fp_num_digs+5) const
    {
      String s;
      const char pc = '.';
      // append coefficients and velocity info
      s += "Solution Analysis:\n";
      s += String("Drag Coefficient (Line)").pad_back(padlen, pc) + ": " + stringify_fp_fix(drag_coeff_line, prec)
        + "   [ Error: " + stringify_fp_sci(drag_err_line, prec) + " ]\n";
      s += String("Lift Coefficient (Line)").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff_line, prec)
        + "   [ Error: " + stringify_fp_sci(lift_err_line, prec) + " ]\n";
      s += String("Drag Coefficient (Vol)").pad_back(padlen, pc) + ": " + stringify_fp_fix(drag_coeff_vol, prec)
        + "   [ Error: " + stringify_fp_sci(drag_err_vol, prec) + " ]\n";
      s += String("Lift Coefficient (Vol)").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff_vol, prec)
        + "   [ Error: " + stringify_fp_sci(lift_err_vol, prec) + " ]\n";
      s += String("Side Coefficient (Vol)").pad_back(padlen, pc) + ": " + stringify_fp_fix(side_coeff_vol, prec)
        + "   [ Error: " + stringify_fp_sci(side_err_vol, prec) + " ]\n";
      s += String("Pressure Difference").pad_back(padlen, pc) + ": " + stringify_fp_fix(pres_diff, prec)
        + "   [ Error: " + stringify_fp_sci(pres_err, prec) + " ]\n";
      s += String("Upper Flux").pad_back(padlen, pc) + ": " + stringify_fp_fix(flux_upper, prec) + "\n";
      s += String("Lower Flux").pad_back(padlen, pc) + ": " + stringify_fp_fix(flux_lower, prec) + "\n";
      s += velo_info.format_string(prec, padlen, pc);

      return s;
    }

    String format_compact(String prefix, int prec = fp_num_digs) const
    {
      String s;
      s += prefix + "DCL/LCL/DC/LC/SC/PD: ";
      s += stringify_fp_fix(drag_coeff_line, prec) + " " + stringify_fp_fix(lift_coeff_line, prec) + " ";
      s += stringify_fp_fix(drag_coeff_vol, prec) + " " + stringify_fp_fix(lift_coeff_vol, prec) + " " + stringify_fp_fix(side_coeff_vol, prec) + " ";
      s += stringify_fp_fix(pres_diff, prec) + "\n";
      s += prefix + "FX/H0/H1: ";
      s += stringify_fp_fix(flux_upper, prec) + " " + stringify_fp_fix(flux_lower, prec) + " ";
      s += stringify_fp_fix(velo_info.norm_h0, prec) + " " + stringify_fp_fix(velo_info.norm_h1, prec) + "\n";
      s += prefix + "VC/DV...: ";
      s += stringify_fp_fix(velo_info.vorticity, prec) + " " + stringify_fp_fix(velo_info.divergence, prec);
      return s;
    }
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

  // trace assembly job for surface body forces computation; this class extends the assembly that is also included
  // in Assembly::TraceAssemblyStokesBodyForceAssemblyJob by the formulae given in the original DFG95 paper, which
  // only work for 2.5D obstacles, whereas the 'new' formulae from the Giles paper also work for 'true' 3D obstacles,
  // however, the resulting coefficients by the two formulae are slightly different, so we print out the coefficients
  // according to both the old and the new fomulae in this application
  template<typename VectorVelo_, typename VectorPres_, typename SpaceVelo_, typename SpacePres_>
  class DFG95SurfaceBodyForceAssemblyJob
  {
  public:
    typedef typename VectorVelo_::DataType DataType;
    typedef typename VectorVelo_::ValueType VeloValueType;
    typedef typename VectorPres_::ValueType PresValueType;

    static constexpr TrafoTags trafo_config = TrafoTags::none;
    static constexpr SpaceTags space_velo_config = SpaceTags::value | SpaceTags::grad;
    static constexpr SpaceTags space_pres_config = SpaceTags::value;
    static constexpr TrafoTags facet_trafo_config = TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::normal;

    static constexpr int dim = SpaceVelo_::shape_dim;

    typedef Tiny::Matrix<DataType, 3, 2> RawNewForces;
    typedef Tiny::Matrix<DataType, 2, 2> RawOldForces;

    class Task :
      public Assembly::TraceAssemblyStokesVectorAnalysisTaskCRTP<Task, VectorVelo_, VectorPres_, SpaceVelo_, SpacePres_, trafo_config, facet_trafo_config, space_velo_config, space_pres_config>
    {
    public:
      /// our base-class typedef
      typedef Assembly::TraceAssemblyStokesVectorAnalysisTaskCRTP<Task, VectorVelo_, VectorPres_, SpaceVelo_, SpacePres_, trafo_config, facet_trafo_config, space_velo_config, space_pres_config> BaseClass;

      /// our assembly traits
      typedef typename BaseClass::AsmTraits AsmTraits;

      // we do not support pairwise assembly
      static constexpr bool assemble_pairwise = false;
      /// this task needs to scatter
      static constexpr bool need_scatter = false;
      /// this task has no combine
      static constexpr bool need_combine = true;

      using typename BaseClass::TrafoEvalData;
      using typename BaseClass::TrafoFacetEvalData;
      using typename BaseClass::VeloData;
      using typename BaseClass::PresData;

    protected:
      RawOldForces& job_raw_old_forces;
      RawNewForces& job_raw_new_forces;
      RawOldForces raw_old_forces;
      RawNewForces raw_new_forces;

    public:
      explicit Task(DFG95SurfaceBodyForceAssemblyJob& job) :
        BaseClass(job.vector_velo, job.vector_pres, job.space_velo, job.space_pres, job.cubature_factory),
        job_raw_old_forces(job.raw_old_forces),
        job_raw_new_forces(job.raw_new_forces),
        raw_old_forces(),
        raw_new_forces()
      {
        raw_old_forces.format();
        raw_new_forces.format();
      }

      void eval(TrafoFacetEvalData& tau_f, TrafoEvalData& DOXY(tau), DataType weight, const VeloData& velo, const PresData& pres)
      {
        this->_eval(weight, tau_f.jac_mat, tau_f.normal, velo.grad, pres.value);
      }

      void scatter()
      {
        // nothing to do here
      }

      void combine()
      {
        job_raw_old_forces += raw_old_forces;
        job_raw_new_forces += raw_new_forces;
      }

    protected:
      /// 2D version
      void _eval(const DataType omega, const Tiny::Matrix<DataType, 2, 1, 2, 1>& jac, const Tiny::Vector<DataType, 2, 2>& n,
        const Tiny::Matrix<DataType, 2, 2, 2, 2>& grad_v, const DataType val_p)
      {
        // 'old' implementation according to original DFG95 paper
        const DataType n2 = DataType(1) / Math::sqrt(jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0));
        const DataType tx = jac(0,0) * n2;
        const DataType ty = jac(1,0) * n2;
        const DataType nx = -ty;
        const DataType ny =  tx;

        Tiny::Matrix<DataType, 2, 2, 2, 2> nt;
        nt(0,0) = tx * nx;
        nt(0,1) = tx * ny;
        nt(1,0) = ty * nx;
        nt(1,1) = ty * ny;

        const DataType dut = Tiny::dot(nt, grad_v);

        raw_old_forces(0, 0) += omega * dut * ny;
        raw_old_forces(1, 0) -= omega * dut * nx;
        raw_old_forces(0, 1) += omega * val_p * nx;
        raw_old_forces(1, 1) += omega * val_p * ny;

        // 'new' implementation according to Giles paper
        raw_new_forces(0, 0) -= omega * (DataType(2) * grad_v(0,0) * n[0] + (grad_v(0, 1) + grad_v(1, 0)) * n[1]);
        raw_new_forces(1, 0) -= omega * (DataType(2) * grad_v(1,1) * n[1] + (grad_v(1, 0) + grad_v(0, 1)) * n[0]);
        raw_new_forces(0, 1) -= omega * val_p * n[0];
        raw_new_forces(1, 1) -= omega * val_p * n[1];
      }

      /// 3D version
      void _eval(const DataType omega, const Tiny::Matrix<DataType, 3, 2, 3, 2>& jac, const Tiny::Vector<DataType, 3, 3>& n,
        const Tiny::Matrix<DataType, 3, 3, 3, 3>& grad_v, const DataType val_p)
      {
        // 'old' implementation according to original DFG95 paper
        const DataType n2 = DataType(1) / Math::sqrt(jac(0,0)*jac(0,0) + jac(1,0)*jac(1,0));
        const DataType tx = jac(0,0) * n2;
        const DataType ty = jac(1,0) * n2;
        const DataType nx = -ty;
        const DataType ny =  tx;

        Tiny::Matrix<DataType, 3, 3, 3, 3> nt;
        nt.format();
        nt(0,0) = tx * nx;
        nt(0,1) = tx * ny;
        nt(1,0) = ty * nx;
        nt(1,1) = ty * ny;

        const DataType dut = Tiny::dot(nt, grad_v);

        raw_old_forces(0, 0) += omega * dut * ny;
        raw_old_forces(1, 0) -= omega * dut * nx;
        raw_old_forces(0, 1) += omega * val_p * nx;
        raw_old_forces(1, 1) += omega * val_p * ny;

        // 'new' implementation according to Giles paper
        raw_new_forces(0, 0) -= omega * (DataType(2) * grad_v(0,0) * n[0] + (grad_v(0, 1) + grad_v(1, 0)) * n[1] + (grad_v(0, 2) + grad_v(2, 0)) * n[2]);
        raw_new_forces(1, 0) -= omega * (DataType(2) * grad_v(1,1) * n[1] + (grad_v(1, 2) + grad_v(2, 1)) * n[2] + (grad_v(1, 0) + grad_v(0, 1)) * n[0]);
        raw_new_forces(2, 0) -= omega * (DataType(2) * grad_v(2,2) * n[2] + (grad_v(2, 0) + grad_v(0, 2)) * n[0] + (grad_v(2, 1) + grad_v(1, 2)) * n[1]);
        raw_new_forces(0, 1) -= omega * val_p * n[0];
        raw_new_forces(1, 1) -= omega * val_p * n[1];
        raw_new_forces(2, 1) -= omega * val_p * n[2];
      }
    }; // class Task

  protected:
    const VectorVelo_& vector_velo;
    const VectorPres_& vector_pres;
    const SpaceVelo_& space_velo;
    const SpacePres_& space_pres;
    Cubature::DynamicFactory cubature_factory;
    RawOldForces raw_old_forces;
    RawNewForces raw_new_forces;

  public:
    explicit DFG95SurfaceBodyForceAssemblyJob(const VectorVelo_& vector_velo_, const VectorPres_& vector_pres_,
      const SpaceVelo_& space_velo_, const SpacePres_& space_pres_, String cubature_):
      vector_velo(vector_velo_),
      vector_pres(vector_pres_),
      space_velo(space_velo_),
      space_pres(space_pres_),
      cubature_factory(cubature_),
      raw_old_forces(),
      raw_new_forces()
    {
      raw_old_forces.format();
      raw_new_forces.format();
    }

    /// Returns the raw drag force coefficient for a given viscosity parameter (formula as in DFG95 paper)

    DataType drag_old(DataType nu) const
    {
      return nu * raw_old_forces(0, 0) - raw_old_forces(0, 1);
    }

    /// Returns the raw lift  force coefficient for a given viscosity parameter (formula as in DFG95 paper)
    DataType lift_old(DataType nu) const
    {
      return nu * raw_old_forces(1, 0) - raw_old_forces(1, 1);
    }

    /// Returns the raw drag force coefficient for a given viscosity parameter
    DataType drag_new(DataType nu) const
    {
      return nu * raw_new_forces(0, 0) - raw_new_forces(0, 1);
    }

    /// Returns the raw lift  force coefficient for a given viscosity parameter
    DataType lift_new(DataType nu) const
    {
      return nu * raw_new_forces(1, 0) - raw_new_forces(1, 1);
    }

    /// Returns the raw side force coefficient for a given viscosity parameter
    DataType side_new(DataType nu) const
    {
      return nu * raw_new_forces(2, 0) - raw_new_forces(2, 1);
    }

    /// Synchronizes the forces over a communicator
    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(&raw_old_forces.v[0][0], &raw_old_forces.v[0][0], std::size_t(4), Dist::op_sum);
      comm.allreduce(&raw_new_forces.v[0][0], &raw_new_forces.v[0][0], std::size_t(6), Dist::op_sum);
    }
  }; // class DFG95SurfaceBodyForceAssemblyJob<...>

  // computes the body forces by the volumetric 'defect vector' approach
  template<typename DT_, typename IT_, int dim_>
  void assemble_bdforces_vol(
    Tiny::Vector<DT_, 3>& forces,
    const LAFEM::DenseVectorBlocked<DT_, IT_, dim_>& vec_def_v,
    const LAFEM::DenseVector<DT_, IT_>& vec_char)
  {
    static_assert(dim_ <= 3, "invalid forces size");

    Tiny::Vector<DT_, dim_, 3> frc;
    frc.format();

    XASSERT(vec_def_v.size() == vec_char.size());

    // get the vector arrays
    const Index n = vec_def_v.size();
    const auto* vdef = vec_def_v.elements();
    const auto* vchr = vec_char.elements();
    for(Index i(0); i < n; ++i)
    {
      frc.axpy(vchr[i], vdef[i]);
    }

    forces.format();
    forces.template copy_n<dim_>(frc);
  }

  /**
   * \brief Navier-Stokes System Level class
   */
  template<
    int dim_,
    typename DataType_ = DFG95::DataType,
    typename IndexType_ = Index,
    typename MatrixBlockA_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, dim_>,
    typename MatrixBlockB_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, dim_, 1>,
    typename MatrixBlockD_ = LAFEM::SparseMatrixBCSR<DataType_, IndexType_, 1, dim_>,
    typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<DataType_, IndexType_>>
  class NavierStokesBlockedSystemLevel :
    public Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_>
  {
  public:
    typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, DataType_, IndexType_, MatrixBlockA_, MatrixBlockB_, MatrixBlockD_, ScalarMatrix_> BaseClass;

    // the noflow velocity filter, used in bench3 simulation
    typename BaseClass::LocalVeloFilter local_velo_filter_noflow;

    // the filtered local system matrix for Vanka
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    // the local velocity mass matrix
    MatrixBlockA_ local_velo_mass_matrix;

    template<typename Trafo_, typename SpaceV_>
    void assemble_velocity_laplace_matrix(Assembly::DomainAssembler<Trafo_>& dom_asm,
      const SpaceV_& space_velo, const String& cubature, const DataType_ nu, bool defo)
    {
      auto& loc_a = this->matrix_a.local();
      loc_a.format();
      Assembly::Common::LaplaceOperatorBlocked<dim_> lapl_op;
      Assembly::Common::DuDvOperatorBlocked<dim_> dudv_op;
      if(defo)
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_a, dudv_op, space_velo, cubature, nu);
      else
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_a, lapl_op, space_velo, cubature, nu);
    }

    template<typename Trafo_, typename SpaceV_>
    void assemble_velocity_mass_matrix(Assembly::DomainAssembler<Trafo_>& dom_asm, const SpaceV_& space_velo, const String& cubature)
    {
      local_velo_mass_matrix = this->matrix_a.local().clone(LAFEM::CloneMode::Weak);
      local_velo_mass_matrix.format();

      Assembly::Common::IdentityOperatorBlocked<dim_> id_op;
      Assembly::assemble_bilinear_operator_matrix_1(dom_asm, local_velo_mass_matrix, id_op, space_velo, cubature);
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
} // namespace DFG95
