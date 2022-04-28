// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once
#ifndef APPLICATIONS_DFG95_BENCH_CCND_COMMON_HPP
#define APPLICATIONS_DFG95_BENCH_CCND_COMMON_HPP 1

#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
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

  // our one and only memory type
  typedef Mem::Main MemType;
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
    DT_ drag_coeff_vol, lift_coeff_vol, drag_err_vol, lift_err_vol;

    // pressure values
    DT_ pres_diff, pres_err;

    // flow through upper/lower region
    DT_ flux_upper, flux_lower;

    // velocity field information
    Assembly::VelocityInfo<DT_, dim_> velo_info;

    static constexpr std::size_t padlen = std::size_t(30);

    BenchmarkSummary() :
      drag_coeff_line(0.0), lift_coeff_line(0.0), drag_err_line(0.0), lift_err_line(0.0),
      drag_coeff_vol(0.0), lift_coeff_vol(0.0), drag_err_vol(0.0), lift_err_vol(0.0),
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
      s += String("Drag Coefficient (Vol)").pad_back(padlen, pc) + ": " + stringify_fp_fix(drag_coeff_vol, prec)
        + "   [ Error: " + stringify_fp_sci(drag_err_vol, prec) + " ]\n";
      s += String("Lift Coefficient (Line)").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff_line, prec)
        + "   [ Error: " + stringify_fp_sci(lift_err_line, prec) + " ]\n";
      s += String("Lift Coefficient (Vol)").pad_back(padlen, pc) + ": " + stringify_fp_fix(lift_coeff_vol, prec)
        + "   [ Error: " + stringify_fp_sci(lift_err_vol, prec) + " ]\n";
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
      s += prefix + "DC/LC/PD: ";
      s += stringify_fp_fix(drag_coeff_line, prec) + " " + stringify_fp_fix(drag_coeff_vol, prec) + " ";
      s += stringify_fp_fix(lift_coeff_line, prec) + " " + stringify_fp_fix(lift_coeff_vol, prec) + " ";
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

  // accumulator for benchmark body forces (i.e. drag and lift)
  // this class is used for the computation of the 'line integration' variants of drag and lift
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
      // question to self: is it even necessary to adjust something?

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
      // question to self: is it even necessary to adjust something?

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

  // accumulator for computation of X-flux
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

  // computes the bbody forces by the 'volume integration' approach
  template<
    typename DataType_,
    int dim_,
    typename VectorTypeV_,
    typename VectorTypeP_,
    typename VectorTypeC_,
    typename SpaceV_,
    typename SpaceP_,
    typename Cubature_>
  void assemble_bdforces_vol(
    //Tiny::Tensor3<DataType_, 3, 2, dim_>& forces,
    Tiny::Matrix<DataType_, 2, dim_>& forces,
    const VectorTypeV_& vec_sol_v,
    const VectorTypeP_& vec_sol_p,
    const VectorTypeC_& vec_char,
    const SpaceV_& space_v,
    const SpaceP_& space_p,
    const Cubature_& cubature_factory)
  {
    forces.format();

    // assembly traits
    typedef Assembly::AsmTraits2<
      DataType_,
      SpaceV_,
      SpaceP_,
      TrafoTags::jac_det,
      SpaceTags::grad,
      SpaceTags::value> AsmTraits;

    // fetch the trafo
    const typename AsmTraits::TrafoType& trafo = space_v.get_trafo();

    // create a trafo evaluator
    typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

    // create space evaluators
    typename AsmTraits::TestEvaluator velo_eval(space_v);
    typename AsmTraits::TrialEvaluator pres_eval(space_p);

    // create dof-mappings
    typename AsmTraits::TestDofMapping velo_dof_mapping(space_v);
    typename AsmTraits::TrialDofMapping pres_dof_mapping(space_p);

    // create trafo evaluation data
    typename AsmTraits::TrafoEvalData trafo_data;

    // create space evaluation data
    typename AsmTraits::TestEvalData velo_data;
    typename AsmTraits::TrialEvalData pres_data;

    // maximum number of dofs
    static constexpr int max_velo_dofs = AsmTraits::max_local_test_dofs;
    static constexpr int max_pres_dofs = AsmTraits::max_local_trial_dofs;

    // create cubature rule
    typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

    // create vector gathers
    typename VectorTypeV_::GatherAxpy gather_velo(vec_sol_v);
    typename VectorTypeP_::GatherAxpy gather_pres(vec_sol_p);
    typename VectorTypeC_::GatherAxpy gather_char(vec_char);

    // create local vector data
    typedef typename VectorTypeV_::ValueType VeloValue;
    typedef typename VectorTypeP_::ValueType PresValue;
    typedef typename VectorTypeC_::ValueType CharValue;
    typedef Tiny::Vector<VeloValue, max_velo_dofs> LocalVectorTypeV;
    typedef Tiny::Vector<PresValue, max_pres_dofs> LocalVectorTypeP;
    typedef Tiny::Vector<CharValue, max_velo_dofs> LocalVectorTypeC;

    LocalVectorTypeV local_velo_dofs;
    LocalVectorTypeP local_pres_dofs;
    LocalVectorTypeC local_char_dofs;//, local_char_dofs2, local_char_dofs3;

    // our local velocity gradient
    Tiny::Matrix<DataType_, dim_, dim_> loc_grad_v;
    PresValue loc_pres;
    Tiny::Vector<DataType_, dim_> loc_nu_v, loc_char1;//, loc_char2, loc_char3;

    loc_grad_v.format();
    loc_pres = DataType_(0);
    loc_char1.format();
    //loc_char2.format();
    //loc_char3.format();

    // loop over all cells of the mesh
    for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
    {
      // prepare trafo evaluator
      trafo_eval.prepare(cell);

      // prepare space evaluator
      velo_eval.prepare(trafo_eval);
      pres_eval.prepare(trafo_eval);

      // initialize dof-mapping
      velo_dof_mapping.prepare(cell);
      pres_dof_mapping.prepare(cell);

      // fetch number of local dofs
      const int num_loc_velo_dofs = velo_eval.get_num_local_dofs();
      const int num_loc_pres_dofs = pres_eval.get_num_local_dofs();

      // gather our local dofs
      local_velo_dofs.format();
      local_pres_dofs.format();
      local_char_dofs.format();
      gather_velo(local_velo_dofs, velo_dof_mapping);
      gather_pres(local_pres_dofs, pres_dof_mapping);
      gather_char(local_char_dofs, velo_dof_mapping);

      // note: the following were experimental variants of the characteristic function that were
      // chose in hope of obtaining higher EOCs. some seemed to work in 2D, but all failed in 3D.
      /*if(dim == 2)
      {
        // build char dofs 2
        local_char_dofs2[0] = local_char_dofs[0];
        local_char_dofs2[1] = local_char_dofs[1];
        local_char_dofs2[2] = local_char_dofs[2];
        local_char_dofs2[3] = local_char_dofs[3];
        local_char_dofs2[4] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[1]));
        local_char_dofs2[5] = Math::sqr(DataType(0.5) * (local_char_dofs[2] + local_char_dofs[3]));
        local_char_dofs2[6] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[2]));
        local_char_dofs2[7] = Math::sqr(DataType(0.5) * (local_char_dofs[1] + local_char_dofs[3]));
        local_char_dofs2[8] = Math::sqr(DataType(0.25) * (local_char_dofs[0] + local_char_dofs[1] + local_char_dofs[2] + local_char_dofs[3]));

        // build char dofs 3
        local_char_dofs3[0] = local_char_dofs[0];
        local_char_dofs3[1] = local_char_dofs[1];
        local_char_dofs3[2] = local_char_dofs[2];
        local_char_dofs3[3] = local_char_dofs[3];
        local_char_dofs3[4] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[1]));
        local_char_dofs3[5] = Math::sqr(DataType(0.5) * (local_char_dofs[2] + local_char_dofs[3]));
        local_char_dofs3[6] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[2]));
        local_char_dofs3[7] = Math::sqr(DataType(0.5) * (local_char_dofs[1] + local_char_dofs[3]));
        local_char_dofs3[8] = 0.0;
      }
      else
      {
        // build char dofs 2
        local_char_dofs2[ 0] = local_char_dofs[0];
        local_char_dofs2[ 1] = local_char_dofs[1];
        local_char_dofs2[ 2] = local_char_dofs[2];
        local_char_dofs2[ 3] = local_char_dofs[3];
        local_char_dofs2[ 4] = local_char_dofs[4];
        local_char_dofs2[ 5] = local_char_dofs[5];
        local_char_dofs2[ 6] = local_char_dofs[6];
        local_char_dofs2[ 7] = local_char_dofs[7];
        local_char_dofs2[ 8] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[1]));
        local_char_dofs2[ 9] = Math::sqr(DataType(0.5) * (local_char_dofs[2] + local_char_dofs[3]));
        local_char_dofs2[10] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[2]));
        local_char_dofs2[11] = Math::sqr(DataType(0.5) * (local_char_dofs[1] + local_char_dofs[3]));
        local_char_dofs2[12] = Math::sqr(DataType(0.5) * (local_char_dofs[4] + local_char_dofs[5]));
        local_char_dofs2[13] = Math::sqr(DataType(0.5) * (local_char_dofs[6] + local_char_dofs[7]));
        local_char_dofs2[14] = Math::sqr(DataType(0.5) * (local_char_dofs[4] + local_char_dofs[6]));
        local_char_dofs2[15] = Math::sqr(DataType(0.5) * (local_char_dofs[5] + local_char_dofs[7]));
        local_char_dofs2[16] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[4]));
        local_char_dofs2[17] = Math::sqr(DataType(0.5) * (local_char_dofs[1] + local_char_dofs[5]));
        local_char_dofs2[18] = Math::sqr(DataType(0.5) * (local_char_dofs[2] + local_char_dofs[6]));
        local_char_dofs2[19] = Math::sqr(DataType(0.5) * (local_char_dofs[3] + local_char_dofs[7]));
        local_char_dofs2[20] = Math::sqr(DataType(0.25) * (local_char_dofs[0] + local_char_dofs[1] + local_char_dofs[2] + local_char_dofs[3]));
        local_char_dofs2[21] = Math::sqr(DataType(0.25) * (local_char_dofs[4] + local_char_dofs[5] + local_char_dofs[6] + local_char_dofs[7]));
        local_char_dofs2[22] = Math::sqr(DataType(0.25) * (local_char_dofs[0] + local_char_dofs[1] + local_char_dofs[4] + local_char_dofs[5]));
        local_char_dofs2[23] = Math::sqr(DataType(0.25) * (local_char_dofs[2] + local_char_dofs[3] + local_char_dofs[6] + local_char_dofs[7]));
        local_char_dofs2[24] = Math::sqr(DataType(0.25) * (local_char_dofs[0] + local_char_dofs[2] + local_char_dofs[4] + local_char_dofs[6]));
        local_char_dofs2[25] = Math::sqr(DataType(0.25) * (local_char_dofs[1] + local_char_dofs[3] + local_char_dofs[5] + local_char_dofs[7]));
        local_char_dofs2[26] = Math::sqr(DataType(0.125) * (local_char_dofs[0] + local_char_dofs[1] + local_char_dofs[2] + local_char_dofs[3] + local_char_dofs[4] + local_char_dofs[5] + local_char_dofs[6] + local_char_dofs[7]));

        // build char dofs 3
        local_char_dofs3[ 0] = local_char_dofs[0];
        local_char_dofs3[ 1] = local_char_dofs[1];
        local_char_dofs3[ 2] = local_char_dofs[2];
        local_char_dofs3[ 3] = local_char_dofs[3];
        local_char_dofs3[ 4] = local_char_dofs[4];
        local_char_dofs3[ 5] = local_char_dofs[5];
        local_char_dofs3[ 6] = local_char_dofs[6];
        local_char_dofs3[ 7] = local_char_dofs[7];
        local_char_dofs3[ 8] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[1]));
        local_char_dofs3[ 9] = Math::sqr(DataType(0.5) * (local_char_dofs[2] + local_char_dofs[3]));
        local_char_dofs3[10] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[2]));
        local_char_dofs3[11] = Math::sqr(DataType(0.5) * (local_char_dofs[1] + local_char_dofs[3]));
        local_char_dofs3[12] = Math::sqr(DataType(0.5) * (local_char_dofs[4] + local_char_dofs[5]));
        local_char_dofs3[13] = Math::sqr(DataType(0.5) * (local_char_dofs[6] + local_char_dofs[7]));
        local_char_dofs3[14] = Math::sqr(DataType(0.5) * (local_char_dofs[4] + local_char_dofs[6]));
        local_char_dofs3[15] = Math::sqr(DataType(0.5) * (local_char_dofs[5] + local_char_dofs[7]));
        local_char_dofs3[16] = Math::sqr(DataType(0.5) * (local_char_dofs[0] + local_char_dofs[4]));
        local_char_dofs3[17] = Math::sqr(DataType(0.5) * (local_char_dofs[1] + local_char_dofs[5]));
        local_char_dofs3[18] = Math::sqr(DataType(0.5) * (local_char_dofs[2] + local_char_dofs[6]));
        local_char_dofs3[19] = Math::sqr(DataType(0.5) * (local_char_dofs[3] + local_char_dofs[7]));
        local_char_dofs3[20] = Math::sqr(DataType(0.25) * (local_char_dofs[0] + local_char_dofs[1] + local_char_dofs[2] + local_char_dofs[3]));
        local_char_dofs3[21] = Math::sqr(DataType(0.25) * (local_char_dofs[4] + local_char_dofs[5] + local_char_dofs[6] + local_char_dofs[7]));
        local_char_dofs3[22] = Math::sqr(DataType(0.25) * (local_char_dofs[0] + local_char_dofs[1] + local_char_dofs[4] + local_char_dofs[5]));
        local_char_dofs3[23] = Math::sqr(DataType(0.25) * (local_char_dofs[2] + local_char_dofs[3] + local_char_dofs[6] + local_char_dofs[7]));
        local_char_dofs3[24] = Math::sqr(DataType(0.25) * (local_char_dofs[0] + local_char_dofs[2] + local_char_dofs[4] + local_char_dofs[6]));
        local_char_dofs3[25] = Math::sqr(DataType(0.25) * (local_char_dofs[1] + local_char_dofs[3] + local_char_dofs[5] + local_char_dofs[7]));
        local_char_dofs3[26] = 0.0;
      }*/

      // loop over all quadrature points and integrate
      for(int point(0); point < cubature_rule.get_num_points(); ++point)
      {
        // compute trafo data
        trafo_eval(trafo_data, cubature_rule.get_point(point));

        // compute basis function data
        velo_eval(velo_data, trafo_data);
        pres_eval(pres_data, trafo_data);

        // pre-compute cubature weight
        const DataType_ weight = trafo_data.jac_det * cubature_rule.get_weight(point);

        // compute local velocity gradient and characteristic gradient
        loc_grad_v.format();
        loc_char1.format();
        //loc_char2.format();
        //loc_char3.format();
        for(int i(0); i < num_loc_velo_dofs; ++i)
        {
          // update velocity gradient
          loc_grad_v.add_outer_product(local_velo_dofs[i], velo_data.phi[i].grad);
          loc_char1.axpy(local_char_dofs[i], velo_data.phi[i].grad);
          //loc_char2.axpy(local_char_dofs2[i], velo_data.phi[i].grad);
          //loc_char3.axpy(local_char_dofs3[i], velo_data.phi[i].grad);
        }

        // compute local pressure
        loc_pres = DataType_(0);
        for(int i(0); i < num_loc_pres_dofs; ++i)
        {
          loc_pres += local_pres_dofs[i] * pres_data.phi[i].value;
        }

        // V-forces
        forces[0].add_mat_vec_mult(loc_grad_v, loc_char1, -weight);
        //forces[0][0].add_mat_vec_mult(loc_grad_v, loc_char1, -weight);
        //forces[1][0].add_mat_vec_mult(loc_grad_v, loc_char2, -weight);
        //forces[2][0].add_mat_vec_mult(loc_grad_v, loc_char3, -weight);

        // P-forces
        forces[1].axpy(weight * loc_pres, loc_char1);
        //forces[0][1].axpy(weight * loc_pres, loc_char1);
        //forces[1][1].axpy(weight * loc_pres, loc_char2);
        //forces[2][1].axpy(weight * loc_pres, loc_char3);

        // continue with next cubature point
      }

      // finish dof mapping
      pres_dof_mapping.finish();
      velo_dof_mapping.finish();

      // finish evaluators
      pres_eval.finish();
      velo_eval.finish();
      trafo_eval.finish();
    }
  }

  /**
   * \brief Navier-Stokes System Level class
   */
  template<
    int dim_,
    typename MemType_ = Mem::Main,
    typename DataType_ = DFG95::DataType,
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

    // the noflow velocity filter, used in bench3 simulation
    typename BaseClass::LocalVeloFilter local_velo_filter_noflow;

    // the filtered local system matrix for Vanka
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    // the local velocity mass matrix
    MatrixBlockA_ local_velo_mass_matrix;

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
  }; // class NavierStokesBlockedSystemLevel
} // namespace DFG95

#endif // APPLICATIONS_DFG95_BENCH_CCND_COMMON_HPP
