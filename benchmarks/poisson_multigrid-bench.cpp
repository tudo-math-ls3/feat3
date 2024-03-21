// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/jacobi_precond.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>
#include <control/statistics.hpp>

namespace PoissonMultigridBench
{
  using namespace FEAT;

  typedef double DataType;
  typedef std::uint64_t IndexType;

  static_assert(sizeof(DataType) == 8u, "invalid data type size");
  static_assert(sizeof(IndexType) == 8u, "invalid index type size");

  static constexpr int dim = 2;

  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange2::Element<TrafoType> SpaceType;

  struct Counts
  {
    static constexpr std::size_t num_ranks = 0u;
    static constexpr std::size_t num_elems = 1u;
    static constexpr std::size_t num_dofs_l = 2u; // local dofs
    static constexpr std::size_t num_dofs_g = 3u; // global dofs [ != sum(num_dofs_l) ]
    static constexpr std::size_t num_nze = 4u;
    static constexpr std::size_t bytes_domain = 5u;
    static constexpr std::size_t bytes_system = 6u;
    static constexpr std::size_t bytes_solver = 7u;
    static constexpr std::size_t elems_mirror = 8u;
    static constexpr std::size_t count = 9u;
  };

  struct Times
  {
    static constexpr std::size_t asm_total = 0u;
    static constexpr std::size_t asm_gate = 1u;
    static constexpr std::size_t asm_muxer = 2u;
    static constexpr std::size_t asm_transfer = 3u;
    static constexpr std::size_t asm_matrix = 4u;
    static constexpr std::size_t gmg_total = 5u;
    static constexpr std::size_t gmg_defect = 6u;
    static constexpr std::size_t gmg_smooth = 7u;
    static constexpr std::size_t gmg_transfer = 8u;
    static constexpr std::size_t gmg_coarse = 9u;
    static constexpr std::size_t count = 10u;
  };

  template<typename T_>
  inline T_ sum(const std::vector<T_>& dv)
  {
    T_ t = T_(0);
    for(auto x : dv) t += x;
    return t;
  }

  template<typename T_, std::size_t n_>
  inline T_ sum(const std::vector<std::array<T_, n_>>& dv, std::size_t i)
  {
    XASSERT(i < n_);
    T_ t = T_(0);
    for(const auto& x : dv) t += x[i];
    return t;
  }

  struct BenchStats
  {
    // per [level][index]
    std::vector<std::array<unsigned long long, Counts::count>> counts, counts_sum, counts_max;

    // per [level][index]
    std::vector<std::array<double, Times::count>> times;

    // (physical, virtual)
    std::array<unsigned long long, 2> mem_use, mem_use_sum, mem_use_max;

    double toe_asm_rhs;

    // final relative defect and final error
    double final_rel_def, final_error;

    explicit BenchStats(std::size_t virt_size) :
      counts(virt_size),
      counts_sum(virt_size),
      counts_max(virt_size),
      times(virt_size),
      toe_asm_rhs(0.0),
      final_rel_def(0.0),
      final_error(0.0)
    {
      for(std::size_t i(0u); i < virt_size; ++i)
      {
        for(std::size_t j(0u); j < Counts::count; ++j)
          counts[i][j] = 0ull;
        for(std::size_t j(0u); j < Times::count; ++j)
          times[i][j] = 0.0;
      }
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(counts.data(), counts_sum.data(), counts.size()*Counts::count, Dist::dt_unsigned_long_long, Dist::op_sum);
      comm.allreduce(counts.data(), counts_max.data(), counts.size()*Counts::count, Dist::dt_unsigned_long_long, Dist::op_max);

      comm.allreduce(mem_use.data(), mem_use_sum.data(), 2u, Dist::op_sum);
      comm.allreduce(mem_use.data(), mem_use_max.data(), 2u, Dist::op_max);
    }

    String format() const
    {
      String s;

      const std::size_t total_elems = counts_sum[0][Counts::num_elems];
      const std::size_t total_gdofs = counts[0][Counts::num_dofs_g];
      const double total_asm_time = sum(times, Times::asm_total);
      const double total_gmg_time = sum(times, Times::gmg_total);
      const double asm_mdofs = (total_asm_time > 1E-5 ? 1E-6 * double(total_gdofs) / total_asm_time : 0.0);
      const double gmg_mdofs = (total_gmg_time > 1E-5 ? 1E-6 * double(total_gdofs) / total_gmg_time : 0.0);

      s += "\nFinal Relative Defect......: " + stringify_fp_sci(final_rel_def);
      s += "\nFinal Normalized Error.....: " + stringify_fp_sci(final_error);

      s += "\n\nMultigrid Timings:\n";
      s += "                Defect /     Smoother /     Transfer /       Coarse /        Total\n";
      s += "Overall : " +
        stringify_fp_fix(sum(times, Times::gmg_defect), 6, 12) + " / " +
        stringify_fp_fix(sum(times, Times::gmg_smooth), 6, 12) + " / " +
        stringify_fp_fix(sum(times, Times::gmg_transfer), 6, 12) + " / " +
        stringify_fp_fix(sum(times, Times::gmg_coarse), 6, 12) + " / " +
        stringify_fp_fix(sum(times, Times::gmg_total), 6, 12) + "\n";
      for(std::size_t i(0); i < times.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(times[i][Times::gmg_defect], 6, 12) + " / " +
          stringify_fp_fix(times[i][Times::gmg_smooth], 6, 12) + " / " +
          stringify_fp_fix(times[i][Times::gmg_transfer], 6, 12) + " / " +
          stringify_fp_fix(times[i][Times::gmg_coarse], 6, 12) + " / " +
          stringify_fp_fix(times[i][Times::gmg_total], 6, 12) + "\n";
      }

      s += "\nAssembly Timings:\n";
      s += "                Gate        Muxer     Transfer       Matrix        Total        RHS\n";
      s += "Overall : " +
        stringify_fp_fix(sum(times, Times::asm_gate), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_muxer), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_transfer), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_matrix), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_total), 6, 10) + " / " +
        stringify_fp_fix(toe_asm_rhs, 6) + "\n";
      for(std::size_t i(0); i < times.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(times[i][Times::asm_gate], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_muxer], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_transfer], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_matrix], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_total], 6, 10) + "\n";
      }

      s += "\nBasic Statistics:\n";
      s += "           Ranks       Elements [    per Patch ]           Dofs [    per Patch ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify(counts_max[i][Counts::num_ranks]).pad_front(6) + " / " +
          stringify(counts_sum[i][Counts::num_elems]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_elems]).pad_front(12) + " ] / " +
          stringify(counts    [i][Counts::num_dofs_g]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_dofs_l]).pad_front(12) + " ]\n";
      }

      s += "\nMirror #elems Statistics:\n";
      s += "Overall : " + stringify(sum(counts_sum, Counts::elems_mirror)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::elems_mirror)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::elems_mirror]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::elems_mirror]).pad_front(15) + " ]\n";
      }

      s += "\nMemory Usage Statistics:\n";
      s += String("Peak Physical") +
        ": " + stringify_fp_fix(double(mem_use_sum[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB " +
        "[ " + stringify_fp_fix(double(mem_use_max[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB ]\n";
      //": " + stringify(mem_use_sum[0]).pad_front(15) +
      //" [" + stringify(mem_use_max[0]).pad_front(15) + " ]\n";
      s += String("Peak Virtual.") +
        ": " + stringify_fp_fix(double(mem_use_sum[1])/(1024.0*1024.0*1024.0), 6, 15) + " GiB " +
        "[ " + stringify_fp_fix(double(mem_use_max[1])/(1024.0*1024.0*1024.0), 6, 15) + " GiB ]\n";
      //": " + stringify(mem_use_sum[1]).pad_front(15) +
      //" [" + stringify(mem_use_max[1]).pad_front(15) + " ]\n";

      s += "\nDomain Level Bytes Statistics:\n";
      s += "Overall : " + stringify(sum(counts_sum, Counts::bytes_domain)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::bytes_domain)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::bytes_domain]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::bytes_domain]).pad_front(15) + " ]\n";
      }

      s += "\nSystem Level Bytes Statistics:\n";
      s += "Overall : " + stringify(sum(counts_sum, Counts::bytes_system)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::bytes_system)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::bytes_system]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::bytes_system]).pad_front(15) + " ]\n";
      }

      unsigned long long domain_bytes = sum(counts_sum, Counts::bytes_domain);
      unsigned long long system_bytes = sum(counts_sum, Counts::bytes_system);
      unsigned long long solver_bytes = sum(counts_sum, Counts::bytes_solver);

      unsigned long long total_bytes = domain_bytes + system_bytes + solver_bytes + 3ull*total_gdofs*sizeof(DataType);

      s += "\n\nOverall Benchmark Summary:";
      s += "\nTotal Number of Elements...: " + stringify(total_elems).pad_front(15);
      s += "\nTotal Number of Global DOFs: " + stringify(total_gdofs).pad_front(15);
      s += "\nTotal Domain Memory Usage..: " + stringify_fp_fix(double(domain_bytes)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(domain_bytes) / double(total_bytes), 2, 6) + "%]";
      s += "\nTotal System Memory Usage..: " + stringify_fp_fix(double(system_bytes)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(system_bytes) / double(total_bytes), 2, 6) + "%]";
      s += "\nTotal Solver Memory Usage..: " + stringify_fp_fix(double(solver_bytes)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(solver_bytes) / double(total_bytes), 2, 6) + "%]";
      s += "\nTotal Combined Memory Usage: " + stringify_fp_fix(double(total_bytes)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += "\nMemory Usage Reported By OS: " + stringify_fp_fix(double(mem_use_sum[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*(double(mem_use_sum[0])/double(total_bytes) - 1.0), 2, 6, true) + "%]";
      s += "\nTotal Assembly Runtime.....: " + stringify_fp_fix(total_asm_time, 6, 15) + " seconds";
      s += "\nTotal Multigrid Runtime....: " + stringify_fp_fix(total_gmg_time, 6, 15) + " seconds";
      s += "\nAssembly Performance.......: " + stringify_fp_fix(asm_mdofs, 6, 15) + " MDOF/s";
      s += "\nSolver Performance.........: " + stringify_fp_fix(gmg_mdofs, 6, 15) + " MDOF/s";

      return s;
    }
  };

  double parse_memory(SimpleArgParser& args)
  {
    String s = args.query("memory")->second.front();
    double f = 1.0;
    if(s.ends_with("k") || s.ends_with("K"))
    {
      f = 1024.0;
      s.pop_back();
    }
    if(s.ends_with("m") || s.ends_with("M"))
    {
      f = 1024.0*1024.0;
      s.pop_back();
    }
    if(s.ends_with("g") || s.ends_with("G"))
    {
      f = 1024.0*1024.0*1024.0;
      s.pop_back();
    }
    if(s.ends_with("t") || s.ends_with("T"))
    {
      f = 1024.0*1024.0*1024.0*1024.0;
      s.pop_back();
    }
    double t = 0.0;
    if(s.parse(t))
      return f*t;
    return 0.0;
  }

  typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;

  void main(int argc, char** argv)
  {
    Dist::Comm comm(Dist::Comm::world());

    Backend::set_preferred_backend(PreferredBackend::generic);

    // enable solver expressions
    Statistics::enable_solver_expressions = true;

    SimpleArgParser args(argc, argv);
    Control::Domain::add_supported_pdc_args(args);

    args.support("memory", "<bytes[K|M|G|T|P]>\n"
      "Specifies the maximum total amount of memory that the benchmark may allocate, summed up over all processes.\n"
      "The size may be suffixed by one of the common powers of 1024, e.g. 50G stands for 50 GibiBytes (aka Gigabytes).\n"
      "If this option is specified, then the number of slices, which is usually set by the --slices option, as well\n"
      "as the refinement levels, which are specified by the --level option, are chosen automatically so that the\n"
      "requested memory usage limit is respected.");
    args.support("slices", "<k>\n"
      "Specifies the number of element slices of the level-0 base mesh in each dimension, so that the base mesh\n"
      "corresponds to a <k>-by-<k> rectilinear mesh; default: 1.");
    args.support("level", "<fine [coarse...]>\n"
      "Specifies the refinement level of the fine mesh on which the benchmark system is solved as well as (optionally)\n"
      "the coarse level for the multigrid hierarchy, potentially including the intermediate levels for the recursive\n"
      "partitioning layers, as required by the Control::PartiDomainControl class. If only the fine level is given,\n"
      "then the multigrid hierarchy is automatically extended down to level 0 on a single MPI process.");
    args.support("cycle", "<V|F|W>\nSpecifies which multigrid cycle to use; default: V");
    args.support("iters", "<n>\nSpecifies the number of PCG-GMG to perform; default: 5");
    args.support("steps", "<n>\nSpecifies the number of Jacobi smoothing steps; default: 5");
    args.support("no-shrink", "\nDon't shrink grid transfer matrices (not recommended).");
    args.support("ext-stats", "\nPrint extended solver and MPI statistics.");
    args.support("test","\nRuns the benchmark in regression test mode.");
    args.support("backend", "<generic|cuda|mkl>\n"
      "Specifies which backend to use for the actual PCG-GMG solution process; default: generic");

    // no arguments given?
    if(args.num_args() <= 1)
    {
      comm.print(
        "FEAT3 2D-Q2 Poisson-Multigrid Solver Benchmark\n"
        "----------------------------------------------\n"
        "This is a small benchmark application that will create a rectilinear 2D grid, then partition and\n"
        "distribute it to an arbitrary number of MPI processes to solve a Poisson PDE discretized by a\n"
        "conforming Q2 finite element scheme and (approximately) solve it by a geometric multigrid\n"
        "preconditioned conjugate gradient solver (quite similar to the famous HPCG benchmark).\n"
        "This tool will then finally report a large number of statistics along with MDOF/s ratings for\n"
        "both the assembly as well as the solution process to provide a means to evaluate the real-world\n"
        "performance of the underlying hardware and software stack.\n\n"
        "At the very least, you have to specify either a maximum refinement level by supplying the\n"
        "--level <n> option or a maximum memory usage limit by supplying the --memory <bound> option.\n\n"
        "This benchmark supports the following command line options:"
      );

      comm.print(args.get_supported_help());

      return;
    }

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());

      // abort
      FEAT::Runtime::abort();
    }

    comm.print(String(120u, '#'));
    comm.print("Dist Comm World Size: " + stringify(comm.size()));

    Index multigrid_iters = 5;
    Index smooth_steps = 5;
    DataType smooth_damp = 0.7;
    String backend("generic");

    Solver::MultiGridCycle multigrid_cycle(Solver::MultiGridCycle::V);

    args.parse("cycle", multigrid_cycle);
    args.parse("iters", multigrid_iters);
    args.parse("steps", smooth_steps);
    bool no_shrink = (args.check("no-shrink") >= 0);

    // choose backend
    args.parse("backend", backend);
    if(backend.compare_no_case("cuda") == 0)
    {
#ifndef FEAT_HAVE_CUDA
      comm.print(std::cerr, "ERROR: 'cuda' backend is only available when FEAT is configured with CUDA support enabled");
      Runtime::abort();
#endif // not FEAT_HAVE_CUDA
    }
    else if(backend.compare_no_case("mkl") == 0)
    {
#ifndef FEAT_HAVE_MKL
      comm.print(std::cerr, "ERROR: 'mkl' backend is only available when FEAT is configured with MKL support enabled");
      Runtime::abort();
#endif // FEAT_HAVE_MKL
    }
    else if(backend.compare_no_case("generic") != 0)
    {
      comm.print(std::cerr, "ERROR: unknown backend '" + backend + "'");
      Runtime::abort();
    }

    // create domain control
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);
    if(!domain.parse_args(args))
      Runtime::abort();

    Index base_mesh_num_slices = 1u;

    // How do we choose the problem size? There are 3 possibilities:
    // 1. run in regression test mode by specifying --test
    // 2. specify a memory usage limit by --memory
    // 3. specify a refinement level by --level
    if(args.check("test") >= 0)
    {
      // run on level 4
      int fine_mesh_level = 4;

      // choose desired levels and base mesh dimensions
      if(comm.size() > 1)
        domain.set_desired_levels(fine_mesh_level, 0, 0);
      else
        domain.set_desired_levels(fine_mesh_level, 0);
    }
    else if(args.check("memory") > 0)
    {
      // parse maximum memory requirement
      double max_memory = parse_memory(args);
      if(max_memory < 1.0)
      {
        comm.print(std::cerr, "ERROR: failed to parse memory requirement");
        FEAT::Runtime::abort();
      }

      // estimate of bytes per element
      const double bytes_per_elem = 3000.0;

      // compute log_4 of memory requirement = level w.r.t. 1x1 base mesh
      double log4_memory = Math::log(max_memory / bytes_per_elem) * 0.72134752;

      // compute fine level w.r.t 2x2 (or 3x3) base mesh
      int fine_mesh_level = int(log4_memory) - 1;

      // compute memory requirement based on 2x2=4 and 3x3=9 slices, resp
      double fine_elems_base = double(1ull << (2*fine_mesh_level));
      double slc2_memory = 4.0 * bytes_per_elem * fine_elems_base;
      double slc3_memory = 9.0 * bytes_per_elem * fine_elems_base;

      // use 3x3 slices for base mesh if possible, otherwise use 2x2 slices
      if(slc3_memory <= max_memory)
        base_mesh_num_slices = 3u;
      else if(slc2_memory <= max_memory)
        base_mesh_num_slices = 2u;
      else // this should never happen
      {
        comm.print(std::cerr, "ERROR: internal error in maximum memory resolution");
        FEAT::Runtime::abort();
      }

      // choose desired levels and base mesh dimensions
      if(comm.size() > 1)
        domain.set_desired_levels(fine_mesh_level, 0, 0);
      else
        domain.set_desired_levels(fine_mesh_level, 0);
    }
    else if(args.check("level") > 1)
    {
      // multiple levels given; pass them directly to the domain control
      domain.set_desired_levels(args.query("level")->second);
    }
    else if(args.check("level") == 1)
    {
      // just fine mesh level given; go down to level 0
      int level = 0;
      args.parse("level", level);
      if(comm.size() > 1)
        domain.set_desired_levels(level, 0, 0);
      else
        domain.set_desired_levels(level, 0);
    }
    else
    {
      comm.print(std::cerr, "ERROR: you must specify either a memory requirement by --memory <size> or a fine mesh level by --level <level>");
      FEAT::Runtime::abort();
    }

    // parse slices; this may override the default setting defined by the --memory option handling above
    args.parse("slices", base_mesh_num_slices);

    // create domain from a rectilinear base mesh
    domain.create_rectilinear(base_mesh_num_slices, base_mesh_num_slices);

    // get the number of slices on the fine mesh
    Index fine_mesh_num_slices = base_mesh_num_slices << (domain.get_chosen_levels().front().first);

    // print partitioning info
    //          12345678901234567890
    comm.print("Base Mesh Dimensions: " + stringify(base_mesh_num_slices) + "^2 = " +
      stringify(base_mesh_num_slices*base_mesh_num_slices) + " elements");
    comm.print("Fine Mesh Dimensions: " + stringify(fine_mesh_num_slices) + "^2 = " +
      stringify(fine_mesh_num_slices*fine_mesh_num_slices) + " elements");
    comm.print("Desired Levels......: " + domain.format_desired_levels());
    comm.print("Chosen Levels.......: " + domain.format_chosen_levels());
    comm.print("Partitioner Info....: " + domain.get_chosen_parti_info());
    comm.print("Preferred Backend...: " + backend);
    comm.print("Test Mode...........: " + String((args.check("test") >= 0) ? "yes" : "no"));
    comm.print("Memory Limit........: " + (args.check("memory") > 0 ? args.query("memory")->second.front() : String("-N/A-")));
    /*if(args.check("memory") > 0)
      comm.print("Memory Limit........: " + args.query("memory")->second.front());
    else
      comm.print("Memory Limit........: -N/A-");*/
    comm.print("Multigrid Iterations: " + stringify(multigrid_iters));
    comm.print("Multigrid Cycle.....: " + stringify(multigrid_cycle));
    comm.print("Smoothing Steps.....: " + stringify(smooth_steps));
    comm.print("Smoother Damping....: " + stringify(smooth_damp));
    comm.print("Shink Matrices......: " + String(no_shrink ? "no" : "yes"));

    // create benchmark statistics
    BenchStats stats(domain.size_virtual());

    {
      std::size_t n = domain.size_virtual() > domain.size_physical() ? domain.size_physical()+1 : domain.size_physical();
      for(std::size_t i(0); i < n; ++i)
      {
        const auto& vdl = domain.at(i);
        if(vdl.is_parent())
          stats.counts[i][Counts::bytes_domain] = vdl.level_c().bytes() + vdl.level_p().bytes();
        else if(vdl.is_child())
          stats.counts[i][Counts::bytes_domain] = vdl.level_c().bytes();
        else
          stats.counts[i][Counts::bytes_domain] = vdl.level().bytes();
      }
    }

    const Index num_levels = domain.size_physical();

    typedef Control::ScalarUnitFilterSystemLevel<DataType, IndexType> SystemLevelType;

    // create system levels
    std::deque<std::shared_ptr<SystemLevelType>> system;

    for (Index i(0); i < num_levels; ++i)
    {
      system.push_back(std::make_shared<SystemLevelType>());
    }

    const String cubature("gauss-legendre:" + stringify(SpaceType::local_degree+1));

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      domain.at(i)->domain_asm.compile_all_elements();
      system.at(i)->assemble_gate(domain.at(i));
      stats.times[i][Times::asm_gate] += ts.elapsed_now();
      if((i+1) < domain.size_virtual())
      {
        TimeStamp ts2;
        system.at(i)->assemble_coarse_muxer(domain.at(i+1));
        stats.times[i][Times::asm_muxer] += ts2.elapsed_now();
        TimeStamp ts3;
        system.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);

        if(!no_shrink)
        {
          system.at(i)->transfer_sys.get_mat_prol().shrink(1E-3);
          system.at(i)->transfer_sys.get_mat_rest().shrink(1E-3);
        }

        stats.times[i][Times::asm_transfer] += ts3.elapsed_now();
      }
      stats.times[i][Times::asm_total] += ts.elapsed_now();
    }

    /* ***************************************************************************************** */

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system.at(i)->assemble_laplace_matrix(domain.at(i)->domain_asm, domain.at(i)->space, cubature);
      double tt = ts.elapsed_now();
      stats.times[i][Times::asm_total] += tt;
      stats.times[i][Times::asm_matrix] += tt;
    }

    /* ***************************************************************************************** */

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system.at(i)->assemble_homogeneous_unit_filter(*domain.at(i), domain.at(i)->space);

      // apply filter
      system.at(i)->filter_sys.local().filter_mat(system.at(i)->matrix_sys.local());
      stats.times[i][Times::asm_total] += ts.elapsed_now();

      stats.counts[i][Counts::num_ranks] = Index(domain.at(i).layer().comm().size());
      stats.counts[i][Counts::num_elems] = domain.at(i)->get_mesh().get_num_elements();
      stats.counts[i][Counts::num_dofs_g] = system.at(i)->matrix_sys.rows();
      stats.counts[i][Counts::num_dofs_l] = system.at(i)->matrix_sys.local().rows();
      stats.counts[i][Counts::num_nze] = system.at(i)->matrix_sys.local().used_elements();
      stats.counts[i][Counts::bytes_system] = system.at(i)->bytes();
      stats.counts[i][Counts::elems_mirror] = 0;
      auto& sys_mirrors =  system.at(i)->gate_sys._mirrors;
      for (auto& mirror : sys_mirrors)
      {
        stats.counts[i][Counts::elems_mirror] += mirror.num_indices();
      }
    }

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system.front();

    TimeStamp ts_rhs;

    // create new vector
    GlobalSystemVector vec_ref = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

    // interpolate sine-bubble to obtain reference solution
    Analytic::Common::SineBubbleFunction<dim> sine_bubble_function;
    Assembly::Interpolator::project(vec_ref.local(), sine_bubble_function, the_domain_level.space);

    // compute RHS from reference solution
    the_system_level.matrix_sys.apply(vec_rhs, vec_ref);

    // format solution vector
    vec_sol.format();

    // and filter it
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_rhs);

    stats.toe_asm_rhs = ts_rhs.elapsed_now();

    Statistics::toe_assembly = stamp_ass.elapsed_now();

    stats.counts.front()[Counts::bytes_system] += 3ull * vec_sol.local().size() * sizeof(DataType);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // from here on, we can switch to another backend if the user desires this
    if(backend.compare_no_case("cuda") == 0)
      Backend::set_preferred_backend(FEAT::PreferredBackend::cuda);
    else if(backend.compare_no_case("mkl") == 0)
      Backend::set_preferred_backend(FEAT::PreferredBackend::mkl);

    auto multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalSystemMatrix,
      typename SystemLevelType::GlobalSystemFilter,
      typename SystemLevelType::GlobalSystemTransfer
      > >(domain.size_virtual());

    // push all levels except the coarse most one
    for (Index i(0); i < num_levels; ++i)
    {
      const SystemLevelType& lvl = *system.at(i);
      const std::size_t lvl_dofs = lvl.matrix_sys.local().rows();

      if((i+1) < domain.size_virtual())
      {
        auto jacobi = Solver::new_jacobi_precond(lvl.matrix_sys, lvl.filter_sys, smooth_damp);
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, 1.0, jacobi);
        smoother->set_min_iter(smooth_steps);
        smoother->set_max_iter(smooth_steps);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);

        // Jacobi: 1 Vector
        // Richardson: 2 Vectors
        // Multigrid: 5 Vectors
        stats.counts.at(i)[Counts::bytes_solver] = 8ull * lvl_dofs * sizeof(DataType);
      }
      else
      {
        auto cgsolver = Solver::new_jacobi_precond(lvl.matrix_sys, lvl.filter_sys, 1.0);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);

        // Jacobi: 1 Vector
        stats.counts.at(i)[Counts::bytes_solver] = lvl_dofs * sizeof(DataType);
      }
    }

    // create PCG-GMG solver
    auto mgv = Solver::new_multigrid(multigrid_hierarchy, multigrid_cycle);
    auto solver = Solver::new_pcg(the_system_level.matrix_sys, the_system_level.filter_sys, mgv);

    // PCG: 3 Vectors
    stats.counts.front()[Counts::bytes_solver] += 3ull * the_system_level.matrix_sys.local().rows() * sizeof(DataType);

    // set tolerance
    solver->set_plot_name("PCG-Multigrid-Jacobi");
    solver->skip_defect_calc(false);
    //solver->set_tol_rel(1E-12);
    //solver->set_max_iter(1000);

    // initialize
    multigrid_hierarchy->init();
    solver->init();

    Statistics::reset();

    comm.print("\nWarming up...");

    // perform a warmup run with 3 iterations only to ensure that all vectors are allocated
    solver->set_min_iter(3);
    solver->set_max_iter(3);
    solver->set_plot_mode(Solver::PlotMode::none);

    // warmup solve
    solver->apply(vec_sol, vec_rhs);

    // reformat solution vector
    vec_sol.format();

    comm.print("Benchmarking...\n");

    // set desired number of iterations and enable plotting
    solver->set_min_iter(multigrid_iters);
    solver->set_max_iter(multigrid_iters);
    solver->set_plot_mode(Solver::PlotMode::summary);

    Statistics::reset();

    TimeStamp at;

    // benchmarking solve
    Solver::Status result = solver->apply(vec_sol, vec_rhs);

    const double solver_toe(at.elapsed_now());

    if (!Solver::status_success(result))
    {
      comm.print("Solver execution FAILED, with status: " + stringify(result));
    }

    // save final relative residual
    stats.final_rel_def = solver->get_def_final() / solver->get_def_initial();

    // release solver
    solver->done();
    multigrid_hierarchy->done();

    // reset backend to generic
    Backend::set_preferred_backend(PreferredBackend::generic);

    // compute and save final relative error
    {
      vec_rhs.axpy(vec_sol, vec_ref, -DataType(1));
      stats.final_error = vec_rhs.norm2() / Math::sqrt(DataType(stats.counts.front()[Counts::num_dofs_g]));
    }

    // get memory info
    {
      MemoryUsage meminfo;
      stats.mem_use[0] = meminfo.get_peak_physical();
      stats.mem_use[1] = meminfo.get_peak_virtual();
    }

    stats.sync(comm);

    // set multigrid timings
    for(Index i(0); i < multigrid_hierarchy->size_physical(); ++i)
    {
      stats.times.at(i)[Times::gmg_defect] = multigrid_hierarchy->get_time_defect(int(i));
      stats.times.at(i)[Times::gmg_smooth] = multigrid_hierarchy->get_time_smooth(int(i));
      stats.times.at(i)[Times::gmg_transfer] = multigrid_hierarchy->get_time_transfer(int(i));
      stats.times.at(i)[Times::gmg_coarse] = multigrid_hierarchy->get_time_coarse(int(i));
      stats.times.at(i)[Times::gmg_total] =
        stats.times.at(i)[Times::gmg_defect] + stats.times.at(i)[Times::gmg_smooth] +
        stats.times.at(i)[Times::gmg_transfer] + stats.times.at(i)[Times::gmg_coarse];
    }

    comm.print(stats.format());

    if(args.check("ext-stats") >= 0)
    {
      FEAT::Control::Statistics::report(solver_toe, 0, MeshType::ShapeType::dimension, system, domain);
      comm.print(FEAT::Statistics::get_formatted_solver_internals());
      comm.print(FEAT::Statistics::get_formatted_solver_tree().trim());
    }

    if(args.check("test") >= 0)
    {
      // check final relative defect and normalized error
      bool bdef = (stats.final_rel_def < 1E-8);
      bool berr = (stats.final_error < 1E-10);
      comm.print(String("\nTEST-MODE: CHECK FINAL DEFECT: ") + (bdef ? "OK" : "FAILED"));
      comm.print(String("TEST-MODE: CHECK FINAL ERROR:  ") + (berr ? "OK" : "FAILED"));
      comm.print((bdef && berr) ? "\nTEST PASSED" : "\nTEST FAILED");
    }
  }
} // namespace PoissonMultigridBench

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  PoissonMultigridBench::main(argc, argv);
  return 0;
}
