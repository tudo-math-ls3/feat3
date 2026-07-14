#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/assembly/bilinear_operator.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <control/stokes_blocked.hpp>
#include <control/domain/domain_level.hpp>
#include <control/domain/parti_domain_control.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/amavanka_richardson.hpp>

// NVTX profiling ranges via the shared RAII wrapper (no-op without FEAT_HAVE_CUDA)
#include <kernel/util/cuda_util.hpp>
#include <kernel/util/arg_parser.hpp>

#include <deque>
#include <numeric>

// Benchmark comparing the standard Richardson(Schwarz(AmaVanka)) smoother stack against the fused
// CUDA AmaVankaRichardson smoother (and the new class's own generic path), analogous to the
// gpu_gmres-bench comparison of GMRES(Schwarz(AmaVanka)) against the fused AmaVankaGMRES.
namespace RichardsonBench
{
  using namespace FEAT;

  constexpr std::size_t pad_len = 24;
  constexpr char pad_char = '.';

  struct ParamHolder
  {
    Index num_iter = 20;
    double omega = 0.5;
    Solver::PlotMode plot_mode = Solver::PlotMode::none;
    /// standard Richardson(Schwarz(AmaVanka)) on the selected backend
    bool std_richardson = true;
    /// AmaVankaRichardson, generic (CPU reference) path
    bool generic_richardson = true;
    /// AmaVankaRichardson, fused CUDA driver
    bool cuda_richardson = true;
    /// AmaVankaRichardson, fused CUDA driver forced onto the parallel (synchronized) kernel
    bool cuda_parallel_richardson = false;
    bool hotrun = false;
    /// print the final defect norm after each solver run
    bool print_defect = false;
  };

  struct Times
  {
    static constexpr std::size_t assembly = 0u;
    static constexpr std::size_t solve_std = 1u;
    static constexpr std::size_t solve_generic = 2u;
    static constexpr std::size_t solve_cuda = 3u;
    static constexpr std::size_t solve_cuda_parallel = 4u;
    static constexpr std::size_t count = 6u;
  };

  struct Sizes
  {
    static constexpr std::size_t number_elements = 0u;
    static constexpr std::size_t nnz = 1u;
    static constexpr std::size_t velo_dofs = 2u;
    static constexpr std::size_t pres_dofs = 3u;
    static constexpr std::size_t dofs = 4u;
    static constexpr std::size_t count = 6u;
  };

  template<typename Solver_, typename VecType_>
  void run_solver(Solver_& solver, VecType_& vec_cor, const VecType_& vec_def, double& time_total)
  {
    TimeStamp ts;
    solver.correct(vec_cor, vec_def);
    time_total = ts.elapsed_now();
  }

  template<int dim_, typename DT_, typename IT_>
  void run_domain(Dist::Comm& comm, const std::array<Index, dim_>& base_size, const ParamHolder& param_holder,
    std::array<double, Times::count>& times, std::array<Index, Sizes::count>& sizes, PreferredBackend backend)
  {
    static_assert(dim_ == 2 || dim_ == 3, "richardson benchmark supports 2D and 3D domains");
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    typedef Control::Domain::PartiDomainControl<DomainLevelType> DomainControl;

    typedef DT_ DataType;
    typedef IT_ IndexType;
    typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, DataType, IndexType> SystemLevelType;

    String name_conv = "Dom";
    for(int i = 0; i < dim_; ++i)
      name_conv += "_" + stringify(base_size[std::size_t(i)]);
    FEAT::Util::ScopedNvtxRange _dom_nvtx((name_conv + "domain").c_str());

    // create domain
    DomainControl domain(comm, false);
    domain.set_desired_levels(0, 0);
    if constexpr(dim_ == 2)
      domain.create_rectilinear(base_size[0], DataType(0), DataType(1), base_size[1], DataType(0), DataType(1));
    else
      domain.create_rectilinear(base_size[0], DataType(0), DataType(1), base_size[1], DataType(0), DataType(1), base_size[2], DataType(0), DataType(1));
    const String cubature("gauss-legendre:" + stringify(SpaceVeloType::local_degree+1));

    std::shared_ptr<SystemLevelType> system;
    {
      system.reset(new SystemLevelType());
      TimeStamp ts;
      domain.front()->domain_asm.compile_all_elements();
      system->assemble_gates(domain.front());
      Assembly::Common::DuDvOperatorBlocked<dim_> dudv_operator;
      system->assemble_grad_div_matrices(domain.front()->domain_asm, domain.front()->space_velo, domain.front()->space_pres, cubature);
      system->assemble_velo_struct(domain.front()->space_velo);
      system->matrix_a.local().format();
      Assembly::assemble_bilinear_operator_matrix_1(domain.front()->domain_asm, system->matrix_a.local(), dudv_operator, domain.front()->space_velo, cubature);
      system->compile_system_matrix();
      times[Times::assembly] = ts.elapsed_now();
    }
    // assemble a Dirichlet velocity filter (pure-Dirichlet keeps the phase-1 device path exact)
    {
      Tiny::Vector<DataType, dim_> velo(DataType(0));
      for(int k = 0; k < dim_; ++k)
        velo[k] = 2.3 * k;
      Analytic::Common::ConstantVectorFunction<dim_> func(velo);
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      for(const char* part : {"bnd:t", "bnd:b", "bnd:l", "bnd:r", "bnd:f", "bnd:n"})
      {
        auto* mesh_part = domain.front()->get_mesh_node()->find_mesh_part(part);
        if(mesh_part)
          unit_asm.add_mesh_part(*mesh_part);
      }
      unit_asm.assemble(system->filter_velo.local(), domain.front()->space_velo, func);
    }
    system->compile_system_filter();
    system->compile_local_matrix_sys_type1();

    {
      auto tv = system->gate_velo._freqs.clone(LAFEM::CloneMode::Deep);
      tv.format(1.0);
      Index velo_dofs = Index(system->gate_velo.dot(tv, tv));
      Index locp_dofs = system->gate_pres._freqs.size();
      Index pres_dofs = Index(system->gate_pres.sum(DataType(locp_dofs)));
      sizes[Sizes::velo_dofs] = velo_dofs/dim_;
      sizes[Sizes::pres_dofs] = pres_dofs;
      sizes[Sizes::dofs] = velo_dofs+pres_dofs;
      sizes[Sizes::nnz] = system->matrix_sys.num_nzes_raw();
    }

    auto vec_cor = system->create_global_vector_sys();
    auto vec_def = system->create_global_vector_sys();
    vec_cor.format(DataType(1));
    system->filter_sys.filter_sol(vec_cor);
    vec_def.format();
    system->matrix_sys.apply(vec_def, vec_cor);
    system->filter_sys.filter_def(vec_def);
    system->filter_sys.filter_cor(vec_cor);

    // final defect norm ||vec_def - A*sol|| (filtered) of a computed solution
    auto final_defect = [&](const auto& sol) -> DataType {
      auto res = vec_def.clone();
      system->matrix_sys.apply(res, sol, vec_def, -DataType(1));
      system->filter_sys.filter_def(res);
      return res.norm2();
    };

    // the AmaVanka preconditioner (shared across the sequential solver runs)
    typedef typename SystemLevelType::LocalSystemMatrix LocalMatrix;
    typedef typename SystemLevelType::LocalSystemFilter LocalFilter;
    auto make_vanka = [&]() {
      auto v = std::make_shared<Solver::AmaVanka<LocalMatrix, LocalFilter>>(system->local_matrix_sys_type1, system->filter_sys.local());
      v->set_skip_singular(true);
      return v;
    };

    // standard Richardson(Schwarz(AmaVanka))
    if(param_holder.std_richardson)
    {
      auto tmp_cor = vec_cor.clone();
      auto vanka = make_vanka();
      auto schwarz = Solver::new_schwarz_precond(vanka, system->filter_sys);
      schwarz->set_ignore_status(true);
      auto rich = Solver::new_richardson(system->matrix_sys, system->filter_sys, DataType(param_holder.omega), schwarz);
      rich->set_plot_mode(Solver::PlotMode::none);
      rich->init();
      Backend::set_preferred_backend(backend);
      rich->set_min_iter(param_holder.num_iter);
      rich->set_max_iter(param_holder.num_iter);
      if(param_holder.hotrun)
        rich->correct(tmp_cor, vec_def);
      tmp_cor.copy(vec_cor);
      rich->set_plot_mode(param_holder.plot_mode);
      FEAT::Util::ScopedNvtxRange _nvtx((name_conv + "std_richardson").c_str());
      run_solver(*rich, tmp_cor, vec_def, times[Times::solve_std]);
      Backend::set_preferred_backend(PreferredBackend::generic);
      if(param_holder.print_defect)
        comm.print("  " + String("StdRich").pad_back(12) + " final |defect| = " + stringify_fp_sci(final_defect(tmp_cor)));
      rich->done();
    }

    // AmaVankaRichardson variants
    struct Variant { bool enabled; bool use_cuda; bool parallel; std::size_t slot; PreferredBackend be; const char* name; };
    const std::array<Variant, 3> variants{{
      {param_holder.generic_richardson, false, false, Times::solve_generic, backend, "NewGeneric"},
      {param_holder.cuda_richardson, true, false, Times::solve_cuda, backend, "NewCuda"},
      {param_holder.cuda_parallel_richardson, true, true, Times::solve_cuda_parallel, backend, "NewCudaPara"},
    }};
    for(const auto& var : variants)
    {
      if(!var.enabled)
        continue;
      auto tmp_cor = vec_cor.clone();
      auto vanka = make_vanka();
      auto rich = Solver::new_amavanka_richardson(system->matrix_sys, system->filter_sys, vanka, DataType(param_holder.omega));
      rich->set_use_cuda(var.use_cuda);
      rich->set_use_parallel_kernel(var.parallel);
      rich->set_plot_mode(Solver::PlotMode::none);
      rich->init();
      Backend::set_preferred_backend(var.be);
      rich->set_min_iter(param_holder.num_iter);
      rich->set_max_iter(param_holder.num_iter);
      if(param_holder.hotrun)
        rich->correct(tmp_cor, vec_def);
      tmp_cor.copy(vec_cor);
      rich->set_plot_mode(param_holder.plot_mode);
      FEAT::Util::ScopedNvtxRange _nvtx((name_conv + "amavanka_richardson").c_str());
      run_solver(*rich, tmp_cor, vec_def, times[var.slot]);
      Backend::set_preferred_backend(PreferredBackend::generic);
      if(param_holder.print_defect)
        comm.print("  " + String(var.name).pad_back(12) + " final |defect| = " + stringify_fp_sci(final_defect(tmp_cor)));
      rich->done();
    }
  }

  /// benchmark command-line interface, built on the typed ArgParser
  struct BenchArgs : public ArgParser
  {
    Parameter<Index> level = parameter(Index(1)).long_flag("--level")
      .help_string("number of size-growth levels");
    Parameter<Index> iters = parameter(Index(20)).long_flag("--iters")
      .help_string("fixed number of Richardson iterations per solve");
    Parameter<double> omega = parameter(double(0.5)).long_flag("--omega")
      .help_string("Richardson damping factor");
    Parameter<double> grow_factor = parameter(double(1.4)).long_flag("--grow-factor")
      .help_string("per-level, per-dimension cell-count growth factor");
    Parameter<std::deque<Index>> base_size = parameter(std::deque<Index>{1u, 1u, 1u}).long_flag("--base-size")
      .help_string("base cell counts; 2 or 3 values select the problem dimension");
    Parameter<String> backend = parameter(String("cuda")).long_flag("--backend")
      .help_string("compute backend: generic | mkl | cuda");
    Parameter<bool> hotrun = parameter(false).long_flag("--hotrun")
      .help_string("perform an untimed warm-up solve before timing");
    Parameter<bool> cuda_parallel = parameter(false).long_flag("--cuda-parallel")
      .help_string("also time the synchronized (parallel) CUDA kernel");
    Parameter<bool> defect = parameter(false).long_flag("--defect")
      .help_string("print the final defect norm after each solver run");
  };

  void main_test(int argc, char** argv)
  {
    Dist::Comm comm(Dist::Comm::world());

    BenchArgs args;
    args.set_description("Compare the standard Richardson(Schwarz(AmaVanka)) smoother against the fused AmaVankaRichardson.");
    if(!args.parse(argc, argv))
    {
      for(const auto& error : args.errors())
        comm.print(error);
      comm.print(args.help_text());
      return;
    }

    ParamHolder param_holder;
    param_holder.num_iter = args.iters.value();
    param_holder.omega = args.omega.value();
    param_holder.hotrun = args.hotrun.value();
    param_holder.cuda_parallel_richardson = args.cuda_parallel.value();
    param_holder.print_defect = args.defect.value();
    const Index level = args.level.value();
    const double grow_factor = args.grow_factor.value();

    PreferredBackend backend = PreferredBackend::cuda;
    {
      const String& backend_s = args.backend.value();
      if(backend_s.compare_no_case("generic") == 0)
        backend = PreferredBackend::generic;
      else if(backend_s.compare_no_case("mkl") == 0)
        backend = PreferredBackend::mkl;
      else if(backend_s.compare_no_case("cuda") == 0)
        backend = PreferredBackend::cuda;
      else
        XABORTM("Backend " + backend_s + " not known");
      comm.print("Using backend: " + backend_s);
    }
    if(backend != PreferredBackend::cuda)
    {
      param_holder.cuda_richardson = false;
      param_holder.cuda_parallel_richardson = false;
    }

    const std::deque<Index>& base_dq = args.base_size.value();
    const int problem_dim = int(base_dq.size());
    XASSERTM(problem_dim == 2 || problem_dim == 3, "--base-size expects either 2 or 3 entries");
    Index base_size[3] = {1u, 1u, 1u};
    for(int i = 0; i < problem_dim; ++i)
      base_size[i] = base_dq[std::size_t(i)];

    String msg = String("Number Processes").pad_back(pad_len, pad_char) + ":" + stringify(comm.size()) + "\n";
    msg += String("Level").pad_back(pad_len, pad_char) + ":" + stringify(level) + "\n";
    msg += String("Grow Factor").pad_back(pad_len, pad_char) + ":" + stringify(grow_factor) + "\n";
    msg += String("Richardson Iterations").pad_back(pad_len, pad_char) + ":" + stringify(param_holder.num_iter) + "\n";
    msg += String("Omega").pad_back(pad_len, pad_char) + ":" + stringify(param_holder.omega) + "\n";
    msg += String("Hotrun").pad_back(pad_len, pad_char) + ":" + (param_holder.hotrun ? String("Yes") : String("No")) + "\n";
    msg += String("Dimension").pad_back(pad_len, pad_char) + ":" + stringify(problem_dim) + "\n";
    comm.print(msg);
    comm.print_flush();

    std::vector<std::array<double, Times::count>> timing_array(std::size_t(level) * std::size_t(problem_dim));
    std::vector<std::array<Index, Sizes::count>> size_array(std::size_t(level) * std::size_t(problem_dim));

    for(Index k = 0; k < level; ++k)
    {
      for(int d = 0; d < problem_dim; ++d)
      {
        const Index idx = k*Index(problem_dim) + Index(d);
        size_array.at(idx)[Sizes::number_elements] = std::accumulate(&base_size[0], &base_size[0] + problem_dim, Index(1), [](Index a, Index b){return a*b;});
        if(problem_dim == 2)
        {
          std::array<Index, 2> run_size{{base_size[0], base_size[1]}};
          run_domain<2, double, std::uint32_t>(comm, run_size, param_holder, timing_array.at(idx), size_array.at(idx), backend);
        }
        else
        {
          std::array<Index, 3> run_size{{base_size[0], base_size[1], base_size[2]}};
          run_domain<3, double, std::uint32_t>(comm, run_size, param_holder, timing_array.at(idx), size_array.at(idx), backend);
        }
        base_size[d] = Index(std::ceil(double(base_size[d]) * grow_factor));
      }
    }
    comm.print("Done");

    constexpr int tp = 15;
    msg = String("cells").pad_back(tp) + String("dofs").pad_back(tp) + String("nnz").pad_back(tp);
    if(param_holder.std_richardson)          msg += String("StdRich").pad_back(tp);
    if(param_holder.generic_richardson)      msg += String("NewGeneric").pad_back(tp);
    if(param_holder.cuda_richardson)         msg += String("NewCuda").pad_back(tp);
    if(param_holder.cuda_parallel_richardson) msg += String("NewCudaPara").pad_back(tp);
    msg += "\n";
    for(std::size_t k = 0; k < timing_array.size(); ++k)
    {
      msg += stringify(size_array.at(k)[Sizes::number_elements]).pad_back(tp);
      msg += stringify(size_array.at(k)[Sizes::dofs]).pad_back(tp);
      msg += stringify(size_array.at(k)[Sizes::nnz]).pad_back(tp);
      if(param_holder.std_richardson)          msg += stringify_fp_fix(timing_array.at(k)[Times::solve_std]).pad_back(tp);
      if(param_holder.generic_richardson)      msg += stringify_fp_fix(timing_array.at(k)[Times::solve_generic]).pad_back(tp);
      if(param_holder.cuda_richardson)         msg += stringify_fp_fix(timing_array.at(k)[Times::solve_cuda]).pad_back(tp);
      if(param_holder.cuda_parallel_richardson) msg += stringify_fp_fix(timing_array.at(k)[Times::solve_cuda_parallel]).pad_back(tp);
      msg += "\n";
    }
    comm.print("Solve wall-clock [s]\n------------------------");
    comm.print(msg);
  }
}

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard guard(argc, argv);
  RichardsonBench::main_test(argc, argv);
}
