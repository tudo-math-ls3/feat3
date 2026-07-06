#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/assembly/burgers_assembly_job.hpp>
#include <kernel/assembly/bilinear_operator.hpp>
#include <kernel/assembly/grad_operator_assembler.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/domain_assembler.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <control/stokes_blocked.hpp>
#include <control/domain/domain_level.hpp>
#include <control/domain/parti_domain_control.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/amavanka_gmres.hpp>

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#else
namespace FEAT::Util{
    inline std::uint64_t cuda_start_range(const char*)
    {
     return 0;
    }
    inline void cuda_stop_range(std::uint64_t)
    {
    }
    inline void cuda_mark_event(const char*)
    {}
  }
#endif

namespace GMRESTest
{

  using namespace FEAT;

  constexpr std::size_t pad_len = 20;
  constexpr char pad_char = '.';

  struct ParamHolder
  {
    Index num_iter = 500;
    Index krylov_dim = 200;
    Solver::PlotMode plot_mode = Solver::PlotMode::summary;
    bool base_gmres = true;
    bool ortho_gmres = true;
    bool dot_gmres = true;
    bool scalar_gmres = true;
    bool cuda_gmres = true;
    bool hotrun = false;
    bool cuda_prefetch = false;
    bool compare_to_parallel_kernel = false;
  };

  struct Times
  {
    static constexpr std::size_t total_time = 0u;
    static constexpr std::size_t assembly = 1u;
    static constexpr std::size_t solve_stnd_gmres = 2u;
    static constexpr std::size_t solve_dot_gmres = 3u;
    static constexpr std::size_t solve_dot_gmres_alt = 4u;
    static constexpr std::size_t solve_scalar_gmres = 5u;
    static constexpr std::size_t solve_cuda_gmres = 6u;
    static constexpr std::size_t solve_cuda_gmres_prefetch = 7u;
    static constexpr std::size_t solve_cuda_gmres_parallel_kernel = 8u;
    /// standard solver/gmres.hpp with the batched classical Gram-Schmidt orthogonalizer
    static constexpr std::size_t solve_stnd_gmres_cgs = 9u;
    /// standard solver/gmres.hpp with the batched CGS2 (reorthogonalizing) orthogonalizer
    static constexpr std::size_t solve_stnd_gmres_cgs2 = 10u;
    /// scalar (CPU) AmaVankaGMRES with the CGS2 (reorthogonalizing) orthogonalizer
    static constexpr std::size_t solve_scalar_gmres_cgs2 = 11u;
    /// cuda AmaVankaGMRES with the CGS2 (reorthogonalizing) orthogonalizer
    static constexpr std::size_t solve_cuda_gmres_cgs2 = 12u;
    static constexpr std::size_t count = 14u;
  };

  struct Sizes
  {
    static constexpr std::size_t number_elements = 0u;
    static constexpr std::size_t nnz = 1u;
    static constexpr std::size_t velo_dofs = 2u;
    static constexpr std::size_t pres_dofs = 3u;
    static constexpr std::size_t dofs = 4u;
    static constexpr std::size_t flops = 5u;
    static constexpr std::size_t count = 10u;
  };

  template<typename Solver_, typename VecType_>
  void run_solver(Solver_& solver, VecType_& vec_cor, VecType_& vec_def, double& time_total)
  {
    TimeStamp ts;
    solver.correct(vec_cor, vec_def);
    time_total = ts.elapsed_now();
  }

  template<typename Solver_, typename VecType_>
  void hotrun_solver(Solver_& solver, VecType_& vec_cor, VecType_& vec_def)
  {
    solver.correct(vec_cor, vec_def);
  }

  template<int dim_, typename DT_, typename IT_>
  void run_domain(Dist::Comm& comm, const std::array<Index, dim_>& base_size, const ParamHolder& param_holder, std::array<double, Times::count>& times, std::array<Index, Sizes::count>& sizes, PreferredBackend backend)
  {
    static_assert(dim_ == 2 || dim_ == 3, "gpu gmres benchmark supports 2D and 3D domains");
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::ElementP1<TrafoType> SpacePresType;
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> DomainLevelType;
    typedef Control::Domain::PartiDomainControl<DomainLevelType> DomainControl;

    String name_conv = "Dom";
    for(int i = 0; i < dim_; ++i)
      name_conv += "_" + stringify(base_size[std::size_t(i)]);
    String dom_name = name_conv + "domain";
    auto dom_id = Util::cuda_start_range(dom_name.c_str());
    typedef DT_ DataType;
    typedef IT_ IndexType;
    typedef Control::StokesBlockedUnitVeloNonePresSystemLevel<dim_, DataType, IndexType> SystemLevelType;
    // create domain
    DomainControl domain(comm, false);
    domain.set_desired_levels(0, 0);
    if constexpr(dim_ == 2)
    {
      domain.create_rectilinear(base_size[0], DataType(0), DataType(1), base_size[1], DataType(0), DataType(1));
    }
    else
    {
      domain.create_rectilinear(base_size[0], DataType(0), DataType(1), base_size[1], DataType(0), DataType(1), base_size[2], DataType(0), DataType(1));
    }
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
    // assemble filter
    {
      TimeStamp ts;
      Tiny::Vector<DataType, dim_> velo(DataType(0));
      for(int k = 0; k < dim_; ++k)
        velo[k] = 2.3 * k;
      Analytic::Common::ConstantVectorFunction<dim_> func(velo);
      Assembly::UnitFilterAssembler<MeshType> unit_asm;
      {
        auto* mesh_part = domain.front()->get_mesh_node()->find_mesh_part("bnd:t");
        if(mesh_part)
          unit_asm.add_mesh_part(*mesh_part);
      }
      {
        auto* mesh_part = domain.front()->get_mesh_node()->find_mesh_part("bnd:l");
        if(mesh_part)
          unit_asm.add_mesh_part(*mesh_part);
      }
      {
        auto* mesh_part = domain.front()->get_mesh_node()->find_mesh_part("bnd:r");
        if(mesh_part)
          unit_asm.add_mesh_part(*mesh_part);
      }
      unit_asm.assemble(system->filter_velo.local(), domain.front()->space_velo, func);

    }
    system->compile_system_filter();
    system->compile_local_matrix_sys_type1();

    {
      // collect some finest-level statistics
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
        sizes[Sizes::flops] = Index(0);
      }

    }

    auto vec_cor = system->create_global_vector_sys();
    auto vec_def = system->create_global_vector_sys();
    vec_cor.format(DataType(1));
    system->filter_sys.filter_sol(vec_cor);
    vec_def.format();
    // calc defect for 1 vector
    system->matrix_sys.apply(vec_def, vec_cor);
    system->filter_sys.filter_def(vec_def);
    system->filter_sys.filter_cor(vec_cor);

    // build vanka matrix
    std::shared_ptr<Solver::AmaVanka<typename SystemLevelType::LocalSystemMatrix, typename SystemLevelType::LocalSystemFilter>> vanka(
       new Solver::AmaVanka<typename SystemLevelType::LocalSystemMatrix, typename SystemLevelType::LocalSystemFilter>(system->local_matrix_sys_type1, system->filter_sys.local()));
    // vanka.init();
    if(param_holder.base_gmres)
    {
      auto tmp_def = vec_def.clone();
      auto tmp_cor = vec_cor.clone();
      // schwarz global to local layer
      auto schwarz = Solver::new_schwarz_precond(vanka, system->filter_sys);
      auto gmres = Solver::new_gmres(system->matrix_sys, system->filter_sys, param_holder.krylov_dim, DataType(1), schwarz);
      gmres->set_plot_mode(Solver::PlotMode::none);

      gmres->init();
      Backend::set_preferred_backend(backend);
      gmres->set_min_iter(param_holder.krylov_dim + 1);
      gmres->set_max_iter(param_holder.krylov_dim + 1);
      if(param_holder.hotrun)
        hotrun_solver(*gmres, tmp_cor, tmp_def);

      gmres->set_plot_mode(param_holder.plot_mode);
      gmres->set_min_iter(param_holder.num_iter);
      gmres->set_max_iter(param_holder.num_iter);
      String name = name_conv + "stnd_gmres";
      auto id = Util::cuda_start_range(name.c_str());
      run_solver(*gmres, tmp_cor, tmp_def, times[Times::solve_stnd_gmres]);
      Util::cuda_stop_range(id);
      Backend::set_preferred_backend(PreferredBackend::generic);
      gmres->done();
    }
    if(param_holder.ortho_gmres)
    {
      // standard solver/gmres.hpp driven with the batched orthogonalizers (classical Gram-Schmidt
      // and CGS2), reproducing the POC's DotGMRES behavior via the KrylovOrthogonalizer interface
      const std::array<std::pair<Solver::KrylovOrthogonalizerType, std::size_t>, 2> variants{{
        {Solver::KrylovOrthogonalizerType::classical_gram_schmidt, Times::solve_stnd_gmres_cgs},
        {Solver::KrylovOrthogonalizerType::cgs2, Times::solve_stnd_gmres_cgs2}
      }};
      for(const auto& variant : variants)
      {
        auto tmp_def = vec_def.clone();
        auto tmp_cor = vec_cor.clone();
        // schwarz global to local layer
        auto schwarz = Solver::new_schwarz_precond(vanka, system->filter_sys);
        auto gmres = Solver::new_gmres(system->matrix_sys, system->filter_sys, param_holder.krylov_dim, DataType(1), schwarz);
        gmres->set_orthogonalizer(variant.first);
        gmres->set_plot_mode(Solver::PlotMode::none);

        gmres->init();
        Backend::set_preferred_backend(backend);
        gmres->set_min_iter(param_holder.krylov_dim + 1);
        gmres->set_max_iter(param_holder.krylov_dim + 1);
        if(param_holder.hotrun)
          hotrun_solver(*gmres, tmp_cor, tmp_def);

        gmres->set_plot_mode(param_holder.plot_mode);
        gmres->set_min_iter(param_holder.num_iter);
        gmres->set_max_iter(param_holder.num_iter);
        String name = name_conv + (variant.first == Solver::KrylovOrthogonalizerType::cgs2 ? "stnd_gmres_cgs2" : "stnd_gmres_cgs");
        auto id = Util::cuda_start_range(name.c_str());
        run_solver(*gmres, tmp_cor, tmp_def, times[variant.second]);
        Util::cuda_stop_range(id);
        Backend::set_preferred_backend(PreferredBackend::generic);
        gmres->done();
      }
    }
    if(param_holder.scalar_gmres)
    {
      // scalar (CPU) AmaVankaGMRES with the classical Gram-Schmidt and CGS2 orthogonalizers
      const std::array<std::pair<Solver::KrylovOrthogonalizerType, std::size_t>, 2> variants{{
        {Solver::KrylovOrthogonalizerType::classical_gram_schmidt, Times::solve_scalar_gmres},
        {Solver::KrylovOrthogonalizerType::cgs2, Times::solve_scalar_gmres_cgs2}
      }};
      for(const auto& variant : variants)
      {
        auto tmp_def = vec_def.clone();
        auto tmp_cor = vec_cor.clone();
        // schwarz global to local layer
        // auto schwarz = Solver::new_schwarz_precond(vanka, system->filter_sys);
        auto gmres = Solver::new_amavanka_gmres(system->matrix_sys, system->filter_sys, vanka, param_holder.krylov_dim, DataType(1));
        gmres->set_plot_mode(Solver::PlotMode::none);
        gmres->set_use_cuda_gmres(false);
        gmres->set_orthogonalizer(variant.first);

        gmres->init();
        Backend::set_preferred_backend(backend);
        gmres->set_min_iter(param_holder.krylov_dim + 1);
        gmres->set_max_iter(param_holder.krylov_dim + 1);
        if(param_holder.hotrun)
          hotrun_solver(*gmres, tmp_cor, tmp_def);
        gmres->set_plot_mode(param_holder.plot_mode);
        gmres->set_min_iter(param_holder.num_iter);
        gmres->set_max_iter(param_holder.num_iter);
        String name = name_conv + (variant.first == Solver::KrylovOrthogonalizerType::cgs2 ? "scalar gmres cgs2" : "scalar gmres");
        auto id = Util::cuda_start_range(name.c_str());
        run_solver(*gmres, tmp_cor, tmp_def, times[variant.second]);
        // std::printf("Flops %lu \n", gmres->get_flops());
        // report the flop baseline from the plain CGS run
        if(variant.first == Solver::KrylovOrthogonalizerType::classical_gram_schmidt)
          sizes[Sizes::flops] = gmres->get_global_flops();
        Util::cuda_stop_range(id);
        Backend::set_preferred_backend(PreferredBackend::generic);
        gmres->done();
      }
    }
    if(param_holder.cuda_gmres)
    {
      // cuda AmaVankaGMRES with the classical Gram-Schmidt and CGS2 orthogonalizers
      const std::array<std::pair<Solver::KrylovOrthogonalizerType, std::size_t>, 2> variants{{
        {Solver::KrylovOrthogonalizerType::classical_gram_schmidt, Times::solve_cuda_gmres},
        {Solver::KrylovOrthogonalizerType::cgs2, Times::solve_cuda_gmres_cgs2}
      }};
      for(const auto& variant : variants)
      {
        auto tmp_def = vec_def.clone();
        auto tmp_cor = vec_cor.clone();
        // schwarz global to local layer
        // auto schwarz = Solver::new_schwarz_precond(vanka, system->filter_sys);
        auto gmres = Solver::new_amavanka_gmres(system->matrix_sys, system->filter_sys, vanka, param_holder.krylov_dim, DataType(1));
        gmres->set_plot_mode(Solver::PlotMode::none);
        gmres->set_use_cuda_gmres(true);
        gmres->set_orthogonalizer(variant.first);

        gmres->init();
        Backend::set_preferred_backend(backend);
        gmres->set_min_iter(param_holder.krylov_dim + 1);
        gmres->set_max_iter(param_holder.krylov_dim + 1);
        if(param_holder.hotrun)
          hotrun_solver(*gmres, tmp_cor, tmp_def);
        gmres->set_plot_mode(param_holder.plot_mode);
        gmres->set_min_iter(param_holder.num_iter);
        gmres->set_max_iter(param_holder.num_iter);
        String name = name_conv + (variant.first == Solver::KrylovOrthogonalizerType::cgs2 ? "cuda gmres cgs2" : "cuda gmres");
        auto id = Util::cuda_start_range(name.c_str());
        run_solver(*gmres, tmp_cor, tmp_def, times[variant.second]);
        Util::cuda_stop_range(id);
        Backend::set_preferred_backend(PreferredBackend::generic);
        gmres->done();
      }
      if(param_holder.compare_to_parallel_kernel)
      {
        auto tmp_def = vec_def.clone();
        auto tmp_cor = vec_cor.clone();
        // schwarz global to local layer
        // auto schwarz = Solver::new_schwarz_precond(vanka, system->filter_sys);
        auto gmres = Solver::new_amavanka_gmres(system->matrix_sys, system->filter_sys, vanka, param_holder.krylov_dim, DataType(1));
        gmres->set_plot_mode(Solver::PlotMode::none);
        gmres->set_use_cuda_gmres(true);
        gmres->set_use_parallel_kernel(true);

        gmres->init();
        Backend::set_preferred_backend(backend);
        gmres->set_min_iter(param_holder.krylov_dim + 1);
        gmres->set_max_iter(param_holder.krylov_dim + 1);
        if(param_holder.hotrun)
          hotrun_solver(*gmres, tmp_cor, tmp_def);
        gmres->set_plot_mode(param_holder.plot_mode);
        gmres->set_min_iter(param_holder.num_iter);
        gmres->set_max_iter(param_holder.num_iter);
        String name = name_conv + "cuda gmres";
        auto id = Util::cuda_start_range(name.c_str());
        run_solver(*gmres, tmp_cor, tmp_def, times[Times::solve_cuda_gmres_parallel_kernel]);
        Util::cuda_stop_range(id);
        Backend::set_preferred_backend(PreferredBackend::generic);
        gmres->done();

      }
      if(param_holder.cuda_prefetch)
      {
        auto tmp_def = vec_def.clone();
        auto tmp_cor = vec_cor.clone();
        // schwarz global to local layer
        // auto schwarz = Solver::new_schwarz_precond(vanka, system->filter_sys);
        auto gmres = Solver::new_amavanka_gmres(system->matrix_sys, system->filter_sys, vanka, param_holder.krylov_dim, DataType(1));
        gmres->set_plot_mode(Solver::PlotMode::none);
        gmres->set_use_cuda_gmres(true);
        gmres->set_prefetch_data(true);

        gmres->init();
        Backend::set_preferred_backend(backend);
        gmres->set_min_iter(param_holder.krylov_dim + 1);
        gmres->set_max_iter(param_holder.krylov_dim + 1);
        if(param_holder.hotrun)
          hotrun_solver(*gmres, tmp_cor, tmp_def);
        gmres->set_plot_mode(param_holder.plot_mode);
        gmres->set_min_iter(param_holder.num_iter);
        gmres->set_max_iter(param_holder.num_iter);
        String name = name_conv + "cuda gmres";
        auto id = Util::cuda_start_range(name.c_str());
        run_solver(*gmres, tmp_cor, tmp_def, times[Times::solve_cuda_gmres_prefetch]);
        Util::cuda_stop_range(id);
        Backend::set_preferred_backend(PreferredBackend::generic);
        gmres->done();
      }
    }


    Util::cuda_stop_range(dom_id);
  }




  void main_test(int argc, char** argv)
  {
    Dist::Comm comm(Dist::Comm::world());

    String msg = String("Number Proccesses").pad_back(pad_len, pad_char) + ":" + stringify(comm.size()) + "\n";

    SimpleArgParser args(argc, argv);
    ParamHolder param_holder;
    args.support("level");
    args.support("iters");
    args.support("krylov-dim");
    args.support("base-size");
    args.support("grow-factor");
    args.support("backend");
    args.support("hotrun");

    Index level = args.parse_default("level", Index(1u));
    double grow_factor = args.parse_default("grow-factor", double(1.4));
    param_holder.num_iter = args.parse_default("iters", Index(500u));
    param_holder.krylov_dim = args.parse_default("krylov-dim", Index(200u));
    PreferredBackend backend = PreferredBackend::cuda;
    {
      String backend_s = args.parse_default("backend", String("cuda"));
      if(backend_s.compare_no_case("generic") == 0)
      {
        backend = PreferredBackend::generic;
      }
      else if(backend_s.compare_no_case("mkl") == 0)
      {
        backend = PreferredBackend::mkl;
      }
      else if(backend_s.compare_no_case("cuda") == 0)
      {
        backend = PreferredBackend::cuda;
      }
      else
      {
        XABORTM("Backend " + backend_s + " not known");
      }
      comm.print("Using backend: " + backend_s);
    }
    if(args.check("hotrun") >= 0)
    {
      param_holder.hotrun = true;
    }

    if(args.check("plot-mode") > 0)
    {
      String plot_mode = args.parse_default("plot-mode", String("summary"));
      if(plot_mode.compare_no_case("summary") == 0)
      {
        param_holder.plot_mode = Solver::PlotMode::summary;
      }
      else if((plot_mode.compare_no_case("none") == 0) || (plot_mode.compare_no_case("off") == 0))
      {
        param_holder.plot_mode = Solver::PlotMode::none;
      }
      else if((plot_mode.compare_no_case("all") == 0))
      {
        param_holder.plot_mode = Solver::PlotMode::all;
      }
      else
      {
        XABORTM(String("Option ") + plot_mode + " not known");
      }
    }

    // if(comm.size() > 1 || backend != PreferredBackend::cuda)
    if(backend != PreferredBackend::cuda)
    {
      param_holder.cuda_gmres = false;
    }


    msg += String("Level").pad_back(pad_len, pad_char) + ":" + stringify(level) + "\n";
    msg += String("Grow Factor per level per dimension").pad_back(pad_len, pad_char) + ":" + stringify(grow_factor) + "\n";
    msg += String("GMRES Iterations").pad_back(pad_len, pad_char) + ":" + stringify(param_holder.num_iter) + "\n";
    msg += String("Krylov-Dim").pad_back(pad_len, pad_char) + ":" + stringify(param_holder.krylov_dim) + "\n";
    msg += String("Hotrun").pad_back(pad_len, pad_char) + ":" + (param_holder.hotrun ?String("Yes") : String("No")) + "\n";

    // setup levels
    int problem_dim = 3;
    Index base_size[3] = {1u, 1u, 1u};

    if(args.check("base-size") >= 0)
    {
      problem_dim = args.check("base-size");
      XASSERTM(problem_dim == 2 || problem_dim == 3, "--base-size expects either 2 or 3 entries");
      bool success = true;
      auto opt = args.query("base-size")->second;
      for(int i = 0; i < problem_dim; ++i)
        success &= opt.at(std::size_t(i)).parse(base_size[i]);
      XASSERTM(success, "Something went wrong while parsing base sizes");
    }

    msg += String("Dimension").pad_back(pad_len, pad_char) + ":" + stringify(problem_dim) + "\n";
    msg += String("Base sizes").pad_back(pad_len, pad_char) + ":";
    for(int i = 0; i < problem_dim; ++i)
    {
      if(i > 0)
        msg += ", ";
      msg += stringify(base_size[i]);
    }
    msg += "\n";
    comm.print(msg);
    comm.print_flush();

    std::vector<std::array<double, Times::count>> timing_array;
    std::vector<std::array<Index, Sizes::count>> size_array;
    timing_array.resize((level) * Index(problem_dim));
    size_array.resize((level) * Index(problem_dim));


    for(Index k = 0; k < level; ++k)
    {
      for(int d = 0; d < problem_dim; ++d)
      {
        msg = String("Running size").pad_back(pad_len, pad_char) + ":";
        for(int i = 0; i < problem_dim; ++i)
        {
          if(i > 0)
            msg += " x ";
          msg += stringify(base_size[i]);
        }
        msg += "\n";
        comm.print(msg);
        comm.print_flush();

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

    constexpr int tabelle_pad = 15;

    // print stnd gmres timings
    {
      msg =String("cells").pad_back(tabelle_pad);
      msg +=String("dofs").pad_back(tabelle_pad);
      msg += String("nnz").pad_back(tabelle_pad);
      if(param_holder.base_gmres)
        msg += String("Stnd").pad_back(15);
      if(param_holder.ortho_gmres)
      {
        msg += String("Stnd-CGS").pad_back(15);
        msg += String("Stnd-CGS2").pad_back(15);
      }
      if(param_holder.scalar_gmres)
      {
        msg += String("Scalar").pad_back(15);
        msg += String("Scalar-CGS2").pad_back(15);
      }
      if(param_holder.cuda_gmres)
      {
        msg += String("Cuda").pad_back(15);
        msg += String("Cuda-CGS2").pad_back(15);
        if(param_holder.compare_to_parallel_kernel)
          msg += String("CudaPara").pad_back(15);
        if(param_holder.cuda_prefetch)
          msg += String("CudaPre").pad_back(15);
      }

      msg += "\n";
      for(std::size_t k = 0; k < timing_array.size(); ++k)
      {
        msg += stringify(size_array.at(k)[Sizes::number_elements]).pad_back(tabelle_pad);
        msg += stringify(size_array.at(k)[Sizes::dofs]).pad_back(tabelle_pad);
        msg += stringify(size_array.at(k)[Sizes::nnz]).pad_back(tabelle_pad);
        if(param_holder.base_gmres)
          msg += stringify_fp_fix(timing_array.at(k)[Times::solve_stnd_gmres]).pad_back(tabelle_pad);
        if(param_holder.ortho_gmres)
        {
          msg += stringify_fp_fix(timing_array.at(k)[Times::solve_stnd_gmres_cgs]).pad_back(tabelle_pad);
          msg += stringify_fp_fix(timing_array.at(k)[Times::solve_stnd_gmres_cgs2]).pad_back(tabelle_pad);
        }
        if(param_holder.scalar_gmres)
        {
          msg += stringify_fp_fix(timing_array.at(k)[Times::solve_scalar_gmres]).pad_back(tabelle_pad);
          msg += stringify_fp_fix(timing_array.at(k)[Times::solve_scalar_gmres_cgs2]).pad_back(tabelle_pad);
        }
        if(param_holder.cuda_gmres)
        {
          msg += stringify_fp_fix(timing_array.at(k)[Times::solve_cuda_gmres]).pad_back(tabelle_pad);
          msg += stringify_fp_fix(timing_array.at(k)[Times::solve_cuda_gmres_cgs2]).pad_back(tabelle_pad);
          if(param_holder.compare_to_parallel_kernel)
            msg += stringify_fp_fix(timing_array.at(k)[Times::solve_cuda_gmres_parallel_kernel]).pad_back(tabelle_pad);
          if(param_holder.cuda_prefetch)
            msg += stringify_fp_fix(timing_array.at(k)[Times::solve_cuda_gmres_prefetch]).pad_back(tabelle_pad);
        }

        msg += "\n";
      }
    }
    comm.print(msg);

    // print flop, if available
    if(size_array.at(0)[Sizes::flops] > 0u)
    {
      // reduce
      msg =String("cells").pad_back(tabelle_pad);
      msg +=String("dofs").pad_back(tabelle_pad);
      msg += String("nnz").pad_back(tabelle_pad);
      if(param_holder.base_gmres)
        msg += String("Stnd").pad_back(15);
      if(param_holder.ortho_gmres)
      {
        msg += String("Stnd-CGS").pad_back(15);
        msg += String("Stnd-CGS2").pad_back(15);
      }
      if(param_holder.scalar_gmres)
      {
        msg += String("Scalar").pad_back(15);
        msg += String("Scalar-CGS2").pad_back(15);
      }
      if(param_holder.cuda_gmres)
      {
        msg += String("Cuda").pad_back(15);
        msg += String("Cuda-CGS2").pad_back(15);
        if(param_holder.compare_to_parallel_kernel)
          msg += String("CudaPara").pad_back(15);
        if(param_holder.cuda_prefetch)
          msg += String("CudaPre").pad_back(15);
      }

      auto calc_flops = [](double time, Index flops){return (double(flops)/1E+9)/time;};

      msg += "\n";
      for(std::size_t k = 0; k < timing_array.size(); ++k)
      {
        msg += stringify(size_array.at(k)[Sizes::number_elements]).pad_back(tabelle_pad);
        msg += stringify(size_array.at(k)[Sizes::dofs]).pad_back(tabelle_pad);
        msg += stringify(size_array.at(k)[Sizes::nnz]).pad_back(tabelle_pad);
        if(param_holder.base_gmres)
          msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_stnd_gmres], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
        if(param_holder.ortho_gmres)
        {
          msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_stnd_gmres_cgs], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
          msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_stnd_gmres_cgs2], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
        }
        if(param_holder.scalar_gmres)
        {
          msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_scalar_gmres], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
          msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_scalar_gmres_cgs2], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
        }
        if(param_holder.cuda_gmres)
        {
          msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_cuda_gmres], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
          msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_cuda_gmres_cgs2], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
          if(param_holder.compare_to_parallel_kernel)
            msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_cuda_gmres_parallel_kernel], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
          if(param_holder.cuda_prefetch)
            msg += stringify_fp_fix(calc_flops(timing_array.at(k)[Times::solve_cuda_gmres_prefetch], size_array.at(k)[Sizes::flops])).pad_back(tabelle_pad);
        }

        msg += "\n";
      }
      comm.print("GFlops/s\n------------------------\n");
      comm.print(msg);
    }

  }

}











int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard guard(argc, argv);

  GMRESTest::main_test(argc, argv);

}
