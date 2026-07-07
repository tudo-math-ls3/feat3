// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/backend.hpp>

#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/util/arg_parser.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/voxel_amavanka.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/gmres.hpp>
#include <kernel/solver/amavanka_gmres.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/richardson.hpp>
#include "logger.hpp"
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/string.hpp>

#include <kernel/solver/direct_sparse_solver.hpp>
#include <kernel/solver/frosch.hpp>
#include <control/scalar_basic.hpp>
#include "format_helper.hpp"
#include "template_helper.hpp"
#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/coloring.hpp>

#include <memory>

namespace Gendie
{
  namespace Intern
  {
    inline bool parse_config_bool(const FEAT::String& option)
    {
      return !(option.empty() || (option.compare_no_case("false") == 0) || (option.compare_no_case("off") == 0) || (option.compare_no_case("0") == 0));
    }

    inline bool is_backend(const FEAT::String& option)
    {
      return (option.compare_no_case("cuda") == 0) || (option.compare_no_case("mkl") == 0) || (option.compare_no_case("generic") == 0);
    }

    inline FEAT::PreferredBackend backend_from_string(const FEAT::String& option)
    {
      if(option.compare_no_case("cuda") == 0)
        return FEAT::PreferredBackend::cuda;
      if(option.compare_no_case("mkl") == 0)
        return FEAT::PreferredBackend::mkl;
      if(option.compare_no_case("generic") == 0)
        return FEAT::PreferredBackend::generic;
      XABORTM("Unknown backend " + option);
      return FEAT::PreferredBackend::generic;
    }

    inline bool is_direct_solver(const FEAT::String& option)
    {
      return (option.compare_no_case("umfpack") == 0) || (option.compare_no_case("mumps") == 0) || (option.compare_no_case("ilu") == 0) || (option.compare_no_case("klu") == 0);
    }

    inline bool is_orthogonalizer(const FEAT::String& option)
    {
      return (option.compare_no_case("mgs") == 0 || option.compare_no_case("cgs") == 0) || (option.compare_no_case("cgs2") == 0);
    }

    inline FEAT::Solver::KrylovOrthogonalizerType orthogonalizer_from_string(const FEAT::String& option)
    {
      if(option.compare_no_case("mgs") == 0)
        return FEAT::Solver::KrylovOrthogonalizerType::modified_gram_schmidt;
      if(option.compare_no_case("cgs") == 0)
        return FEAT::Solver::KrylovOrthogonalizerType::classical_gram_schmidt;
      if(option.compare_no_case("cgs2") == 0)
        return FEAT::Solver::KrylovOrthogonalizerType::cgs2;
      XABORTM("Unknown orthogonalizer " + option);
      return FEAT::Solver::KrylovOrthogonalizerType::cgs2;
    }

    inline bool is_frosch_precond(const FEAT::String& option)
    {
      return (option.compare_no_case("onelevel") == 0) || (option.compare_no_case("twolevel") == 0) || (option.compare_no_case("twolevelblock") == 0);
    }

    inline bool is_frosch_partition_type(const FEAT::String& option)
    {
      return (option.compare_no_case("PHG") == 0) || (option.compare_no_case("PARMETIS") == 0) || (option.compare_no_case("SCOTCH") == 0) || (option.compare_no_case("PTSCOTCH") == 0) || (option.compare_no_case("BLOCK") == 0);
    }

    inline bool is_frosch_partition_approach(const FEAT::String& option)
    {
      return (option.compare_no_case("partition") == 0) || (option.compare_no_case("repartition") == 0);
    }

    inline bool is_frosch_combine_overlap(const FEAT::String& option)
    {
      return (option.compare_no_case("restricted") == 0) || (option.compare_no_case("full") == 0) || (option.compare_no_case("averaging") == 0);
    }

    inline bool is_frosch_combine_level(const FEAT::String& option)
    {
      return (option.compare_no_case("additive") == 0) || (option.compare_no_case("multiplicative") == 0);
    }

    inline bool is_frosch_ipou(const FEAT::String& option)
    {
      return (option.compare_no_case("gdsw") == 0) || (option.compare_no_case("gdswstar") == 0) || (option.compare_no_case("rgdsw") == 0);
    }

    template<typename DT_>
    struct SolverParser : public FEAT::ArgParser
    {
      using DataType = DT_;

      FEAT::Parameter<FEAT::String> nl_solver = parameter(FEAT::String("alpine"))
                                                  .property("nl-solver")
                                                  .name("Nonlinear Solver")
                                                  .validator([](const auto& s) { return s.compare_no_case("alpine") == 0 || s.compare_no_case("backtrace-newton") == 0 || s.compare_no_case("pseudo-ts") == 0 || s.compare_no_case("newton") == 0 || s.compare_no_case("picard") == 0; });
      FEAT::Parameter<FEAT::String> backend = parameter(FEAT::String("generic"))
                                                .property("backend")
                                                .name("Solver Backend")
                                                .validator(is_backend);

      FEAT::Parameter<FEAT::Index> min_mg_iter = parameter(FEAT::Index(1)).property("min-mg-iter").name("Min MG Iter");
      FEAT::Parameter<FEAT::Index> max_mg_iter = parameter(FEAT::Index(50)).property("max-mg-iter").name("Max MG Iter");
      FEAT::Parameter<FEAT::Index> smooth_gmres = parameter(FEAT::Index(0)).property("smooth-gmres").name("Smoother GMRES Dimension");
      FEAT::Parameter<FEAT::Index> solve_gmres = parameter(FEAT::Index(4)).property("solve-gmres").name("Solver GMRES Dimension");
      FEAT::Parameter<FEAT::Index> smooth_steps = parameter(FEAT::Index(32)).property("smooth-steps").name("Smoother Steps");
      FEAT::Parameter<DataType> smooth_damp = parameter(DataType(0.5)).property("smooth-damp").name("Smoother Damping");
      FEAT::Parameter<DataType> mg_tol_rel = parameter(DataType(0)).property("mg-tol-rel").name("MG Relative Tolerance");
      FEAT::Parameter<bool> coarse_gmres = parameter(true).property("coarse-gmres").name("Use Coarse GMRES").parser(parse_config_bool);
      FEAT::Parameter<bool> coarse_frosch = parameter(false).property("coarse-frosch").name("Use Coarse FROSch").parser(parse_config_bool);
      FEAT::Parameter<bool> coarse_info = parameter(false).property("coarse-info").name("Print Coarse Solver Info").parser(parse_config_bool);

      FEAT::Parameter<FEAT::Index> coarse_gmres_dim = parameter(FEAT::Index(200)).property("coarse-gmres-parameters/gmres-dim").name("Coarse GMRES Dimension");
      FEAT::Parameter<FEAT::Index> coarse_min_iter = parameter(FEAT::Index(1)).property("coarse-gmres-parameters/min-iter").name("Coarse Min Iter");
      FEAT::Parameter<FEAT::Index> coarse_max_iter = parameter(FEAT::Index(1000)).property("coarse-gmres-parameters/max-iter").name("Coarse Max Iter");
      FEAT::Parameter<DataType> coarse_tol_rel = parameter(DataType(1E-2)).property("coarse-gmres-parameters/tol-rel").name("Coarse Relative Tolerance");
      FEAT::Parameter<DataType> coarse_tol_abs = parameter(DataType(1E+10)).property("coarse-gmres-parameters/tol-abs").name("Coarse Absolute Tolerance");
      FEAT::Parameter<DataType> coarse_inner_rescale = parameter(DataType(1E+1)).property("coarse-gmres-parameters/inner-rescale").name("Coarse Inner Rescale");
      FEAT::Parameter<FEAT::String> coarse_orthogonalizer = parameter(FEAT::String("cgs2")).property("coarse-gmres-parameters/orthogonalizer").name("Coarse GMRES Orthogonalizer").validator(is_orthogonalizer);

      FEAT::Parameter<FEAT::Index> min_nl_iter = parameter(FEAT::Index(1)).property("min-nl-iter").name("Min Nonlinear Iter");
      FEAT::Parameter<FEAT::Index> max_nl_iter = parameter(FEAT::Index(50)).property("max-nl-iter").name("Max Nonlinear Iter");
      FEAT::Parameter<DataType> nl_tol_abs = parameter(DataType(1E-3)).property("nl-tol-abs").name("Nonlinear Absolute Tolerance");
      FEAT::Parameter<DataType> nl_tol_rel = parameter(DataType(1E-5)).property("nl-tol-rel").name("Nonlinear Relative Tolerance");
      FEAT::Parameter<DataType> nl_tol_low_abs = parameter(DataType(1E-8)).property("nl-tol-low-abs").name("Nonlinear Low Absolute Tolerance");
      FEAT::Parameter<DataType> nl_stagrate = parameter(DataType(0.98)).property("nl-stagrate").name("Nonlinear Stagnation Rate");
      FEAT::Parameter<DataType> fixed_ls_tol = parameter(DataType(-1)).property("fixed-ls-tol").name("Fixed Linear Solver Tolerance");
      FEAT::Parameter<FEAT::Index> max_stagnations = parameter(FEAT::Index(4)).property("max-stagnations").name("Max Stagnations");
      FEAT::Parameter<bool> init_with_stokes = parameter(false).property("init-with-stokes").name("Initialize With Stokes").parser(parse_config_bool);
      FEAT::Parameter<bool> reinit_fbm_interface_filter = parameter(false).property("reinit-fbm-interface-filter").name("Reinit FBM Interface Filter").parser(parse_config_bool);
      FEAT::Parameter<FEAT::Index> min_picard_steps = parameter(FEAT::Index(3)).property("min-picard-steps").name("Min Picard Steps");

      FEAT::Parameter<bool> adaptive_step_size = parameter(true).property("adaptive-step-size").name("Adaptive Time Step Size").parser(parse_config_bool);
      FEAT::Parameter<DataType> pseudo_ts_size = parameter(DataType(1E-3)).property("pseudo-ts-size").name("Pseudo Time Step Size");
      FEAT::Parameter<DataType> adapt_factor = parameter(DataType(4)).property("adapt-factor").name("Adaptivity Factor");
      FEAT::Parameter<DataType> stability_tol = parameter(DataType(0.05)).property("stability-tol").name("Stability Tolerance");
      FEAT::Parameter<FEAT::Index> stability_sample_size = parameter(FEAT::Index(8)).property("stability-sample-size").name("Stability Sample Size");
      FEAT::Parameter<DataType> stagnation_increase = parameter(DataType(100)).property("stagnation-increase").name("Stagnation Increase");

    #ifdef FEAT_HAVE_TRILINOS
      FEAT::Parameter<FEAT::Index> frosch_gmres_dim = parameter(FEAT::Index(32)).property("frosch-gmres-parameters/gmres-dim").name("FROSch GMRES Dimension");
      FEAT::Parameter<FEAT::Index> frosch_min_iter = parameter(FEAT::Index(3)).property("frosch-gmres-parameters/min-iter").name("FROSch Min Iter");
      FEAT::Parameter<FEAT::Index> frosch_max_iter = parameter(FEAT::Index(128)).property("frosch-gmres-parameters/max-iter").name("FROSch Max Iter");
      FEAT::Parameter<DataType> frosch_tol_rel = parameter(DataType(1E-4)).property("frosch-gmres-parameters/tol-rel").name("FROSch Relative Tolerance");
      FEAT::Parameter<DataType> frosch_tol_abs = parameter(DataType(1E+10)).property("frosch-gmres-parameters/tol-abs").name("FROSch Absolute Tolerance");
      FEAT::Parameter<DataType> frosch_inner_rescale = parameter(DataType(1E+1)).property("frosch-gmres-parameters/inner-rescale").name("FROSch Inner Rescale");

      FEAT::Parameter<FEAT::String> frosch_xml = parameter(FEAT::String("")).property("frosch-parameters/frosch-xml").name("FROSch XML");
      FEAT::Parameter<std::deque<int>> frosch_subregions = parameter(std::deque<int>{1}).property("frosch-parameters/frosch-subregions").name("FROSch Subregions");
      FEAT::Parameter<int> frosch_overlap = parameter(1).property("frosch-parameters/frosch-overlap").name("FROSch Overlap");
      FEAT::Parameter<FEAT::String> frosch_direct_solver_coarse = parameter(FEAT::String("umfpack")).property("frosch-parameters/frosch-direct-solver-coarse").name("FROSch Direct Solver Coarse").validator(is_direct_solver);
      FEAT::Parameter<FEAT::String> frosch_direct_solver_overlapping = parameter(FEAT::String("umfpack")).property("frosch-parameters/frosch-direct-solver-overlapping").name("FROSch Direct Solver Overlapping").validator(is_direct_solver);
      FEAT::Parameter<FEAT::String> frosch_direct_solver_extension = parameter(FEAT::String("umfpack")).property("frosch-parameters/frosch-direct-solver-extension").name("FROSch Direct Solver Extension").validator(is_direct_solver);
      FEAT::Parameter<FEAT::String> frosch_precond_type = parameter(FEAT::String("twolevelblock")).property("frosch-parameters/frosch-precond-type").name("FROSch Preconditioner Type").validator(is_frosch_precond);
      FEAT::Parameter<FEAT::String> frosch_partition_type = parameter(FEAT::String("parmetis")).property("frosch-parameters/frosch-partition-type").name("FROSch Partition Type").validator(is_frosch_partition_type);
      FEAT::Parameter<FEAT::String> frosch_partition_approach = parameter(FEAT::String("repartition")).property("frosch-parameters/frosch-partition-approach").name("FROSch Partition Approach").validator(is_frosch_partition_approach);
      FEAT::Parameter<FEAT::String> frosch_combine_overlap = parameter(FEAT::String("restricted")).property("frosch-parameters/frosch-combine-overlap").name("FROSch Combine Overlap").validator(is_frosch_combine_overlap);
      FEAT::Parameter<FEAT::String> frosch_combine_level = parameter(FEAT::String("additive")).property("frosch-parameters/frosch-combine-level").name("FROSch Combine Level").validator(is_frosch_combine_level);
      FEAT::Parameter<FEAT::String> frosch_ipou_velocity = parameter(FEAT::String("gdswstar")).property("frosch-parameters/frosch-ipou-velocity").name("FROSch IPOU Velocity").validator(is_frosch_ipou);
      FEAT::Parameter<FEAT::String> frosch_ipou_pressure = parameter(FEAT::String("rgdsw")).property("frosch-parameters/frosch-ipou-pressure").name("FROSch IPOU Pressure").validator(is_frosch_ipou);
      FEAT::Parameter<bool> frosch_verbose = parameter(false).property("frosch-parameters/frosch-verbose").name("FROSch Verbose").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_exclude_velocity_pressure = parameter(false).property("frosch-parameters/frosch-exclude-velocity-pressure").name("FROSch Exclude Velocity Pressure").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_exclude_pressure_velocity = parameter(false).property("frosch-parameters/frosch-exclude-pressure-velocity").name("FROSch Exclude Pressure Velocity").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_reuse_sf_ao = parameter(false).property("frosch-parameters/frosch-reuse-sf-algebraic-overlapping-solver").name("FROSch Reuse SF AO").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_reuse_sf_cm = parameter(false).property("frosch-parameters/frosch-reuse-sf-coarse-matrix").name("FROSch Reuse SF CM").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_reuse_sf_ext = parameter(false).property("frosch-parameters/frosch-reuse-sf-extension-solver").name("FROSch Reuse SF Extension").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_reuse_cm = parameter(false).property("frosch-parameters/frosch-reuse-coarse-matrix").name("FROSch Reuse Coarse Matrix").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_reuse_cb = parameter(false).property("frosch-parameters/frosch-reuse-coarse-basis").name("FROSch Reuse Coarse Basis").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_use_timers = parameter(false).property("frosch-parameters/frosch-use-timers").name("FROSch Use Timers").parser(parse_config_bool);
      FEAT::Parameter<bool> frosch_print_internal = parameter(false).property("frosch-parameters/frosch-print-internal").name("FROSch Print Internal").parser(parse_config_bool);
    #endif
    };

    template<typename DomainLevel_>
    inline static FEAT::Adjacency::Coloring get_element_coloring(const DomainLevel_& dom)
    {
      // create coloring
      const auto& verts_at_elem = dom.get_mesh().template get_index_set<DomainLevel_::MeshType::shape_dim, 0>();
      FEAT::Adjacency::Graph elems_at_vert(FEAT::Adjacency::RenderType::transpose, verts_at_elem);
      FEAT::Adjacency::Graph elems_at_elem(FEAT::Adjacency::RenderType::injectify, verts_at_elem, elems_at_vert);
      return FEAT::Adjacency::Coloring(elems_at_elem);
    }
  }
  // CRTP interface for Flow System Assembler
  template<typename Derived_, typename SystemLevel_>
  class SteadyStokesFlowSolverBaseCRTP
  {
  public:
    typedef SystemLevel_ SystemLevelType;
    typedef std::deque<std::shared_ptr<SystemLevelType>> SystemLevels;
    typedef typename SystemLevelType::DataType DataType;

    typedef typename SystemLevelType::GlobalSystemVector DefectVectorType;
    typedef typename SystemLevelType::GlobalSystemVector VectorType;
    typedef std::shared_ptr<FEAT::Solver::SolverBase<DefectVectorType>> BaseSolver;
    typedef std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>> BaseIterativeSolver;

    SystemLevels system_levels;
    BaseSolver base_solver;
    BaseIterativeSolver iter_solver;
    mutable FEAT::StopWatch watch_total, watch_apply, watch_init_symbolic, watch_init_numeric;

    FEAT::PreferredBackend solver_backend = FEAT::PreferredBackend::generic;

    static constexpr int padlen = 30;
    static constexpr char pc = '.';

  protected:
    /**
      * \brief Casts \c this to its true type
      *
      * \returns A Derived_ reference to \c this
      */
    Derived_& cast() {return static_cast<Derived_&>(*this);}

    /// \copydoc cast()
    const Derived_& cast() const {return static_cast<const Derived_&>(*this);}

    bool _parse(Intern::SolverParser<DataType>& DOXY(parser))
    {
      return true;
    }

    void _set_tol_abs(DataType tol)
    {
      if(iter_solver) iter_solver->set_tol_abs(tol);
    }

    DataType _get_tol_abs() const
    {
      return iter_solver ? iter_solver->get_tol_abs() : DataType(0);
    }

    void _set_tol_rel(DataType tol)
    {
      if(iter_solver) iter_solver->set_tol_rel(tol);
    }

    DataType _get_tol_rel() const
    {
      return iter_solver ? iter_solver->get_tol_rel() : DataType(0);
    }

    void _set_max_iter(FEAT::Index max_iter)
    {
      if(iter_solver) iter_solver->set_max_iter(max_iter);
    }

    FEAT::Index _get_max_iter() const
    {
      return iter_solver ? iter_solver->get_max_iter() : FEAT::Index(0);
    }

    FEAT::Index _get_num_iter() const
    {
      return iter_solver ? iter_solver->get_num_iter() :0u;
    }

    DataType _get_def_initial() const
    {
      return iter_solver ? iter_solver->get_def_initial() : DataType(0);
    }

    DataType _get_def_final() const
    {
      return iter_solver ? iter_solver->get_def_final() : DataType(0);
    }

    FEAT::Solver::Status _apply(DefectVectorType& vec_cor, const DefectVectorType& vec_def) const
    {
      return base_solver->apply(vec_cor, vec_def);
    }

    FEAT::Solver::Status _correct(DefectVectorType& vec_cor, const DefectVectorType& vec_def) const
    {
      XASSERTM(this->iter_solver, "Correct without iterative solver initialized");
      return iter_solver->correct(vec_cor, vec_def);
    }

    FEAT::Solver::Status _solve(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs) const
    {
      if(iter_solver)
        return iter_solver->correct(vec_sol, vec_rhs);
      else
        return FEAT::Solver::solve(*base_solver, vec_sol, vec_rhs, this->system_levels.front()->matrix_sys, this->system_levels.front()->filter_sys);
    }

    DefectVectorType _create_vector() const
    {
      XABORTM("Create vector has to be specialized by child class");
      return DefectVectorType();
    }

    void _init_symbolic()
    {
      base_solver->init_symbolic();
    }

    void _done_symbolic()
    {
      base_solver->done_symbolic();
    }

    void _init_numeric()
    {
      base_solver->init_numeric();
    }

    void _done_numeric()
    {
      base_solver->done_numeric();
    }

    void _set_backend()
    {
      FEAT::Backend::set_preferred_backend(this->solver_backend);
    }

    FEAT::String _format_multigrid_timings() const
    {
      return FEAT::String{};
    }

    FEAT::String _format_parameters() const
    {
      return FEAT::String{};
    }

  public:

    SteadyStokesFlowSolverBaseCRTP() = default;
    ~SteadyStokesFlowSolverBaseCRTP() = default;
    SteadyStokesFlowSolverBaseCRTP(const SteadyStokesFlowSolverBaseCRTP&) = default;
    SteadyStokesFlowSolverBaseCRTP(SteadyStokesFlowSolverBaseCRTP&&) = default;
    SteadyStokesFlowSolverBaseCRTP& operator=(const SteadyStokesFlowSolverBaseCRTP&) = default;
    SteadyStokesFlowSolverBaseCRTP& operator=(SteadyStokesFlowSolverBaseCRTP&&) = default;

    SteadyStokesFlowSolverBaseCRTP(const SystemLevels& sys_levels)
      : system_levels(sys_levels)
    {}

    SteadyStokesFlowSolverBaseCRTP(SystemLevels&& sys_levels)
      : system_levels(std::move(sys_levels))
    {}

    void init_symbolic()
    {
      watch_total.start();
      watch_init_symbolic.start();
      this->cast()._init_symbolic();
      watch_init_symbolic.stop();
      watch_total.stop();
    }

    void init_numeric()
    {
      watch_total.start();
      watch_init_numeric.start();
      this->cast()._init_numeric();
      watch_init_numeric.stop();
      watch_total.stop();
    }

    void done_symbolic()
    {
      watch_total.start();
      this->cast()._done_symbolic();
      watch_total.stop();
    }

    void done_numeric()
    {
      watch_total.start();
      this->cast()._done_numeric();
      watch_total.stop();
    }

    SystemLevels& get_systems()
    {
      return this->system_levels;
    }

    const SystemLevels& get_systems() const
    {
      return this->system_levels;
    }

    void set_tol_abs(DataType tol)
    {
      this->cast()._set_tol_abs(tol);
    }

    DataType get_tol_abs() const
    {
      return this->cast()._get_tol_abs();
    }

    void set_tol_rel(DataType tol)
    {
      this->cast()._set_tol_rel(tol);
    }

    DataType get_tol_rel() const
    {
      return this->cast()._get_tol_rel();
    }

    void set_max_iter(FEAT::Index max_iter)
    {
      this->cast()._set_max_iter(max_iter);
    }

    FEAT::Index get_max_iter() const
    {
      return this->cast()._get_max_iter();
    }

    FEAT::Index get_num_iter() const
    {
      return this->cast()._get_num_iter();
    }

    DataType get_def_final() const
    {
      return this->cast()._get_def_final();
    }

    DataType get_def_initial() const
    {
      return this->cast()._get_def_initial();
    }

    void set_plot_mode(FEAT::Solver::PlotMode plotmode)
    {
      if(this->iter_solver) iter_solver->set_plot_mode(plotmode);
    }

    FEAT::Solver::Status apply(DefectVectorType& vec_cor, const DefectVectorType& vec_def) const
    {
      watch_total.start();
      watch_apply.start();
      auto status = this->cast()._apply(vec_cor, vec_def);
      watch_apply.stop();
      watch_total.stop();
      return status;
    }

    FEAT::Solver::Status correct(DefectVectorType& vec_cor, const DefectVectorType& vec_def) const
    {
      watch_total.start();
      watch_apply.start();
      auto status = this->cast()._correct(vec_cor, vec_def);
      watch_apply.stop();
      watch_total.stop();
      return status;
    }

    FEAT::Solver::Status solve(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs) const
    {
      watch_total.start();
      watch_apply.start();
      auto status = this->cast()._solve(vec_sol, vec_rhs);
      watch_apply.stop();
      watch_total.stop();
      return status;
    }

    DefectVectorType create_defect() const
    {
      watch_total.start();
      DefectVectorType tmp = this->cast()._create_vector();
      watch_total.stop();
      return tmp;
    }

    DefectVectorType create_correction() const
    {
      watch_total.start();
      DefectVectorType tmp = this->cast()._create_vector();
      watch_total.stop();
      return tmp;
    }

    void set_backend()
    {
      return this->cast()._set_backend();
    }

    FEAT::String format_timings() const
    {
      return this->cast()._format_timings();
    }

    FEAT::String format_multigrid_timings() const
    {
      return this->cast()._format_multigrid_timings();
    }

    FEAT::String format_string() const
    {
      return this->cast()._format_parameters();
    }

    double get_total_time_elapsed() const
    {
      return this->watch_total.elapsed();
    }

  }; // class SteadyFlowSolverBaseCRTP<...>

  template<typename SystemMatrix_, typename SystemTransfer_, typename SystemFilter_>
  class MultigridVankaComponent
  {
  public:

    typedef SystemMatrix_ SystemMatrix;
    typedef typename SystemMatrix_::LocalMatrixType LocalSystemMatrix;
    typedef SystemTransfer_ SystemTransfer;
    typedef SystemFilter_ SystemFilter;
    typedef typename SystemFilter_::LocalFilterType LocalSystemFilter;
    typedef typename SystemMatrix_::VectorTypeL DefectVectorType;

    typedef typename LocalSystemMatrix::DataType DataType;
    typedef typename LocalSystemMatrix::IndexType IndexType;

    // static constexpr int dim = SystemLevelType::dim;
    typedef std::shared_ptr<FEAT::Solver::SolverBase<DefectVectorType>> BaseSolver;
    typedef std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>> BaseIterativeSolver;

    BaseSolver base_solver;
    BaseIterativeSolver iter_solver;

    std::shared_ptr<FEAT::Solver::MultiGridHierarchy<SystemMatrix, SystemFilter, SystemTransfer>> multigrid_hierarchy;
    // hold a deque of the local vanka solvers
    std::deque<std::shared_ptr<FEAT::Solver::AmaVanka<LocalSystemMatrix, LocalSystemFilter>>> vanka_solver;
    std::deque<std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>>> smoother;
    std::shared_ptr<FEAT::Solver::SolverBase<DefectVectorType>> coarse_solver;
    std::shared_ptr<FEAT::Solver::MultiGrid<SystemMatrix, SystemFilter, SystemTransfer>> multigrid;
    #ifdef FEAT_HAVE_TRILINOS
    std::shared_ptr<FEAT::Solver::StokesFROSchPreconditioner<SystemMatrix, SystemFilter>> frosch_precond;
    std::shared_ptr<FEAT::Solver::Trilinos::FROSchParameterList> frosch_params;
    std::vector<int> _frosch_subregions;
    int _frosch_overlap;
    FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver _frosch_directsolver_coarse;
    FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver _frosch_directsolver_overlapping;
    FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver _frosch_directsolver_extension;
    FEAT::Solver::Trilinos::FROSchParameterList::CombineLevel _frosch_combine_lvl;
    FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner _frosch_precond_type;
    FEAT::Solver::Trilinos::FROSchParameterList::PartitionType _frosch_parti;
    FEAT::Solver::Trilinos::FROSchParameterList::PartitionApproach _frosch_approach;
    FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap _frosch_combine_overlap;
    FEAT::Solver::Trilinos::FROSchParameterList::IPOU _frosch_ipou_velo;
    FEAT::Solver::Trilinos::FROSchParameterList::IPOU _frosch_ipou_pres;
    bool _frosch_verbose;
    bool _frosch_exclude_velocity_pressure;
    bool _frosch_exclude_pressure_velocity;
    // reuse symbolic factorization of overlapping operator
    bool _frosch_reuse_sf_ao;
    // reuse symbolic factorization of coarse matrix
    bool _frosch_reuse_sf_cm;
    // reuse symbolic factorization of extension solver
    bool _frosch_reuse_sf_ext;
    // reuse coarse matrix
    bool _frosch_reuse_cm;
    // reuse coarse basis
    bool _frosch_reuse_cb;
    FEAT::String _frosch_xml;
    bool _frosch_use_timers;
    bool _frosch_print_internal;
    #endif
    struct CoarseParams
    {
      FEAT::Index _gmres_dim = FEAT::Index(200);
      FEAT::Index _miniter = FEAT::Index(1);
      FEAT::Index _maxiter = FEAT::Index(500);
      DataType _tol_rel = DataType(1E-2);
      DataType _tol_abs = DataType(1E+10);
      DataType _inner_rescale = DataType(1E+1);
      FEAT::Solver::KrylovOrthogonalizerType _ortho = FEAT::Solver::KrylovOrthogonalizerType::cgs2;
    };
    CoarseParams coarse_gmres_params;
  #ifdef FEAT_HAVE_TRILINOS
    CoarseParams frosch_gmres_params;
  #endif

    std::size_t domain_virtual_size;
    DataType smoother_damp;
    DataType mg_tol_rel;
    FEAT::Index smooth_gmres_dim;
    FEAT::Index solve_gmres_dim;
    FEAT::Index smooth_steps;
    FEAT::Index min_iter;
    FEAT::Index max_iter;
    FEAT::Index min_stag_iter;
    FEAT::PreferredBackend solver_backend;
    FEAT::Solver::MultiGridCycle cycle;
    mutable FEAT::StopWatch watch_mg_hirarch;
    bool coarse_frosch;
    bool coarse_gmres;
    bool voxel_vanka;
    bool coarse_solver_info;

  public:

    bool parse(Intern::SolverParser<DataType>& parser)
    {
      min_iter = *parser.min_mg_iter;
      max_iter = *parser.max_mg_iter;
      smooth_gmres_dim = *parser.smooth_gmres;
      solve_gmres_dim = *parser.solve_gmres;
      smooth_steps = *parser.smooth_steps;
      smoother_damp = *parser.smooth_damp;
      mg_tol_rel = *parser.mg_tol_rel;
      solver_backend = Intern::backend_from_string(*parser.backend);
      coarse_gmres |= *parser.coarse_gmres;
      coarse_frosch |= *parser.coarse_frosch;
      coarse_solver_info |= *parser.coarse_info;

      this->_parse_coarse_params(coarse_gmres_params, parser);

    #ifdef FEAT_HAVE_TRILINOS
      this->_parse_frosch_coarse_params(frosch_gmres_params, parser);
      this->_parse_frosch_params(parser);
    #else
      coarse_frosch = false;
    #endif // FEAT_HAVE_TRILINOS
      return true;
    }

    bool parse(const FEAT::PropertyMap* prop)
    {
      if(prop == nullptr)
        return true;
      Intern::SolverParser<DataType> parser;
      if(!parser.parse(*prop))
        return false;
      return this->parse(parser);
    }

    #ifdef FEAT_HAVE_TRILINOS
    bool _parse_frosch_params(const Intern::SolverParser<DataType>& parser)
    {
      using FROSchList = FEAT::Solver::Trilinos::FROSchParameterList;

      auto direct_solver = [](const FEAT::String& option)
      {
        if(option.compare_no_case("umfpack") == 0) return FROSchList::DirectSolver::UMFPACK;
        if(option.compare_no_case("mumps") == 0) return FROSchList::DirectSolver::MUMPS;
        if(option.compare_no_case("ilu") == 0) return FROSchList::DirectSolver::ILU;
        if(option.compare_no_case("klu") == 0) return FROSchList::DirectSolver::KLU;
        XABORTM("Unknown direct solver " + option);
        return FROSchList::DirectSolver::UMFPACK;
      };

      _frosch_xml = *parser.frosch_xml;
      _frosch_subregions.assign(parser.frosch_subregions->begin(), parser.frosch_subregions->end());
      XASSERTM(!_frosch_subregions.empty() && _frosch_subregions.back() == 1, "Last number of subregions needs to be == 1");
      _frosch_overlap = *parser.frosch_overlap;

      _frosch_directsolver_coarse = direct_solver(*parser.frosch_direct_solver_coarse);
      _frosch_directsolver_overlapping = direct_solver(*parser.frosch_direct_solver_overlapping);
      _frosch_directsolver_extension = direct_solver(*parser.frosch_direct_solver_extension);

      if(parser.frosch_precond_type->compare_no_case("onelevel") == 0)
        _frosch_precond_type = FROSchList::Preconditioner::ONELEVEL;
      else if(parser.frosch_precond_type->compare_no_case("twolevel") == 0)
        _frosch_precond_type = FROSchList::Preconditioner::TWOLEVEL;
      else if(parser.frosch_precond_type->compare_no_case("twolevelblock") == 0)
        _frosch_precond_type = FROSchList::Preconditioner::TWOLEVELBLOCK;

      if(parser.frosch_partition_type->compare_no_case("PHG") == 0)
        _frosch_parti = FROSchList::PartitionType::PHG;
      else if(parser.frosch_partition_type->compare_no_case("PARMETIS") == 0)
        _frosch_parti = FROSchList::PartitionType::PARMETIS;
      else if(parser.frosch_partition_type->compare_no_case("SCOTCH") == 0)
        _frosch_parti = FROSchList::PartitionType::SCOTCH;
      else if(parser.frosch_partition_type->compare_no_case("PTSCOTCH") == 0)
        _frosch_parti = FROSchList::PartitionType::PTSCOTCH;
      else if(parser.frosch_partition_type->compare_no_case("BLOCK") == 0)
        _frosch_parti = FROSchList::PartitionType::BLOCK;

      _frosch_approach = (parser.frosch_partition_approach->compare_no_case("partition") == 0) ? FROSchList::PartitionApproach::PARTITION : FROSchList::PartitionApproach::REPARTITION;

      if(parser.frosch_combine_overlap->compare_no_case("restricted") == 0)
        _frosch_combine_overlap = FROSchList::CombineOverlap::RESTRICTED;
      else if(parser.frosch_combine_overlap->compare_no_case("full") == 0)
        _frosch_combine_overlap = FROSchList::CombineOverlap::FULL;
      else if(parser.frosch_combine_overlap->compare_no_case("averaging") == 0)
        _frosch_combine_overlap = FROSchList::CombineOverlap::AVERAGING;

      _frosch_combine_lvl = (parser.frosch_combine_level->compare_no_case("additive") == 0) ? FROSchList::CombineLevel::ADDITIVE : FROSchList::CombineLevel::MULTIPLICATIVE;

      auto ipou = [](const FEAT::String& option)
      {
        if(option.compare_no_case("gdsw") == 0) return FROSchList::IPOU::GDSW;
        if(option.compare_no_case("gdswstar") == 0) return FROSchList::IPOU::GDSWSTAR;
        if(option.compare_no_case("rgdsw") == 0) return FROSchList::IPOU::RGDSW;
        XABORTM("Unknown IPOU " + option);
        return FROSchList::IPOU::GDSW;
      };
      _frosch_ipou_velo = ipou(*parser.frosch_ipou_velocity);
      _frosch_ipou_pres = ipou(*parser.frosch_ipou_pressure);

      _frosch_verbose = *parser.frosch_verbose;
      _frosch_exclude_velocity_pressure = *parser.frosch_exclude_velocity_pressure;
      _frosch_exclude_pressure_velocity = *parser.frosch_exclude_pressure_velocity;
      _frosch_reuse_sf_ao = *parser.frosch_reuse_sf_ao;
      _frosch_reuse_sf_cm = *parser.frosch_reuse_sf_cm;
      _frosch_reuse_sf_ext = *parser.frosch_reuse_sf_ext;
      _frosch_reuse_cm = *parser.frosch_reuse_cm;
      _frosch_reuse_cb = *parser.frosch_reuse_cb;
      _frosch_use_timers = *parser.frosch_use_timers;
      _frosch_print_internal = *parser.frosch_print_internal;

      return true;
    }

    bool _setup_frosch_params()
    {
      // for now local variables, will be used in class as configuration needs arise
      // for now, always use _frosch_nlvels = 1, else we would have to use vectors for the interface
      // if we have an xml file. we always only parse this
      if(!_frosch_xml.empty())
      {
        //if(_logger) _logger->print("Use frosch xml file");
        this->frosch_params->read_from_xml_file(_frosch_xml);
        return true;
      }

      this->frosch_params->set_nlevels(int(_frosch_subregions.size()));
      this->frosch_params->set_subregions(_frosch_subregions);
      this->frosch_params->set_precond(_frosch_precond_type);
      this->frosch_params->set_coarse_solver(_frosch_directsolver_coarse);
      this->frosch_params->set_solvers(_frosch_directsolver_overlapping, _frosch_directsolver_extension);
      this->frosch_params->set_overlaps(_frosch_overlap);
      this->frosch_params->set_combine_overlap(_frosch_combine_overlap);
      this->frosch_params->set_combine_lvl(_frosch_combine_lvl);
      this->frosch_params->set_parti_types(_frosch_parti);
      this->frosch_params->set_parti_approach(_frosch_approach);
      this->frosch_params->set_ipous(_frosch_ipou_velo, _frosch_ipou_pres);
      this->frosch_params->set_excludes_velo_pres(_frosch_exclude_velocity_pressure);
      this->frosch_params->set_excludes_pres_velo(_frosch_exclude_pressure_velocity);
      this->frosch_params->set_reuse_sf(_frosch_reuse_sf_ao, _frosch_reuse_sf_cm);
      this->frosch_params->set_reuse_coarse(_frosch_reuse_cm, _frosch_reuse_cb);

      this->frosch_params->set_print(_frosch_print_internal);
      this->frosch_params->set_use_timer(_frosch_use_timers);
      this->frosch_params->set_verbose(_frosch_verbose);

      return true;
    }
    #endif

    bool _parse_coarse_params(CoarseParams& params, const Intern::SolverParser<DataType>& parser)
    {
      params._gmres_dim = *parser.coarse_gmres_dim;
      params._miniter = *parser.coarse_min_iter;
      params._maxiter = *parser.coarse_max_iter;
      params._tol_rel = *parser.coarse_tol_rel;
      params._tol_abs = *parser.coarse_tol_abs;
      params._inner_rescale = *parser.coarse_inner_rescale;
      params._ortho = Intern::orthogonalizer_from_string(*parser.coarse_orthogonalizer);
      return true;
    }

  #ifdef FEAT_HAVE_TRILINOS
    bool _parse_frosch_coarse_params(CoarseParams& params, const Intern::SolverParser<DataType>& parser)
    {
      params._gmres_dim = *parser.frosch_gmres_dim;
      params._miniter = *parser.frosch_min_iter;
      params._maxiter = *parser.frosch_max_iter;
      params._tol_rel = *parser.frosch_tol_rel;
      params._tol_abs = *parser.frosch_tol_abs;
      params._inner_rescale = *parser.frosch_inner_rescale;
      return true;
    }
  #endif

    void clear_solver()
    {
      #ifdef FEAT_HAVE_TRILINOS
      this->frosch_precond.reset();
      this->frosch_params.reset();
      #endif
      this->base_solver.reset();
      this->iter_solver.reset();
      this->multigrid.reset();
      this->multigrid_hierarchy.reset();
      this->vanka_solver.clear();
      this->coarse_solver.reset();
      this->smoother.clear();
    }

    template<typename Domain_, typename System_, typename... TArgs, template<typename...> typename Container_>
    void set_solver(const Container_<std::shared_ptr<System_>, TArgs...>& system_levels, const Domain_& domain)
    {
      // typedef System_ SystemLevelType;
      // XASSERTM(this->_init_params, "You have to initialize the system before calling this function");
      XASSERTM(system_levels.size() > 0, "System levels not set");
      this->clear_solver();
      const auto& matrix_sys = system_levels.front()->matrix_sys;
      const auto& filter_sys = system_levels.front()->filter_sys;

      domain_virtual_size = domain.size_virtual();

      // if we use frosch, initialize frosch params now
      #ifdef FEAT_HAVE_TRILINOS
      if(coarse_frosch && (system_levels.size() == domain_virtual_size))
      {
        // create on same communicater our coarse level uses
        frosch_params = std::make_shared<FEAT::Solver::Trilinos::FROSchParameterList>(domain.back().layer().comm(), 3, FEAT::Solver::Trilinos::FROSchParameterList::SADDLEPOINT);
        // if(_logger) _logger->print("Parsing from " + frosch_xml, info);

        // frosch_params->read_from_xml_file(frosch_xml);
        this->_setup_frosch_params();
        frosch_params->create_core();
      }

      #endif

      // case only one level
      if(domain_virtual_size == std::size_t(1))
      {
        // create direct solver or FMGRES?
        if(!coarse_frosch && !coarse_gmres && FEAT::Solver::have_direct_sparse_solver(matrix_sys, filter_sys))
        {
          auto backendt = FEAT::Backend::get_preferred_backend();
          FEAT::Backend::set_preferred_backend(solver_backend);
          this->base_solver = FEAT::Solver::new_direct_sparse_solver(matrix_sys, filter_sys);
          FEAT::Backend::set_preferred_backend(backendt);
        }
#ifdef FEAT_HAVE_TRILINOS
        else if(coarse_frosch)
        {
          this->frosch_precond = FEAT::Solver::new_stokes_frosch(matrix_sys, filter_sys, *frosch_params);

          std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>> solver_iterative;

          if(frosch_gmres_params._gmres_dim > 0)
          {
            solver_iterative = Solver::new_fgmres(matrix_sys, filter_sys, frosch_gmres_params._gmres_dim, frosch_gmres_params._inner_rescale, this->frosch_precond);
          }
          else
          {
            solver_iterative = Solver::new_richardson(matrix_sys, filter_sys, DataType(1), this->frosch_precond);
          }
          solver_iterative->set_tol_rel(mg_tol_rel);
          // solver_iterative->set_tol_abs(tol_abs);
          solver_iterative->set_tol_abs_low(DataType(1e-14));
          solver_iterative->set_min_iter(min_iter);
          solver_iterative->set_max_iter(max_iter);
          // solver_iterative->set_plot_mode(FEAT::Solver::PlotMode::summary);
          solver_iterative->set_plot_mode(coarse_solver_info ? FEAT::Solver::PlotMode::summary : FEAT::Solver::PlotMode::none);

          this->base_solver = solver_iterative;
        }
#endif // FEAT_HAVE_TRILINOS
        else
        {
          if(voxel_vanka)
          {
            if constexpr(Intern::supports_element_coloring<Domain_>)
              this->vanka_solver.push_back(FEAT::Solver::new_voxel_amavanka(system_levels.front()->local_matrix_sys, filter_sys.local(), domain.front()->element_coloring));
            else
            {
              this->vanka_solver.push_back(FEAT::Solver::new_voxel_amavanka(system_levels.front()->local_matrix_sys, filter_sys.local(), Intern::get_element_coloring(*domain.front())));
            }
          }
          else
            this->vanka_solver.push_back(FEAT::Solver::new_amavanka(system_levels.front()->local_matrix_sys, filter_sys.local()));
          this->vanka_solver.front()->set_skip_singular(true);
          auto schwarz = FEAT::Solver::new_schwarz_precond(this->vanka_solver.front(), filter_sys);
          if(solve_gmres_dim > FEAT::Index(0)) // todo: is this smart?
            this->base_solver = this->iter_solver = FEAT::Solver::new_gmres(matrix_sys, filter_sys, solve_gmres_dim, coarse_gmres_params._inner_rescale, schwarz);
          else
            this->base_solver = this->iter_solver = FEAT::Solver::new_gmres(matrix_sys, filter_sys, 16, coarse_gmres_params._inner_rescale, schwarz);

          // configure solver
          this->iter_solver->set_plot_name("FGMRES-AmaVanka");
          this->iter_solver->set_min_iter(min_iter);
          this->iter_solver->set_max_iter(max_iter);
          this->iter_solver->set_tol_rel(mg_tol_rel);
          this->iter_solver->set_min_stag_iter(min_stag_iter);
          this->iter_solver->set_plot_mode(coarse_solver_info ? FEAT::Solver::PlotMode::summary : FEAT::Solver::PlotMode::none);
        }
        return;
      }

      multigrid_hierarchy = std::make_shared<FEAT::Solver::MultiGridHierarchy<
                              SystemMatrix, SystemFilter, SystemTransfer>>(domain_virtual_size);

      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        auto& lvl = *system_levels.at(i);

        if((i+1) < domain_virtual_size)
        {
          // create Schwarz-AmaVanka
          std::shared_ptr<FEAT::Solver::AmaVanka<LocalSystemMatrix, LocalSystemFilter>> vanka(nullptr);
          if(voxel_vanka)
          {
            if constexpr(Intern::supports_element_coloring<Domain_>)
              vanka = FEAT::Solver::new_voxel_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local(), domain.at(i)->element_coloring);
            else
              vanka = FEAT::Solver::new_voxel_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local(), Intern::get_element_coloring(*domain.at(i)));
          }
          else
            vanka = FEAT::Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
          vanka->set_skip_singular(true);
          // ama_vankas.push_back(vanka);
          auto schwarz = FEAT::Solver::new_schwarz_precond(vanka, lvl.filter_sys);
          schwarz->set_ignore_status(true);

          // create smoother: either GMRES or Richardson
          std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>> smoother_l;
          if(smooth_gmres_dim > FEAT::Index(0))
            smoother_l = FEAT::Solver::new_gmres(lvl.matrix_sys, lvl.filter_sys, smooth_gmres_dim, 0.0, schwarz);
          else
            smoother_l = FEAT::Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smoother_damp, schwarz);
          smoother_l->set_min_iter(smooth_steps);
          smoother_l->set_max_iter(smooth_steps);
          // smoother_l->set_plot_mode(FEAT::Solver::PlotMode::all);
          // save our smoother and vanka
          this->multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother_l, smoother_l, smoother_l);
          this->smoother.push_back(std::move(smoother_l));
          this->vanka_solver.push_back(std::move(vanka));
        }
        else if(!coarse_frosch && !coarse_gmres && FEAT::Solver::have_direct_sparse_solver(lvl.matrix_sys, lvl.filter_sys))
        {
          auto backendt = FEAT::Backend::get_preferred_backend();
          FEAT::Backend::set_preferred_backend(solver_backend);
          // create UMFPACK coarse grid solver
          this->coarse_solver = FEAT::Solver::new_direct_sparse_solver(lvl.matrix_sys, lvl.filter_sys);
          this->multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, this->coarse_solver);
          FEAT::Backend::set_preferred_backend(backendt);
        }
#ifdef FEAT_HAVE_TRILINOS
        else if(coarse_frosch)
        {
          this->frosch_precond = Solver::new_stokes_frosch(lvl.matrix_sys, lvl.filter_sys, *frosch_params);

          std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>> solver_iterative;

          if(frosch_gmres_params._gmres_dim > 0)
          {
            solver_iterative = Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, frosch_gmres_params._gmres_dim, frosch_gmres_params._inner_rescale, this->frosch_precond);
          }
          else
          {
            solver_iterative = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, DataType(1), this->frosch_precond);
          }
          solver_iterative->set_tol_rel(frosch_gmres_params._tol_rel);
          solver_iterative->set_tol_abs(frosch_gmres_params._tol_abs);
          solver_iterative->set_tol_abs_low(DataType(1e-14));
          solver_iterative->set_min_iter(frosch_gmres_params._miniter);
          solver_iterative->set_max_iter(frosch_gmres_params._maxiter);
          // solver_iterative->set_plot_mode(FEAT::Solver::PlotMode::summary);
          solver_iterative->set_plot_mode(coarse_solver_info ? FEAT::Solver::PlotMode::summary : FEAT::Solver::PlotMode::none);

          this->coarse_solver = solver_iterative;
          this->multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, this->coarse_solver);
        }
#endif
        else
        {
          // create FGMRES-AmaVanka coarse grid solver
          std::shared_ptr<FEAT::Solver::AmaVanka<LocalSystemMatrix, LocalSystemFilter>> vanka(nullptr);
          if(voxel_vanka)
          {
            if constexpr(Intern::supports_element_coloring<Domain_>)
              vanka = FEAT::Solver::new_voxel_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local(), domain.at(i)->element_coloring);
            else
              vanka = FEAT::Solver::new_voxel_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local(), Intern::get_element_coloring(*domain.at(i)));
          }
          else
            vanka = FEAT::Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
          vanka->set_skip_singular(true);
          // ama_vankas.push_back(vanka);
          // the coarse grid solver: on the CUDA backend we use the specialized, kernel-fused
          // AmaVankaGMRES (which uses the raw AmaVanka directly), otherwise the standard GMRES
          // with a Schwarz-wrapped AmaVanka preconditioner. Both default to the cgs2 orthogonalizer.
          std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>> coarse_solver_t;
#ifdef FEAT_HAVE_CUDA
          if(solver_backend == FEAT::PreferredBackend::cuda)
          {
            auto amavanka_gmres = FEAT::Solver::new_amavanka_gmres(lvl.matrix_sys, lvl.filter_sys, vanka, coarse_gmres_params._gmres_dim, coarse_gmres_params._inner_rescale);
            amavanka_gmres->set_use_cuda_gmres(true);
            amavanka_gmres->set_orthogonalizer(coarse_gmres_params._ortho);
            coarse_solver_t = amavanka_gmres;
          }
          else
#endif // FEAT_HAVE_CUDA
          {
            auto schwarz = FEAT::Solver::new_schwarz_precond(vanka, lvl.filter_sys);
            schwarz->set_ignore_status(true);
            //auto coarse_solver = FEAT::Solver::new_bicgstab(lvl.matrix_sys, lvl.filter_sys, schwarz);
            auto gmres = FEAT::Solver::new_gmres(lvl.matrix_sys, lvl.filter_sys, coarse_gmres_params._gmres_dim, coarse_gmres_params._inner_rescale, schwarz);
            gmres->set_orthogonalizer(coarse_gmres_params._ortho);
            coarse_solver_t = gmres;
          }
          coarse_solver_t->set_tol_rel(coarse_gmres_params._tol_rel);
          coarse_solver_t->set_tol_abs(coarse_gmres_params._tol_abs);
          coarse_solver_t->set_tol_abs_low(DataType(1e-14));
          coarse_solver_t->set_min_iter(coarse_gmres_params._miniter);
          coarse_solver_t->set_max_iter(coarse_gmres_params._maxiter);
          coarse_solver_t->set_plot_mode(coarse_solver_info ? FEAT::Solver::PlotMode::summary : FEAT::Solver::PlotMode::none);

          this->multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, coarse_solver_t);
          this->coarse_solver = std::move(coarse_solver_t);
          this->vanka_solver.push_back(std::move(vanka));
        }
      }

      this->multigrid = FEAT::Solver::new_multigrid(this->multigrid_hierarchy, cycle);
      // create our solver
      if(solve_gmres_dim > FEAT::Index(0))
        this->base_solver = this->iter_solver = FEAT::Solver::new_fgmres(matrix_sys, filter_sys, solve_gmres_dim, DataType(0.7), this->multigrid);
      else
        this->base_solver = this->iter_solver = FEAT::Solver::new_richardson(matrix_sys, filter_sys, DataType(1.0), this->multigrid);

      // configure iterative solver
      if(this->iter_solver)
      {
        // set solver name
        if(solve_gmres_dim > FEAT::Index(0))
          this->iter_solver->set_plot_name("FGMRES-MG[" + FEAT::stringify(solve_gmres_dim) + "]");
        else
          this->iter_solver->set_plot_name("Multigrid");

        // configure solver
        this->iter_solver->set_min_iter(min_iter);
        this->iter_solver->set_max_iter(max_iter);
        this->iter_solver->set_tol_rel(mg_tol_rel);
        this->iter_solver->set_min_stag_iter(min_stag_iter);
      }
    }

    void init_symbolic()
    {
      if(multigrid_hierarchy) multigrid_hierarchy->init_symbolic();
    }

    void done_symbolic()
    {
      if(multigrid_hierarchy) multigrid_hierarchy->done_symbolic();
    }

    void init_numeric()
    {
      // FEAT::PreferredBackend prev_backend = FEAT::Backend::get_preferred_backend();
      FEAT::Backend::set_preferred_backend(solver_backend);
      watch_mg_hirarch.start();
      if(multigrid_hierarchy) multigrid_hierarchy->init_numeric();
      watch_mg_hirarch.stop();
    }

    void done_numeric()
    {
      // FEAT::PreferredBackend prev_backend = FEAT::Backend::get_preferred_backend();
      FEAT::Backend::set_preferred_backend(solver_backend);
      watch_mg_hirarch.start();
      if(multigrid_hierarchy) multigrid_hierarchy->done_numeric();
      watch_mg_hirarch.stop();
      // FEAT::Backend::set_preferred_backend(prev_backend);
    }

    FEAT::String _format_coarse_gmres_info(unsigned long long padlen) const
    {
      FEAT::String s{};
      s += FEAT::String("---------- Coarse GMRES --------------\n");
      s += FEAT::String("GMRES Dim").pad_back(padlen, '.') + ":" + FEAT::stringify(coarse_gmres_params._gmres_dim) + "\n";
      s += FEAT::String("Min Iters").pad_back(padlen, '.') + ":" + FEAT::stringify(coarse_gmres_params._miniter) + "\n";
      s += FEAT::String("Max Iters").pad_back(padlen, '.') + ":" + FEAT::stringify(coarse_gmres_params._maxiter) + "\n";
      s += FEAT::String("Tol Rel").pad_back(padlen, '.') + ":" + FEAT::stringify_fp_sci(coarse_gmres_params._tol_rel, 3, 5) + "\n";
      s += FEAT::String("Tol Abs").pad_back(padlen, '.') + ":" + FEAT::stringify_fp_sci(coarse_gmres_params._tol_abs, 3, 5) + "\n";
      s += FEAT::String("Inner Rescale").pad_back(padlen, '.') + ":" + FEAT::stringify_fp_fix(coarse_gmres_params._inner_rescale, 3, 5) + "\n";
      s += FEAT::String("Orthogonalizer").pad_back(padlen, '.') + ":" + (coarse_gmres_params._ortho == FEAT::Solver::KrylovOrthogonalizerType::cgs2 ? FEAT::String("cgs2") : FEAT::String("cgs")) + "\n";

      return s;
    }

    FEAT::String format_parameters(unsigned long long padlen) const
    {
      FEAT::String s{};

      s += FEAT::String("-----------------------------------\n Multigrid Component \n-----------------------------------\n");

      s += FEAT::String("Smoother Damp").pad_back(padlen, '.') + ":" + FEAT::stringify_fp_fix(smoother_damp, 3, 5) + "\n";
      s += FEAT::String("Smoother Steps").pad_back(padlen, '.') + ":" + FEAT::stringify(smooth_steps) + "\n";
      s += FEAT::String("Smoother GMRES Dim").pad_back(padlen, '.') + ":" + FEAT::stringify(smooth_gmres_dim) + "\n";
      s += FEAT::String("Solve GMRES Dim").pad_back(padlen, '.') + ":" + FEAT::stringify(solve_gmres_dim) + "\n";
      s += FEAT::String("Min MG Iter").pad_back(padlen, '.') + ":" + FEAT::stringify(min_iter) + "\n";
      s += FEAT::String("Max MG Iter").pad_back(padlen, '.') + ":" + FEAT::stringify(max_iter) + "\n";
      s += FEAT::String("MG Tol Rel").pad_back(padlen, '.') + ":" + FEAT::stringify_fp_sci(mg_tol_rel, 3, 5) + "\n";
      s += FEAT::String("Solver Backend").pad_back(padlen, '.') + ":" + FEAT::stringify(solver_backend) + "\n";
      if(coarse_gmres)
      {
        s += _format_coarse_gmres_info(padlen);
      }

      return s;
    }

    FEAT::String format_timings(const FEAT::Dist::Comm& comm, double total_time, unsigned long long padlen) const
    {
      enum class _TimeID : std::size_t
      {
        mg_hirarch = 0,
        vanka_symbolic_init,
        vanka_numeric_init,
        num_entries
      };
      FEAT::String s;
      double timings_max[to_underlying(_TimeID::num_entries)], timings_min[to_underlying(_TimeID::num_entries)], timings[to_underlying(_TimeID::num_entries)];
      timings_min[to_underlying(_TimeID::mg_hirarch)] = timings_max[to_underlying(_TimeID::mg_hirarch)] = timings[to_underlying(_TimeID::mg_hirarch)] = this->watch_mg_hirarch.elapsed();
      timings_min[to_underlying(_TimeID::vanka_symbolic_init)] = timings_max[to_underlying(_TimeID::vanka_symbolic_init)] = timings[to_underlying(_TimeID::vanka_symbolic_init)] = std::accumulate(this->vanka_solver.begin(), this->vanka_solver.end(), double(0),
                              [](double t, const auto& van){return t + van->time_init_symbolic();});
      timings_min[to_underlying(_TimeID::vanka_numeric_init)] = timings_max[to_underlying(_TimeID::vanka_numeric_init)] = timings[to_underlying(_TimeID::vanka_numeric_init)] = std::accumulate(this->vanka_solver.begin(), this->vanka_solver.end(), double(0),
                              [](double t, const auto& van){return t + van->time_init_numeric();});
      // sync our timings
      comm.allreduce(timings_max, timings_max, to_underlying(_TimeID::num_entries), FEAT::Dist::op_max);
      comm.allreduce(timings_min, timings_min, to_underlying(_TimeID::num_entries), FEAT::Dist::op_min);

      s += format_subtime_mm("Create MG Hirarch", timings[to_underlying(_TimeID::mg_hirarch)], total_time, timings_min[to_underlying(_TimeID::mg_hirarch)], timings_max[to_underlying(_TimeID::mg_hirarch)], padlen);
      s += format_subtime_mm("Vanka Symbolic Init", timings[to_underlying(_TimeID::vanka_symbolic_init)], total_time, timings_min[to_underlying(_TimeID::vanka_symbolic_init)], timings_max[to_underlying(_TimeID::vanka_symbolic_init)], padlen);
      s += format_subtime_mm("Vanka Numeric Init", timings[to_underlying(_TimeID::vanka_numeric_init)], total_time, timings_min[to_underlying(_TimeID::vanka_numeric_init)], timings_max[to_underlying(_TimeID::vanka_numeric_init)], padlen);
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += this->_format_multigrid_timings();
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");

      return s;
    }

    FEAT::String _format_multigrid_timings() const
    {
      using String = FEAT::String;
      String s;

      if(!multigrid_hierarchy)
        return s;

      s = "Multigrid Timings:\n";
      s += "              Defect /   Smoother /   Transfer /     Coarse\n";
      s += "Overall : " +
        FEAT::stringify_fp_fix(multigrid_hierarchy->get_time_defect(), 3, 10) + " / " +
        FEAT::stringify_fp_fix(multigrid_hierarchy->get_time_smooth(), 3, 10) + " / " +
        FEAT::stringify_fp_fix(multigrid_hierarchy->get_time_transfer(), 3, 10) + " / " +
        FEAT::stringify_fp_fix(multigrid_hierarchy->get_time_coarse(), 3, 10) + "\n";
      for(int i(0); i < int(multigrid_hierarchy->size_physical()); ++i)
      {
        s += "Level " + FEAT::stringify(i).pad_front(2) + ": " +
          FEAT::stringify_fp_fix(multigrid_hierarchy->get_time_defect(i), 3, 10) + " / " +
          FEAT::stringify_fp_fix(multigrid_hierarchy->get_time_smooth(i), 3, 10) + " / " +
          FEAT::stringify_fp_fix(multigrid_hierarchy->get_time_transfer(i), 3, 10) + " / " +
          FEAT::stringify_fp_fix(multigrid_hierarchy->get_time_coarse(i), 3, 10) + "\n";
      }

      return s;
    }


  public:
    MultigridVankaComponent() :
#ifdef FEAT_HAVE_TRILINOS
        _frosch_subregions({1}),
        _frosch_overlap(1),
        _frosch_directsolver_coarse(FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK),
        _frosch_directsolver_overlapping(FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK),
        _frosch_directsolver_extension(FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK),
        _frosch_combine_lvl(FEAT::Solver::Trilinos::FROSchParameterList::CombineLevel::ADDITIVE),
        _frosch_precond_type(FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner::TWOLEVELBLOCK),
        _frosch_parti(FEAT::Solver::Trilinos::FROSchParameterList::PartitionType::PARMETIS),
        _frosch_approach(FEAT::Solver::Trilinos::FROSchParameterList::PartitionApproach::REPARTITION),
        _frosch_combine_overlap(FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap::RESTRICTED),
        //_frosch_ipou_velo(dim == 3 ? FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSWSTAR : FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSW),
        _frosch_ipou_velo(FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSWSTAR),
        _frosch_ipou_pres(FEAT::Solver::Trilinos::FROSchParameterList::IPOU::RGDSW),
        _frosch_verbose(false),
        _frosch_exclude_velocity_pressure(false),
        _frosch_exclude_pressure_velocity(false),
        _frosch_reuse_sf_ao(false),
        _frosch_reuse_sf_cm(false),
        _frosch_reuse_sf_ext(false),
        _frosch_reuse_cm(false),
        _frosch_reuse_cb(false),
        _frosch_xml(),
        _frosch_use_timers(false),
        _frosch_print_internal(false),
#endif
        domain_virtual_size(std::size_t(0)),
        smoother_damp(DataType(0.5)),
        mg_tol_rel(DataType(0)),
        smooth_gmres_dim(0u),
        solve_gmres_dim(4u),
        smooth_steps(12u),
        min_iter(1u),
        max_iter(50u),
        min_stag_iter(3u),
        solver_backend(FEAT::PreferredBackend::generic),
        cycle(FEAT::Solver::MultiGridCycle::V),
        coarse_frosch(false),
        coarse_gmres(false),
        voxel_vanka(false),
        coarse_solver_info(false)
    {}

    ~MultigridVankaComponent() = default;
    MultigridVankaComponent(const MultigridVankaComponent&) = default;
    MultigridVankaComponent(MultigridVankaComponent&&) = default;
    MultigridVankaComponent& operator=(const MultigridVankaComponent&) = default;
    MultigridVankaComponent& operator=(MultigridVankaComponent&&) = default;

    template<typename Domain_, typename System_, template<typename> typename Container_>
    MultigridVankaComponent(const Container_<std::shared_ptr<System_>>& sys_levels, const Domain_& domain, const FEAT::PropertyMap* prop) :
#ifdef FEAT_HAVE_TRILINOS
        _frosch_subregions(1, 1),
        _frosch_overlap(1),
        _frosch_directsolver_coarse(FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK),
        _frosch_directsolver_overlapping(FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK),
        _frosch_directsolver_extension(FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK),
        _frosch_combine_lvl(FEAT::Solver::Trilinos::FROSchParameterList::CombineLevel::ADDITIVE),
        _frosch_precond_type(FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner::TWOLEVELBLOCK),
        _frosch_parti(FEAT::Solver::Trilinos::FROSchParameterList::PartitionType::PARMETIS),
        _frosch_approach(FEAT::Solver::Trilinos::FROSchParameterList::PartitionApproach::REPARTITION),
        _frosch_combine_overlap(FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap::RESTRICTED),
        //_frosch_ipou_velo(dim == 3 ? FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSWSTAR : FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSW),
        _frosch_ipou_velo(FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSWSTAR),
        _frosch_ipou_pres(FEAT::Solver::Trilinos::FROSchParameterList::IPOU::RGDSW),
        _frosch_verbose(false),
        _frosch_exclude_velocity_pressure(false),
        _frosch_exclude_pressure_velocity(false),
        _frosch_reuse_sf_ao(false),
        _frosch_reuse_sf_cm(false),
        _frosch_reuse_sf_ext(false),
        _frosch_reuse_cm(false),
        _frosch_reuse_cb(false),
        _frosch_xml(),
        _frosch_use_timers(false),
        _frosch_print_internal(false),
#endif
        domain_virtual_size(domain.size_virtual()),
        smoother_damp(DataType(0.5)),
        mg_tol_rel(DataType(0)),
        smooth_gmres_dim(0u),
        solve_gmres_dim(4u),
        smooth_steps(12u),
        min_iter(1u),
        max_iter(50u),
        min_stag_iter(3u),
        solver_backend(FEAT::PreferredBackend::generic),
        cycle(FEAT::Solver::MultiGridCycle::V),
        coarse_frosch(false),
        coarse_gmres((domain.back_layer().comm().size() > 1)),
        voxel_vanka(false),
        coarse_solver_info(false)
    {
      this->parse(prop);
      this->set_solver(sys_levels, domain);
    }

  }; //class MultigridVankaComponent

  template<typename LevelType_, bool use_fbm_ = false>
  class MultigridVankaFlowSolver : public SteadyStokesFlowSolverBaseCRTP<MultigridVankaFlowSolver<LevelType_, use_fbm_>, LevelType_>
  {
  public:
    typedef SteadyStokesFlowSolverBaseCRTP<MultigridVankaFlowSolver, LevelType_> BaseClass;
    typedef typename BaseClass::SystemLevelType SystemLevelType;
    typedef typename BaseClass::SystemLevels SystemLevels;
    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::DefectVectorType DefectVectorType;
    typedef typename BaseClass::VectorType VectorType;

    friend BaseClass;

    typedef typename SystemLevelType::GlobalSystemMatrix SystemMatrix;
    typedef typename SystemLevelType::LocalSystemMatrix LocalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemTransfer SystemTransfer;
    typedef typename SystemLevelType::GlobalSystemFilter SystemFilter;
    typedef typename SystemLevelType::LocalSystemFilter LocalSystemFilter;

    static constexpr int dim = SystemLevelType::dim;
    static constexpr bool use_fbm = use_fbm_;

    MultigridVankaComponent<SystemMatrix, SystemTransfer, SystemFilter> multigrid_solver;

    std::size_t domain_virtual_size;
    FEAT::PreferredBackend solver_backend;
    mutable FEAT::StopWatch watch_create_solver, watch_system_filter;

  protected:
    const Gendie::Logger* _logger;
    bool _init_params;

    DefectVectorType _create_vector() const
    {
      return this->system_levels.front()->create_global_vector_sys();
    }

    bool _parse(Intern::SolverParser<DataType>& parser)
    {
      bool success = true;
      success &= BaseClass::_parse(parser);
      solver_backend = Intern::backend_from_string(*parser.backend);
      success &= multigrid_solver.parse(parser);

      return success;
    }

    void _clear_solver()
    {
      this->base_solver.reset();
      this->iter_solver.reset();
      multigrid_solver.clear_solver();
    }

    template<typename Domain_>
    void _set_solver(const Domain_& domain)
    {
      watch_create_solver.start();
      XASSERTM(this->_init_params, "You have to initialize the system before calling this function");
      XASSERTM(this->system_levels.size() > 0, "System levels not set");
      this->_clear_solver();
      multigrid_solver.template set_solver<>(this->system_levels, domain);
      // multigrid_solver.template set_solver<decltype(this->system_levels), decltype(domain)>(this->system_levels, domain);
      this->base_solver = multigrid_solver.base_solver;
      this->iter_solver = multigrid_solver.iter_solver;
      watch_create_solver.stop();
    }

    void _init_symbolic()
    {
      // in any case, we need to compile the local type1 matrix for our preconditioners
      for(std::size_t i = 0; i < this->system_levels.size(); ++i)
      {
        auto& system = *this->system_levels.at(i);
        system.compile_system_matrix();
        system.compile_local_matrix();
      }
      multigrid_solver.init_symbolic();
      BaseClass::_init_symbolic();
    }

    void _done_symbolic()
    {
      BaseClass::_done_symbolic();
      multigrid_solver.done_symbolic();
      // free up type 1 matrix
      for(std::size_t i = 0; i < this->system_levels.size(); ++i)
      {
        auto& system = *this->system_levels.at(i);
        system.local_matrix_sys = typename SystemLevelType::LocalSystemMatrix();
      }

    }

    void _init_numeric()
    {
      // FEAT::PreferredBackend prev_backend = FEAT::Backend::get_preferred_backend();
      FEAT::Backend::set_preferred_backend(solver_backend);
      // before the numeric initialization, we have to guarantee that our system is compiled
      watch_system_filter.start();
      for(std::size_t i = 0; i < this->system_levels.size(); ++i)
      {
        auto& system = *this->system_levels.at(i);
        // this solver always assumes fbm
        if constexpr(use_fbm_)
          system.filter_interface_fbm.filter_weak_matrix_rows(system.matrix_a.local(), system.velo_mass_matrix.local());
        system.compile_system_matrix();
        system.compile_local_matrix();
      }
      watch_system_filter.stop();
      multigrid_solver.init_numeric();
      BaseClass::_init_numeric();
      // FEAT::Backend::set_preferred_backend(prev_backend);
    }

    void _done_numeric()
    {
      // FEAT::PreferredBackend prev_backend = FEAT::Backend::get_preferred_backend();
      FEAT::Backend::set_preferred_backend(solver_backend);
      BaseClass::_done_numeric();
      multigrid_solver.done_numeric();
      // FEAT::Backend::set_preferred_backend(prev_backend);
    }

    FEAT::Solver::Status _apply(DefectVectorType& vec_cor, const DefectVectorType& vec_def) const
    {
      // FEAT::PreferredBackend prev_backend = FEAT::Backend::get_preferred_backend();

      FEAT::Backend::set_preferred_backend(solver_backend);
      // todo: explicitly transfer vec_cor and vec_def memory to our device memory?
      FEAT::Solver::Status status = BaseClass::_apply(vec_cor, vec_def);
      // FEAT::Backend::set_preferred_backend(prev_backend);
      return status;
    }

    FEAT::String _format_timings() const
    {
      enum class _TimeID : std::size_t
      {
        total = 0,
        solver_apply,
        symbolic_init,
        create_solver,
        numeric_init,
        system_filter,
        num_entries
      };
      FEAT::String s;
      double timings_max[to_underlying(_TimeID::num_entries)], timings_min[to_underlying(_TimeID::num_entries)], timings[to_underlying(_TimeID::num_entries)];
      timings_min[to_underlying(_TimeID::total)] = timings_max[to_underlying(_TimeID::total)] = timings[to_underlying(_TimeID::total)] = this->watch_total.elapsed();
      timings_min[to_underlying(_TimeID::solver_apply)] = timings_max[to_underlying(_TimeID::solver_apply)] = timings[to_underlying(_TimeID::solver_apply)] = this->watch_apply.elapsed();
      timings_min[to_underlying(_TimeID::symbolic_init)] = timings_max[to_underlying(_TimeID::symbolic_init)] = timings[to_underlying(_TimeID::symbolic_init)] = this->watch_init_symbolic.elapsed();
      timings_min[to_underlying(_TimeID::create_solver)] = timings_max[to_underlying(_TimeID::create_solver)] = timings[to_underlying(_TimeID::create_solver)] = watch_create_solver.elapsed();
      timings_min[to_underlying(_TimeID::numeric_init)] = timings_max[to_underlying(_TimeID::numeric_init)] = timings[to_underlying(_TimeID::numeric_init)] = this->watch_init_numeric.elapsed();
      timings_min[to_underlying(_TimeID::system_filter)] = timings_max[to_underlying(_TimeID::system_filter)] = timings[to_underlying(_TimeID::system_filter)] = watch_system_filter.elapsed();
      // sync our timings
      _logger->comm.allreduce(timings_max, timings_max, to_underlying(_TimeID::num_entries), FEAT::Dist::op_max);
      _logger->comm.allreduce(timings_min, timings_min, to_underlying(_TimeID::num_entries), FEAT::Dist::op_min);

      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Timings MG Vanka Solver---------------------------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += format_subtime_mm("Filter Local System", timings[to_underlying(_TimeID::system_filter)], timings_max[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::system_filter)], timings_max[to_underlying(_TimeID::system_filter)], this->padlen);
      s += format_subtime_mm("Create Solver", timings[to_underlying(_TimeID::create_solver)], timings_max[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::create_solver)], timings_max[to_underlying(_TimeID::create_solver)], this->padlen);
      s += format_subtime_mm("Symbolic Init", timings[to_underlying(_TimeID::symbolic_init)], timings_max[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::symbolic_init)], timings_max[to_underlying(_TimeID::symbolic_init)], this->padlen);
      s += format_subtime_mm("Numeric Init", timings[to_underlying(_TimeID::numeric_init)], timings_max[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::numeric_init)], timings_max[to_underlying(_TimeID::numeric_init)], this->padlen);
      s += format_subtime_mm("Apply Linear Solver", timings[to_underlying(_TimeID::solver_apply)], timings_max[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::solver_apply)], timings_max[to_underlying(_TimeID::solver_apply)], this->padlen);
      s += format_subtime_mm("Linear Solver Total Time", timings_max[to_underlying(_TimeID::total)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::total)], timings_max[to_underlying(_TimeID::total)], this->padlen);
      s += multigrid_solver.format_timings(_logger->comm, timings_max[to_underlying(_TimeID::total)], this->padlen);

      return s;
    }

    FEAT::String _format_parameters() const
    {
      FEAT::String s{};
      return this->multigrid_solver.format_parameters(20);
    }

  public:
    MultigridVankaFlowSolver()
     :  BaseClass(),
        multigrid_solver(),
        domain_virtual_size(std::size_t(0)),
        solver_backend(FEAT::PreferredBackend::generic),
        _logger(nullptr),
        _init_params(false)
    {}

    ~MultigridVankaFlowSolver() = default;
    MultigridVankaFlowSolver(const MultigridVankaFlowSolver&) = default;
    MultigridVankaFlowSolver(MultigridVankaFlowSolver&&) = default;
    MultigridVankaFlowSolver& operator=(const MultigridVankaFlowSolver&) = default;
    MultigridVankaFlowSolver& operator=(MultigridVankaFlowSolver&&) = default;

    template<typename Domain_>
    MultigridVankaFlowSolver(SystemLevels sys_levels, const Domain_& domain, const FEAT::PropertyMap* prop, const Logger* logger = nullptr)
      : BaseClass(std::move(sys_levels)),
        multigrid_solver(),
        domain_virtual_size(domain.size_virtual()),
        solver_backend(FEAT::PreferredBackend::generic),
        _logger(logger),
        _init_params(true)
    {
      if(prop != nullptr)
      {
        Intern::SolverParser<DataType> parser;
        XASSERTM(parser.parse(*prop), "Could not parse solver section");
        _parse(parser);
      }
      this->_set_solver(domain);
    }

    template<typename Domain_>
    void set_system(SystemLevels sys_levels, const Domain_& domain)
    {
      this->system_levels = std::move(sys_levels);
      domain_virtual_size = domain.size_virtual();
      _init_params = true;
      this->_set_solver(domain);
    }

    void set_logger(const Logger* log)
    {
      _logger = log;
    }
  }; //class MultigridVankaFlowSolver

}
