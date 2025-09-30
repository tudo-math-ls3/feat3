// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/backend.hpp>

#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/solver/voxel_amavanka.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/richardson.hpp>
#include "logger.hpp"
#include <applications/gendie/parsing_helper.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/string.hpp>

#include <kernel/solver/direct_stokes_solver.hpp>
// #include <kernel/solver/frosch.hpp>
#include <control/scalar_basic.hpp>
#include "scalexa_gendie_scalarize_helper.hpp"
#include "format_helper.hpp"

#include <memory>

namespace Gendie
{
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

    bool _parse(const FEAT::PropertyMap* DOXY(prop))
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

    double get_total_time_elapsed() const
    {
      return this->watch_total.elapsed();
    }
  }; // class SteadyFlowSolverBaseCRTP<...>

  template<typename LevelType_>
  class MultigridVankaFBMFlowSolver : public SteadyStokesFlowSolverBaseCRTP<MultigridVankaFBMFlowSolver<LevelType_>, LevelType_>
  {
  public:
    typedef SteadyStokesFlowSolverBaseCRTP<MultigridVankaFBMFlowSolver, LevelType_> BaseClass;
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

    std::shared_ptr<FEAT::Solver::MultiGridHierarchy<SystemMatrix, SystemFilter, SystemTransfer>> multigrid_hierarchy;
    // hold a deque of the local vanka solvers
    std::deque<std::shared_ptr<FEAT::Solver::AmaVanka<LocalSystemMatrix, LocalSystemFilter>>> vanka_solver;
    std::deque<std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>>> smoother;
    std::shared_ptr<FEAT::Solver::SolverBase<DefectVectorType>> coarse_solver;
    std::shared_ptr<FEAT::Solver::MultiGrid<SystemMatrix, SystemFilter, SystemTransfer>> multigrid;
    #ifdef FEAT_HAVE_TRILINOS
    // define our scalarized solver level for the coarse grid solver
    typedef FEAT::Control::ScalarUnitFilterSystemLevel<double, typename LocalSystemMatrix::IndexType> ScalarizedSystemLevelType;
    // create scalarize helper
    std::shared_ptr<Gendie::GendieScalarizeHelper<SystemLevelType, ScalarizedSystemLevelType>> scalarize_helper;
    std::shared_ptr<FEAT::Solver::StokesFROSchPreconditioner<typename ScalarizedSystemLevelType::GlobalSystemMatrix, typename ScalarizedSystemLevelType::GlobalSystemFilter>> frosch_precond;
    std::shared_ptr<FEAT::Solver::Trilinos::FROSchParameterList> frosch_params;
    std::vector<int> _frosch_subregions;
    FEAT::Index _frosch_gmres_dim;
    FEAT::Index _frosch_miniter;
    FEAT::Index _frosch_maxiter;
    int _frosch_nlevels;
    int _frosch_overlap;
    FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver _frosch_directsolver;
    FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner _frosch_precond_type;
    FEAT::Solver::Trilinos::FROSchParameterList::PartitionType _frosch_parti;
    FEAT::Solver::Trilinos::FROSchParameterList::PartitionApproach _frosch_approach;
    FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap _frosch_combine_overlap;
    FEAT::Solver::Trilinos::FROSchParameterList::IPOU _frosch_ipou_velo;
    FEAT::Solver::Trilinos::FROSchParameterList::IPOU _frosch_ipou_pres;
    FEAT::String _frosch_xml;
    bool _frosch_use_timers;
    bool _frosch_print_internal;
    #endif

    std::size_t domain_virtual_size;
    DataType smoother_damp;
    DataType mg_tol_rel;
    DataType coarse_tol_rel;
    FEAT::Index smooth_gmres_dim;
    FEAT::Index solve_gmres_dim;
    FEAT::Index smooth_steps;
    FEAT::Index coarse_max_steps;
    FEAT::Index min_iter;
    FEAT::Index max_iter;
    FEAT::Index min_stag_iter;
    FEAT::PreferredBackend solver_backend;
    FEAT::Solver::MultiGridCycle cycle;
    mutable FEAT::StopWatch watch_create_solver, watch_system_filter, watch_mg_hirarch, watch_vanka_symbolic_init, watch_vanka_numeric_init;
    bool coarse_frosch;
    bool coarse_gmres;
    bool voxel_vanka;
    bool coarse_solver_info;

  protected:
    const Gendie::Logger* _logger;
    bool _init_params;

    DefectVectorType _create_vector() const
    {
      return this->system_levels.front()->create_global_vector_sys();
    }

    bool _parse(const FEAT::PropertyMap* prop)
    {
      bool success = true;
      success &= BaseClass::_parse(prop);
      if(prop == nullptr)
        return success;

      success &= prop->parse_entry("min-mg-iter", min_iter);
      success &= prop->parse_entry("max-mg-iter", max_iter);
      success &= prop->parse_entry("smooth-gmres", smooth_gmres_dim);
      success &= prop->parse_entry("solve-gmres", solve_gmres_dim);
      success &= prop->parse_entry("smooth-steps", smooth_steps);
      success &= prop->parse_entry("smooth-damp", smoother_damp);
      solver_backend = Gendie::parse_backend(prop->query("backend", "generic"));
      coarse_gmres |=  Gendie::check_for_config_option(prop->query("coarse-gmres"));
      coarse_frosch |=  Gendie::check_for_config_option(prop->query("coarse-frosch"));
      coarse_solver_info |=  Gendie::check_for_config_option(prop->query("coarse-info"));

    #ifdef FEAT_HAVE_TRILINOS
      success &= this->_parse_frosch_params(prop->get_sub_section("frosch-parameters"));
    #else
      coarse_frosch = false;
    #endif // FEAT_HAVE_TRILINOS
      return success;
    }

    #ifdef FEAT_HAVE_TRILINOS
    bool _parse_frosch_params(const FEAT::PropertyMap* prop)
    {
      bool success = true;
      if(prop == nullptr)
        return success;

      // _frosch_subregions(1, 1),
      success &= prop->parse_entry("frosch-gmres-dim", _frosch_gmres_dim);
      success &= prop->parse_entry("frosch-min-iter", _frosch_miniter);
      success &= prop->parse_entry("frosch-max-iter", _frosch_maxiter);
      if(Gendie::check_for_config_option(prop->query("frosch-xml")))
      {
        _frosch_xml = prop->query("frosch-xml", "");
      }
      success &= prop->parse_entry("frosch-nlevels", _frosch_nlevels);
      XASSERTM(_frosch_nlevels == 1, "For now, only nlevels == 1 allowed");
      success &= prop->parse_entry("frosch-overlap", _frosch_overlap);

      // parse directsolver
      {
        const auto& val = prop->query("frosch-direct-solver");
        if(val.second)
        {
          if(val.first.compare_no_case("umfpack"))
          {
            _frosch_directsolver = FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK;
          }
          else if(val.first.compare_no_case("mumps"))
          {
            _frosch_directsolver = FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::MUMPS;
          }
          else if(val.first.compare_no_case("ilu"))
          {
            _frosch_directsolver = FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::ILU;
          }
          else if(val.first.compare_no_case("klu"))
          {
            _frosch_directsolver = FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::KLU;
          }
          else
          {
            XABORTM("Unknown parameter " + val.first);
          }
        }
      }
      // parse precond type
      {
        const auto& val = prop->query("frosch-precond-type");
        if(val.second)
        {
          if(val.first.compare_no_case("onelevel"))
          {
            _frosch_precond_type = FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner::ONELEVEL;
          }
          else if(val.first.compare_no_case("twolevel"))
          {
            _frosch_precond_type = FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner::TWOLEVEL;
          }
          else if(val.first.compare_no_case("twolevelblock"))
          {
            _frosch_precond_type = FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner::TWOLEVELBLOCK;
          }
          else
          {
            XABORTM("Unknown parameter " + val.first);
          }
        }
      }
      // parse partition type
      {
        const auto& val = prop->query("frosch-partition-type");
        if(val.second)
        {
          if(val.first.compare_no_case("PHG"))
          {
            _frosch_parti = FEAT::Solver::Trilinos::FROSchParameterList::PartitionType::PHG;
          }
          else if(val.first.compare_no_case("PARMETIS"))
          {
            _frosch_parti = FEAT::Solver::Trilinos::FROSchParameterList::PartitionType::PARMETIS;
          }
          else if(val.first.compare_no_case("BLOCK"))
          {
            _frosch_parti = FEAT::Solver::Trilinos::FROSchParameterList::PartitionType::BLOCK;
          }
          else
          {
            XABORTM("Unknown parameter " + val.first);
          }
        }
      }
      // parse partition approach
      {
        const auto& val = prop->query("frosch-partition-approach");
        if(val.second)
        {
          if(val.first.compare_no_case("partition"))
          {
            _frosch_approach = FEAT::Solver::Trilinos::FROSchParameterList::PartitionApproach::PARTITION;
          }
          else if(val.first.compare_no_case("repartition"))
          {
            _frosch_approach = FEAT::Solver::Trilinos::FROSchParameterList::PartitionApproach::REPARTITION;
          }
          else
          {
            XABORTM("Unknown parameter " + val.first);
          }
        }
      }
      // parse combine overlap
      {
        const auto& val = prop->query("frosch-combine-overlap");
        if(val.second)
        {
          if(val.first.compare_no_case("restricted"))
          {
            _frosch_combine_overlap = FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap::RESTRICTED;
          }
          else if(val.first.compare_no_case("full"))
          {
            _frosch_combine_overlap = FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap::FULL;
          }
          else if(val.first.compare_no_case("averaging"))
          {
            _frosch_combine_overlap = FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap::AVERAGING;
          }
          else
          {
            XABORTM("Unknown parameter " + val.first);
          }
        }
      }
      // parse ipou velocity
      {
        const auto& val = prop->query("frosch-ipou-velocity");
        if(val.second)
        {
          if(val.first.compare_no_case("gdsw"))
          {
            _frosch_ipou_velo = FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSW;
          }
          else if(val.first.compare_no_case("gdswstar"))
          {
            _frosch_ipou_velo = FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSWSTAR;
          }
          else if(val.first.compare_no_case("rgdsw"))
          {
            _frosch_ipou_velo = FEAT::Solver::Trilinos::FROSchParameterList::IPOU::RGDSW;
          }
          else
          {
            XABORTM("Unknown parameter " + val.first);
          }
        }
      }
      // parse ipou pres
      {
        const auto& val = prop->query("frosch-ipou-pressure");
        if(val.second)
        {
          if(val.first.compare_no_case("gdsw"))
          {
            _frosch_ipou_pres = FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSW;
          }
          else if(val.first.compare_no_case("gdswstar"))
          {
            _frosch_ipou_pres = FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSWSTAR;
          }
          else if(val.first.compare_no_case("rgdsw"))
          {
            _frosch_ipou_pres = FEAT::Solver::Trilinos::FROSchParameterList::IPOU::RGDSW;
          }
          else
          {
            XABORTM("Unknown parameter " + val.first);
          }
        }
      }
      // print internal timers and debug output?
      {
        _frosch_use_timers = Gendie::check_for_config_option(prop->query("frosch-use-timers"));
        _frosch_print_internal = Gendie::check_for_config_option(prop->query("frosch-print-internal"));
      }

      return success;
    }

    bool _setup_frosch_params()
    {
      // for now local variables, will be used in class as configuration needs arise
      // for now, always use _frosch_nlvels = 1, else we would have to use vectors for the interface
      // if we have an xml file. we always only parse this
      if(!_frosch_xml.empty())
      {
        if(_logger) _logger->print("Use frosch xml file");
        this->frosch_params->read_from_xml_file(_frosch_xml);
        return true;
      }

      this->frosch_params->set_nlevels(_frosch_nlevels);
      this->frosch_params->set_precond(_frosch_precond_type);
      this->frosch_params->set_coarse_solver(_frosch_directsolver);
      this->frosch_params->set_solvers(_frosch_directsolver, _frosch_directsolver);
      this->frosch_params->set_overlaps(_frosch_overlap);
      this->frosch_params->set_combine_overlap(_frosch_combine_overlap);
      this->frosch_params->set_parti_types(_frosch_parti);
      this->frosch_params->set_parti_approach(_frosch_approach);
      this->frosch_params->set_ipous(_frosch_ipou_velo, _frosch_ipou_pres);
      this->frosch_params->set_subregions(_frosch_subregions);

      this->frosch_params->set_print(_frosch_use_timers);
      this->frosch_params->set_use_timer(_frosch_print_internal);

      return true;
    }
    #endif

    void _clear_solver()
    {
      #ifdef FEAT_HAVE_TRILINOS
      this->frosch_precond.reset();
      this->frosch_params.reset();
      this->scalarize_helper.reset();
      #endif
      this->base_solver.reset();
      this->iter_solver.reset();
      this->multigrid.reset();
      this->multigrid_hierarchy.reset();
      this->vanka_solver.clear();
      this->coarse_solver.reset();
      this->smoother.clear();
    }

    template<typename Domain_>
    void _set_solver(const Domain_& domain)
    {
      watch_create_solver.start();
      XASSERTM(this->_init_params, "You have to initialize the system before calling this function");
      XASSERTM(this->system_levels.size() > 0, "System levels not set");
      this->_clear_solver();
      const auto& matrix_sys = this->system_levels.front()->matrix_sys;
      const auto& filter_sys = this->system_levels.front()->filter_sys;

      // if we use frosch, initialize frosch params now
      #ifdef FEAT_HAVE_TRILINOS
      if(coarse_frosch && (this->system_levels.size() == domain_virtual_size))
      {
        scalarize_helper = std::make_shared<Gendie::GendieScalarizeHelper<SystemLevelType, ScalarizedSystemLevelType>>();
        scalarize_helper->create(*this->system_levels.back());
        // create on same communicater our coarse level uses
        frosch_params = std::make_shared<FEAT::Solver::Trilinos::FROSchParameterList>(domain.back().layer().comm(), 3, FEAT::Solver::Trilinos::FROSchParameterList::SADDLEPOINT);
        // if(_logger) _logger->print("Parsing from " + frosch_xml, info);

        // frosch_params->read_from_xml_file(frosch_xml);
        this->_setup_frosch_params();
        frosch_params->create_core();
      }

      #endif

      if(domain_virtual_size == std::size_t(1))
      {
        // create direct solver or FMGRES?
#if defined(FEAT_HAVE_UMFPACK) || defined(FEAT_HAVE_CUDSS)
        if(!coarse_frosch && !coarse_gmres)
        {
          auto backendt = FEAT::Backend::get_preferred_backend();
          FEAT::Backend::set_preferred_backend(solver_backend);
          if(_logger) _logger->print("Only 1 level chosen; creating single grid solver: DirectStokesSolver", info);
          this->base_solver = FEAT::Solver::new_direct_stokes_solver(matrix_sys, filter_sys);
          FEAT::Backend::set_preferred_backend(backendt);
        }
#endif //FEAT_HAVE_UMFPACK
#ifdef FEAT_HAVE_TRILINOS
#if defined(FEAT_HAVE_UMFPACK) || defined(FEAT_HAVE_CUDSS)
        else if(coarse_frosch)
#else
        if(coarse_frosch)
#endif
        {
          if(_logger) _logger->print("Only 1 level chosen; creating single grid solver: FroschSolver", info);
          Index num_owned_pres_dofs = this->system_levels.front()->gate_pres.get_num_local_dofs();
          this->frosch_precond = FEAT::Solver::new_stokes_frosch(scalarize_helper->scalarized_level.matrix_sys,
            scalarize_helper->scalarized_level.filter_sys, num_owned_pres_dofs, *frosch_params);

          std::shared_ptr<FEAT::Solver::IterativeSolver<typename ScalarizedSystemLevelType::GlobalSystemVector>> solver_iterative;

          if(_frosch_gmres_dim > 0)
            solver_iterative = Solver::new_fgmres(scalarize_helper->scalarized_level.matrix_sys,
              scalarize_helper->scalarized_level.filter_sys, _frosch_gmres_dim, DataType(1.), this->frosch_precond);
          else
            solver_iterative = Solver::new_richardson(scalarize_helper->scalarized_level.matrix_sys,
              scalarize_helper->scalarized_level.filter_sys, DataType(1), this->frosch_precond);
          solver_iterative->set_tol_rel(mg_tol_rel);
          // solver_iterative->set_tol_abs(tol_abs);
          solver_iterative->set_tol_abs_low(1e-14);
          solver_iterative->set_min_iter(min_iter);
          solver_iterative->set_max_iter(max_iter);
          // solver_iterative->set_plot_mode(FEAT::Solver::PlotMode::summary);
          solver_iterative->set_plot_mode(coarse_solver_info ? FEAT::Solver::PlotMode::summary : FEAT::Solver::PlotMode::none);

          this->base_solver = std::make_shared<Gendie::GendieScalarizeWrapper<SystemLevelType, ScalarizedSystemLevelType>>(*scalarize_helper, solver_iterative);
        }
        else
#endif // FEAT_HAVE_TRILINOS
        {
          if(_logger) _logger->print("Only 1 level chosen; creating single grid solver: FGMRES-AmaVanka");
          if(voxel_vanka)
            this->vanka_solver.push_back(FEAT::Solver::new_voxel_amavanka(this->system_levels.front()->local_matrix_sys, filter_sys.local(), domain.front()->element_coloring));
          else
            this->vanka_solver.push_back(FEAT::Solver::new_amavanka(this->system_levels.front()->local_matrix_sys, filter_sys.local()));
          this->vanka_solver.front()->set_skip_singular(true);
          auto schwarz = FEAT::Solver::new_schwarz_precond(this->vanka_solver.front(), filter_sys);
          if(solve_gmres_dim > FEAT::Index(0)) // todo: is this smart?
            this->base_solver = this->iter_solver = FEAT::Solver::new_fgmres(matrix_sys, filter_sys, solve_gmres_dim, 0.4, schwarz);
          else
            this->base_solver = this->iter_solver = FEAT::Solver::new_fgmres(matrix_sys, filter_sys, 16, 0.4, schwarz);

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

      for(std::size_t i = 0; i < this->system_levels.size(); ++i)
      {
        auto& lvl = *this->system_levels.at(i);

        if((i+1) < domain_virtual_size)
        {
          // create Schwarz-AmaVanka
          std::shared_ptr<FEAT::Solver::AmaVanka<LocalSystemMatrix, LocalSystemFilter>> vanka(nullptr);
          if(voxel_vanka)
            vanka = FEAT::Solver::new_voxel_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local(), domain.at(i)->element_coloring);
          else
            vanka = FEAT::Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
          vanka->set_skip_singular(true);
          // ama_vankas.push_back(vanka);
          auto schwarz = FEAT::Solver::new_schwarz_precond(vanka, lvl.filter_sys);
          schwarz->set_ignore_status(true);

          // create smoother: either GMRES or Richardson
          std::shared_ptr<FEAT::Solver::IterativeSolver<DefectVectorType>> smoother_l;
          if(smooth_gmres_dim > FEAT::Index(0))
            smoother_l = FEAT::Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, smooth_gmres_dim, 0.0, schwarz);
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
#if defined(FEAT_HAVE_UMFPACK) || defined(FEAT_HAVE_CUDSS)
        else if(!coarse_frosch && !coarse_gmres)
        {
          if(_logger) _logger->print("Choose direct coarse solver", debug);
          auto backendt = FEAT::Backend::get_preferred_backend();
          FEAT::Backend::set_preferred_backend(solver_backend);
          // create UMFPACK coarse grid solver
          this->coarse_solver = FEAT::Solver::new_direct_stokes_solver(lvl.matrix_sys, lvl.filter_sys);
          this->multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, this->coarse_solver);
          FEAT::Backend::set_preferred_backend(backendt);
        }
#endif //  FEAT_HAVE_UMFPACK
#ifdef FEAT_HAVE_TRILINOS
        else if(coarse_frosch)
        {
          if(_logger) _logger->print("Chose coarse grid solver: FroschSolver", info);
          Index num_owned_pres_dofs = lvl.gate_pres.get_num_local_dofs();
          this->frosch_precond = Solver::new_stokes_frosch(scalarize_helper->scalarized_level.matrix_sys,
            scalarize_helper->scalarized_level.filter_sys, num_owned_pres_dofs, *frosch_params);

          std::shared_ptr<FEAT::Solver::IterativeSolver<typename ScalarizedSystemLevelType::GlobalSystemVector>> solver_iterative;

          if(_frosch_gmres_dim > 0)
            solver_iterative = Solver::new_fgmres(scalarize_helper->scalarized_level.matrix_sys,
              scalarize_helper->scalarized_level.filter_sys, _frosch_gmres_dim, DataType(1.), this->frosch_precond);
          else
            solver_iterative = Solver::new_richardson(scalarize_helper->scalarized_level.matrix_sys,
              scalarize_helper->scalarized_level.filter_sys, DataType(1), this->frosch_precond);

          solver_iterative->set_tol_rel(coarse_tol_rel);
          // solver_iterative->set_tol_abs(tol_abs);
          solver_iterative->set_tol_abs_low(1e-14);
          solver_iterative->set_min_iter(1);
          solver_iterative->set_max_iter(coarse_max_steps);
          solver_iterative->set_plot_mode(FEAT::Solver::PlotMode::summary);
          solver_iterative->set_plot_mode(coarse_solver_info ? FEAT::Solver::PlotMode::summary : FEAT::Solver::PlotMode::none);

          this->coarse_solver = std::make_shared<Gendie::GendieScalarizeWrapper<SystemLevelType, ScalarizedSystemLevelType>>(*scalarize_helper, solver_iterative);
          this->multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, this->coarse_solver);
        }
#endif
        else
        {
          // create FGMRES-AmaVanka coarse grid solver
          std::shared_ptr<FEAT::Solver::AmaVanka<LocalSystemMatrix, LocalSystemFilter>> vanka(nullptr);
          if(voxel_vanka)
            vanka = FEAT::Solver::new_voxel_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local(), domain.at(i)->element_coloring);
          else
            vanka = FEAT::Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());
          vanka->set_skip_singular(true);
          // ama_vankas.push_back(vanka);
          auto schwarz = FEAT::Solver::new_schwarz_precond(vanka, lvl.filter_sys);
          schwarz->set_ignore_status(true);
          //auto coarse_solver = FEAT::Solver::new_bicgstab(lvl.matrix_sys, lvl.filter_sys, schwarz);
          auto coarse_solver_t = FEAT::Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, 16, 0.0, schwarz);
          coarse_solver_t->set_max_iter(coarse_max_steps);
          coarse_solver_t->set_tol_rel(coarse_tol_rel);
          coarse_solver_t->set_plot_mode(coarse_solver_info ? FEAT::Solver::PlotMode::summary : FEAT::Solver::PlotMode::none);

          this->multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, coarse_solver_t);
          this->coarse_solver = std::move(coarse_solver_t);
          this->vanka_solver.push_back(std::move(vanka));
        }
      }

      this->multigrid = FEAT::Solver::new_multigrid(this->multigrid_hierarchy, cycle);
      // create our solver
      if(solve_gmres_dim > FEAT::Index(0))
        this->base_solver = this->iter_solver = FEAT::Solver::new_fgmres(matrix_sys, filter_sys, solve_gmres_dim, 0.7, this->multigrid);
      else
        this->base_solver = this->iter_solver = FEAT::Solver::new_richardson(matrix_sys, filter_sys, 1.0, this->multigrid);

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
      if(multigrid_hierarchy) multigrid_hierarchy->init_symbolic();
      BaseClass::_init_symbolic();
    }

    void _done_symbolic()
    {
      BaseClass::_done_symbolic();
      if(multigrid_hierarchy) multigrid_hierarchy->done_symbolic();
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
        system.filter_interface_fbm.filter_weak_matrix_rows(system.matrix_a.local(), system.velo_mass_matrix.local());
        system.compile_system_matrix();
        system.compile_local_matrix();
      }
      watch_system_filter.stop();
      watch_mg_hirarch.start();
      if(multigrid_hierarchy) multigrid_hierarchy->init_numeric();
      watch_mg_hirarch.stop();
      BaseClass::_init_numeric();
      // FEAT::Backend::set_preferred_backend(prev_backend);
    }

    void _done_numeric()
    {
      // FEAT::PreferredBackend prev_backend = FEAT::Backend::get_preferred_backend();
      FEAT::Backend::set_preferred_backend(solver_backend);
      BaseClass::_done_numeric();
      watch_mg_hirarch.start();
      if(multigrid_hierarchy) multigrid_hierarchy->done_numeric();
      watch_mg_hirarch.stop();
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
      enum _TimeID
      {
        total = 0,
        solver_apply = 1,
        symbolic_init = 2,
        create_solver = 3,
        numeric_init = 4,
        system_filter = 5,
        mg_hirarch = 6,
        vanka_symbolic_init = 7,
        vanka_numeric_init = 8,
        num_entries = 9
      };
      FEAT::String s;
      double timings_max[_TimeID::num_entries], timings_min[_TimeID::num_entries], timings[_TimeID::num_entries];
      timings_min[_TimeID::total] = timings_max[_TimeID::total] = timings[_TimeID::total] = this->watch_total.elapsed();
      timings_min[_TimeID::solver_apply] = timings_max[_TimeID::solver_apply] = timings[_TimeID::solver_apply] = this->watch_apply.elapsed();
      timings_min[_TimeID::symbolic_init] = timings_max[_TimeID::symbolic_init] = timings[_TimeID::symbolic_init] = this->watch_init_symbolic.elapsed();
      timings_min[_TimeID::create_solver] = timings_max[_TimeID::create_solver] = timings[_TimeID::create_solver] = watch_create_solver.elapsed();
      timings_min[_TimeID::numeric_init] = timings_max[_TimeID::numeric_init] = timings[_TimeID::numeric_init] = this->watch_init_numeric.elapsed();
      timings_min[_TimeID::system_filter] = timings_max[_TimeID::system_filter] = timings[_TimeID::system_filter] = watch_system_filter.elapsed();
      timings_min[_TimeID::mg_hirarch] = timings_max[_TimeID::mg_hirarch] = timings[_TimeID::mg_hirarch] = watch_mg_hirarch.elapsed();
      timings_min[_TimeID::vanka_symbolic_init] = timings_max[_TimeID::vanka_symbolic_init] = timings[_TimeID::vanka_symbolic_init] = std::accumulate(this->vanka_solver.begin(), this->vanka_solver.end(), double(0),
                              [](double t, const auto& van){return t + van->time_init_symbolic();});
      timings_min[_TimeID::vanka_numeric_init] = timings_max[_TimeID::vanka_numeric_init] = timings[_TimeID::vanka_numeric_init] = std::accumulate(this->vanka_solver.begin(), this->vanka_solver.end(), double(0),
                              [](double t, const auto& van){return t + van->time_init_numeric();});
      // sync our timings
      _logger->comm.allreduce(timings_max, timings_max, _TimeID::num_entries, FEAT::Dist::op_max);
      _logger->comm.allreduce(timings_min, timings_min, _TimeID::num_entries, FEAT::Dist::op_min);

      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Timings MG Vanka Solver---------------------------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += format_subtime_mm("Create MG Hirarch", timings[_TimeID::mg_hirarch], timings[_TimeID::total], timings_min[_TimeID::mg_hirarch], timings_max[_TimeID::mg_hirarch], this->padlen);
      s += format_subtime_mm("Filter Local System", timings[_TimeID::system_filter], timings[_TimeID::total], timings_min[_TimeID::system_filter], timings_max[_TimeID::system_filter], this->padlen);
      s += format_subtime_mm("Create Solver", timings[_TimeID::create_solver], timings[_TimeID::total], timings_min[_TimeID::create_solver], timings_max[_TimeID::create_solver], this->padlen);
      s += format_subtime_mm("Symbolic Init", timings[_TimeID::symbolic_init], timings[_TimeID::total], timings_min[_TimeID::symbolic_init], timings_max[_TimeID::symbolic_init], this->padlen);
      s += format_subtime_mm("Vanka Symbolic Init", timings[_TimeID::vanka_symbolic_init], timings[_TimeID::total], timings_min[_TimeID::vanka_symbolic_init], timings_max[_TimeID::vanka_symbolic_init], this->padlen);
      s += format_subtime_mm("Numeric Init", timings[_TimeID::numeric_init], timings[_TimeID::total], timings_min[_TimeID::numeric_init], timings_max[_TimeID::numeric_init], this->padlen);
      s += format_subtime_mm("Vanka Numeric Init", timings[_TimeID::vanka_numeric_init], timings[_TimeID::total], timings_min[_TimeID::vanka_numeric_init], timings_max[_TimeID::vanka_numeric_init], this->padlen);
      s += format_subtime_mm("Apply Linear Solver", timings[_TimeID::solver_apply], timings[_TimeID::total], timings_min[_TimeID::solver_apply], timings_max[_TimeID::solver_apply], this->padlen);
      s += format_subtime_mm("Linear Solver Total Time", timings[_TimeID::total], timings[_TimeID::total], timings_min[_TimeID::total], timings_max[_TimeID::total], this->padlen);
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
    MultigridVankaFBMFlowSolver()
     :  BaseClass(),
#ifdef FEAT_HAVE_TRILINOS
        _frosch_subregions(1, 1),
        _frosch_gmres_dim(FEAT::Index(16)),
        _frosch_nlevels(1),
        _frosch_overlap(1),
        _frosch_directsolver(FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK),
        _frosch_precond_type(FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner::TWOLEVELBLOCK),
        _frosch_parti(FEAT::Solver::Trilinos::FROSchParameterList::PartitionType::PARMETIS),
        _frosch_approach(FEAT::Solver::Trilinos::FROSchParameterList::PartitionApproach::REPARTITION),
        _frosch_combine_overlap(FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap::RESTRICTED),
        _frosch_ipou_velo(dim == 3 ? FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSWSTAR : FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSW),
        _frosch_ipou_pres(FEAT::Solver::Trilinos::FROSchParameterList::IPOU::RGDSW),
        _frosch_xml(),
        _frosch_use_timers(false),
        _frosch_print_internal(false),
#endif
        domain_virtual_size(std::size_t(0)),
        smoother_damp(DataType(0.5)),
        mg_tol_rel(DataType(0)),
        coarse_tol_rel(DataType(1E-3)),
        smooth_gmres_dim(0u),
        solve_gmres_dim(4u),
        smooth_steps(12u),
        coarse_max_steps(500u),
        min_iter(1u),
        max_iter(50u),
        min_stag_iter(3u),
        solver_backend(FEAT::PreferredBackend::generic),
        cycle(FEAT::Solver::MultiGridCycle::V),
        coarse_frosch(false),
        coarse_gmres(false),
        voxel_vanka(false),
        coarse_solver_info(false),
        _logger(nullptr),
        _init_params(false)
    {}

    ~MultigridVankaFBMFlowSolver() = default;
    MultigridVankaFBMFlowSolver(const MultigridVankaFBMFlowSolver&) = default;
    MultigridVankaFBMFlowSolver(MultigridVankaFBMFlowSolver&&) = default;
    MultigridVankaFBMFlowSolver& operator=(const MultigridVankaFBMFlowSolver&) = default;
    MultigridVankaFBMFlowSolver& operator=(MultigridVankaFBMFlowSolver&&) = default;

    template<typename Domain_>
    MultigridVankaFBMFlowSolver(SystemLevels sys_levels, const Domain_& domain, const FEAT::PropertyMap* prop, const Logger* logger = nullptr)
      : BaseClass(std::move(sys_levels)),
#ifdef FEAT_HAVE_TRILINOS
        _frosch_subregions(1, 1),
        _frosch_gmres_dim(FEAT::Index(16)),
        _frosch_nlevels(1),
        _frosch_overlap(1),
        _frosch_directsolver(FEAT::Solver::Trilinos::FROSchParameterList::DirectSolver::UMFPACK),
        _frosch_precond_type(FEAT::Solver::Trilinos::FROSchParameterList::Preconditioner::TWOLEVELBLOCK),
        _frosch_parti(FEAT::Solver::Trilinos::FROSchParameterList::PartitionType::PARMETIS),
        _frosch_approach(FEAT::Solver::Trilinos::FROSchParameterList::PartitionApproach::REPARTITION),
        _frosch_combine_overlap(FEAT::Solver::Trilinos::FROSchParameterList::CombineOverlap::RESTRICTED),
        _frosch_ipou_velo(dim == 3 ? FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSWSTAR : FEAT::Solver::Trilinos::FROSchParameterList::IPOU::GDSW),
        _frosch_ipou_pres(FEAT::Solver::Trilinos::FROSchParameterList::IPOU::RGDSW),
        _frosch_xml(),
        _frosch_use_timers(false),
        _frosch_print_internal(false),
#endif
        domain_virtual_size(domain.size_virtual()),
        smoother_damp(DataType(0.5)),
        mg_tol_rel(DataType(0)),
        coarse_tol_rel(DataType(1E-3)),
        smooth_gmres_dim(0u),
        solve_gmres_dim(4u),
        smooth_steps(12u),
        coarse_max_steps(500u),
        min_iter(1u),
        max_iter(50u),
        min_stag_iter(3u),
        solver_backend(FEAT::PreferredBackend::generic),
        cycle(FEAT::Solver::MultiGridCycle::V),
        coarse_frosch(false),
        coarse_gmres((domain.back_layer().comm().size() > 1)),
        voxel_vanka(false),
        coarse_solver_info(false),
        _logger(logger),
        _init_params(true)
    {
      _parse(prop);
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
  }; //class MultigridVankaFBMFlowSolver

}
