// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include <kernel/solver/base.hpp>
#include <kernel/util/property_map.hpp>
#include "logger.hpp"
#include "format_helper.hpp"
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/string.hpp>
#include <applications/gendie/gendie_common.hpp>

#include <memory>

namespace Gendie
{
  typedef FEAT::Index Index;
  /**
   * \brief BaseClass for steady-solver interface
   */
  template<typename VectorType_>
  class NonlinearSteadyFlowSolverBase
  {
  public:
    typedef VectorType_ DefectVectorType;

    NonlinearSteadyFlowSolverBase() = default;
    virtual ~NonlinearSteadyFlowSolverBase() = default;
    NonlinearSteadyFlowSolverBase(const NonlinearSteadyFlowSolverBase&) = delete;
    NonlinearSteadyFlowSolverBase(NonlinearSteadyFlowSolverBase&&) = default;
    NonlinearSteadyFlowSolverBase& operator=(const NonlinearSteadyFlowSolverBase&) = delete;
    NonlinearSteadyFlowSolverBase& operator=(NonlinearSteadyFlowSolverBase&&) = default;

    /**
     * \brief Parse property map to configure the solver
     *
     * \param[in] prop Property map to be parsed, interpreted based on actual implementation
     *
     * \returns True if parsing was sucessfull
     */
    virtual bool parse(const FEAT::PropertyMap* prop) = 0;

    virtual FEAT::String format_string() const = 0;

    virtual void init() = 0;
    // virtual void init_symbolic() = 0;
    // virtual void init_numeric() = 0;
    virtual void done() = 0;
    // virtual void done_symbolic() = 0;
    // virtual void done_numeric() = 0;

    virtual void reset() = 0;

    /**
     * \brief Applies the solver to a rhs handside
     *
     * This applies a nonlinear flow solver to the rhs. The solution steps are applied in a defect correction manner,
     * therefore rhs and solution do need to be initialized with the correct boundary values.
     *
     * \param[inout] vec_sol The solution vector, has to be assembled to correct size and boundary values (i.e. type 1). The initial value is
     *                 NOT ignored, so it is advantegeous if a good initial guess is given.
     * \param[in] vec_rhs The rhs vector with correct boundary values. Has to be synchronized, i.e. type 1
     */
    virtual FEAT::Solver::Status apply(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs) = 0;

    virtual bool init_sol(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs) = 0;

    virtual FEAT::String format_timings() const = 0;

    virtual double get_assembly_time() const = 0;

    virtual double get_linear_solver_time() const = 0;
  };

  namespace Intern
  {
    template<typename DT_>
    struct DataTypeId;

    template<>
    struct DataTypeId<double>
    {
      static constexpr char id[] = "fp_64";
    };

    template<>
    struct DataTypeId<float>
    {
      static constexpr char id[] = "fp_32";
    };
  }

  template<typename Derived_, typename StokesSolver_, typename DefectAssembler_, typename SystemAssembler_, typename DefectFilter_>
  class NonlinearSteadyFlowSolverCRTP : public NonlinearSteadyFlowSolverBase<typename DefectAssembler_::VectorType>
  {
  public:
    typedef NonlinearSteadyFlowSolverBase<typename DefectAssembler_::VectorType> BaseClass;
    typedef typename BaseClass::DefectVectorType DefectVectorType;
    typedef StokesSolver_ SteadyStokesSolver;
    typedef DefectAssembler_ DefectAssembler;
    typedef SystemAssembler_ SystemAssembler;
    typedef DefectFilter_ DefectFilter;
    typedef typename SteadyStokesSolver::VectorType SolverVectorType;
    typedef DefectFilter FilterType;
    typedef typename SteadyStokesSolver::VectorType::DataType SolverDataType;
    typedef typename DefectVectorType::DataType DefectDataType;

    static constexpr bool no_convert = std::is_same_v<DefectVectorType, SolverVectorType>;
    const String solver_datatype_id = String(Intern::DataTypeId<SolverDataType>::id);
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

    bool _parse(const FEAT::PropertyMap* prop)
    {
      bool success = true;
      {
        if(prop == nullptr)
        {
          if(_logger)
            _logger->print("Propertymap is null!", error);
          return false;
        }

        success &= prop->parse_entry("min-nl-iter", _min_nonlin_iter);
        success &= prop->parse_entry("max-nl-iter", _max_nonlin_iter);
        success &= prop->parse_entry("nl-tol-abs", _tol_abs);
        success &= prop->parse_entry("nl-tol-rel", _tol_rel);
        success &= prop->parse_entry("nl-tol-low-abs", _tol_low_abs);
        success &= prop->parse_entry("nl-stagrate", _stag_rate);
        success &= prop->parse_entry("fixed-ls-tol", _fixed_ls_tol);
        _solve_stokes = Gendie::check_for_config_option(prop->query("init-with-stokes"));
      }

      return success;
    }

    FEAT::String _format_string() const
    {
      FEAT::String s;

      s += FEAT::String("Solver Datatype ").pad_back(padlen, pc) + ": " + solver_datatype_id + "\n";
      s += FEAT::String("Min Nonlin Iter ").pad_back(padlen, pc) + ": " + FEAT::stringify(_min_nonlin_iter) + "\n";
      s += FEAT::String("Max Nonlin Iter ").pad_back(padlen, pc) + ": " + FEAT::stringify(_max_nonlin_iter) + "\n";
      s += FEAT::String("Nonlin Tol Abs ").pad_back(padlen, pc) + ": " + FEAT::stringify_fp_sci(_tol_abs, 2, 4) + "\n";
      s += FEAT::String("Nonlin Tol Rel ").pad_back(padlen, pc) + ": " + FEAT::stringify_fp_sci(_tol_rel, 2, 4) + "\n";
      s += FEAT::String("Nonlin Tol Low Abs ").pad_back(padlen, pc) + ": " + FEAT::stringify_fp_sci(_tol_low_abs, 2, 4) + "\n";
      s += FEAT::String("Nonlin Stagrate ").pad_back(padlen, pc) + ": " + FEAT::stringify_fp_sci(_stag_rate, 2, 4) + "\n";
      s += FEAT::String("Fixed Linear Solver tolerance ").pad_back(padlen, pc) + ": " +
           ((_fixed_ls_tol < DefectDataType(0)) ? FEAT::String("adaptive") : FEAT::stringify_fp_sci(_fixed_ls_tol, 2, 4)) + "\n";


      return s;

    }

    void _init()
    {
      this->flow_solver->init_symbolic();
      this->defect_asm->init();
      this->system_asm->init();

      // is this really necessary?
      // XASSERTM(&(this->flow_solver->get_domain()) == &(this->defect_asm->get_domain()), "Solver and assembler do not share domain");
      // XASSERTM(&(this->flow_solver->get_domain()) == &(this->system_asm->get_domain()), "Solver and assembler do not share domain");
      this->_vec_def = this->defect_asm->create_vector();
      this->_vec_cor = this->defect_asm->create_vector();

      // probably not required...
      // this->_system_matrix = this->system_asm.create_system_matrix();

      // is the solver vectortype the same as the defect vectortype?
      if constexpr(no_convert)
      {
        this->_vec_solver_def = this->_vec_def.clone(FEAT::LAFEM::CloneMode::Shallow);
        this->_vec_solver_cor = this->_vec_cor.clone(FEAT::LAFEM::CloneMode::Shallow);
      }
      else
      {
        this->_vec_solver_def = this->flow_solver->create_defect();
        this->_vec_solver_cor = this->flow_solver->create_correction();
      }
    }

    void _done()
    {
      this->_vec_def.clear();
      this->_vec_cor.clear();
      if constexpr(!no_convert)
      {
        this->_vec_solver_cor.clear();
        this->_vec_solver_def.clear();
      }
      // this->_system_matrix.clear();
      this->system_asm->done();
      this->defect_asm->done();
      this->flow_solver->done_symbolic();
    }

    FEAT::String _format_timings() const
    {
      enum _TimeID
      {
        total_t = 0,
        scaling_t = 1,
        init_sol_t = 2,
        num_entries_t = 3
      };
      FEAT::String s;
      double timings_max[num_entries_t], timings_min[num_entries_t], timings[num_entries_t];
      timings_min[_TimeID::total_t] = timings_max[_TimeID::total_t] = timings[_TimeID::total_t] = watch_total.elapsed();
      timings_min[_TimeID::scaling_t] = timings_max[_TimeID::scaling_t] = timings[_TimeID::scaling_t] = watch_internal_scaling.elapsed();
      timings_min[_TimeID::init_sol_t] = timings_max[_TimeID::init_sol_t] = timings[_TimeID::init_sol_t] = watch_init_sol.elapsed();
      // sync our timings
      _logger->comm.allreduce(timings_max, timings_max, _TimeID::num_entries_t, FEAT::Dist::op_max);
      _logger->comm.allreduce(timings_min, timings_min, _TimeID::num_entries_t, FEAT::Dist::op_min);


      s += this->defect_asm->format_timings();
      s += this->system_asm->format_timings();
      s += this->flow_solver->format_timings();
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Timings Flow Solver-------------------------------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");

      s += format_subtime_mm("Internal Pressure Scaling", timings[_TimeID::scaling_t], timings[_TimeID::total_t], timings_min[_TimeID::scaling_t], timings_max[_TimeID::scaling_t], this->padlen);
      s += format_subtime_mm("Init Solution", timings[_TimeID::init_sol_t], timings[_TimeID::total_t], timings_min[_TimeID::init_sol_t], timings_max[_TimeID::init_sol_t], this->padlen);
      s += format_subtime_mm("Total Time", timings[_TimeID::total_t], timings[_TimeID::total_t], timings_min[_TimeID::total_t], timings_max[_TimeID::total_t], this->padlen);

      return s;
    }

    //for now does nothing
    void _reset()
    {
    }

    bool _init_sol(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs)
    {
      // in any case, format vec_sol
      bool success = true;
      vec_sol.format();
      if(this->_solve_stokes)
      {
        if(this->_logger)
        {
          this->_logger->print("\nInit solution vector by solving stokes", info);
        }

        this->system_asm->set_backend();
        // we solve "stokes" by simply assembling with a zero convection vector field
        this->system_asm->assemble_matrices_isotrop(this->flow_solver->get_systems());

        // // we use a homogenous rhs
        this->_vec_cor.local().template at<0>().scale(vec_rhs.local().template at<0>(), DefectDataType(1)/this->_scaling_factor);
        this->_vec_cor.local().template at<1>().scale(vec_rhs.local().template at<1>(), DefectDataType(1));

        //filter
        this->filter.filter_rhs(this->_vec_cor);
        this->filter.filter_sol(vec_sol);


        this->flow_solver->set_plot_mode(FEAT::Solver::PlotMode::iter);

        // this->flow_solver->set_backend();

        this->flow_solver->init_numeric();

        // we only want to solve by 3 orders of magnitude or at maximum 10 solver iterations
        this->flow_solver->set_tol_abs(1E+100);
        this->flow_solver->set_tol_rel(1E-3);
        auto prev_max_iter = this->flow_solver->get_max_iter();
        this->flow_solver->set_max_iter(8);

        FEAT::Solver::Status status;
        if constexpr(no_convert)
        {
          status = this->flow_solver->solve(vec_sol, this->_vec_cor);
        }
        else
        {
          this->_vec_solver_cor.local().convert(vec_sol.local());
          this->_vec_solver_def.local().convert(this->_vec_cor.local());
          status = this->flow_solver->solve(this->_vec_solver_cor, this->_vec_solver_def);
          vec_sol.local().convert(this->_vec_solver_cor.local());
        }
        // rescale pressure part
        this->watch_internal_scaling.start();
        vec_sol.local().template at<1>().scale(vec_sol.local().template at<1>(), this->_scaling_factor);
        this->watch_internal_scaling.stop();

        this->flow_solver->set_max_iter(prev_max_iter);

        if(!FEAT::Solver::status_success(status))
        {
          this->_logger->print("Stokes solver diverged, refomatting with zero", warning);
          vec_sol.format();
          this->filter.filter_sol(vec_sol);
          success = false;
        }

        this->flow_solver->set_plot_mode(FEAT::Solver::PlotMode::none);

        this->flow_solver->done_numeric();
      }
      else
      {
        this->filter.filter_sol(vec_sol);
      }

      return success;
    }

    double _get_assembly_time() const
    {
      return system_asm->get_total_time_elapsed() + defect_asm->get_total_time_elapsed();
    }

    double _get_linear_solver_time() const
    {
      return flow_solver->get_total_time_elapsed();
    }

    const Gendie::Logger* _logger;
    DefectVectorType _vec_def, _vec_cor;
    SolverVectorType _vec_solver_def, _vec_solver_cor;
    DefectDataType _tol_abs, _tol_rel, _tol_low_abs;
    DefectDataType _stag_rate;
    DefectDataType _fixed_ls_tol;
    DefectDataType _scaling_factor;
    Index _min_nonlin_iter, _max_nonlin_iter, _cur_iter, _max_stokes_iter;

    std::deque<DefectDataType> _defects;

    FEAT::StopWatch watch_total, watch_init_sol, watch_internal_scaling;

    bool _solve_stokes;


  public:
    std::shared_ptr<SteadyStokesSolver> flow_solver;
    std::shared_ptr<DefectAssembler> defect_asm;
    std::shared_ptr<SystemAssembler> system_asm;
    DefectFilter filter;

    // rule of 5
    NonlinearSteadyFlowSolverCRTP() = default;
    virtual ~NonlinearSteadyFlowSolverCRTP() = default;
    NonlinearSteadyFlowSolverCRTP(const NonlinearSteadyFlowSolverCRTP&) = delete;
    NonlinearSteadyFlowSolverCRTP& operator=(const NonlinearSteadyFlowSolverCRTP&) = delete;
    NonlinearSteadyFlowSolverCRTP(NonlinearSteadyFlowSolverCRTP&&) = default;
    NonlinearSteadyFlowSolverCRTP& operator=(NonlinearSteadyFlowSolverCRTP&&) = default;

    /**
     * \brief Constructs Solver out of shared and unique ptr to our internal solver and assembler
     *
     * \param[in] f_solver_ A shared ptr provding a linear solver interface
     * \param[in] defect_asm_ A shared ptr to an assembler interface that can assemble the defect for a given primal and convective vector
     * \param[in] system_asm_ A shared ptr to an assembler interface that can assemble the system matrix for a given convective vector
     * \param[in] filter_ A reference to the defect filter, which is shallow cloned.
     */
    NonlinearSteadyFlowSolverCRTP(std::shared_ptr<SteadyStokesSolver> f_solver_, std::shared_ptr<DefectAssembler> defect_asm_, std::shared_ptr<SystemAssembler> system_asm_,
                                  DefectFilter& filter_, const Gendie::Logger* logger = nullptr)
     : BaseClass(),
       _logger(logger),
      //  _system_matrix(),
       _vec_def(),
       _vec_cor(),
       _vec_solver_def(),
       _vec_solver_cor(),
       _tol_abs(DefectDataType(1E-1)),
       _tol_rel(DefectDataType(1E-4)),
       _tol_low_abs(DefectDataType(1E-8)),
       _stag_rate(DefectDataType(0.98)),
       _fixed_ls_tol(DefectDataType(-1)),
       _scaling_factor(DefectDataType(1)),
       _min_nonlin_iter(Index(1)),
       _max_nonlin_iter(Index(20)),
       _cur_iter(Index(0)),
       _max_stokes_iter(Index(8)),
       _defects(),
       _solve_stokes(false),
       flow_solver(std::move(f_solver_)),
       defect_asm(std::move(defect_asm_)),
       system_asm(std::move(system_asm_)),
       filter(filter_.clone(FEAT::LAFEM::CloneMode::Shallow))
    {
    }

    virtual bool parse(const FEAT::PropertyMap* prop) override final
    {
      return this->cast()._parse(prop);
    }

    virtual FEAT::String format_string() const override final
    {
      return this->cast()._format_string();
    }

    virtual void init() override final
    {
      this->watch_total.start();
      this->cast()._init();
      this->watch_total.stop();
    }

    // // TODO: necessary to split this? Not really, since we only really need symbolic init
    // virtual void init_symbolic() override final
    // {
    //   this->cast()._symbolic_init();
    // }

    // virtual void init_numeric() override final
    // {
    //   this->cast()._numeric_init();
    // }

    virtual void done() override final
    {
      this->watch_total.start();
      this->cast()._done();
      this->watch_total.stop();
    }

    // virtual void done_symbolic() override final
    // {
    //   this->cast()._symbolic_done();
    // }

    // virtual void done_numeric() override final
    // {
    //   this->cast()._numeric_done();
    // }

    virtual void reset() override final
    {
      this->cast()._reset();
    }

    virtual FEAT::Solver::Status apply(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs) override final
    {
      this->watch_total.start();
      auto status = this->cast()._apply(vec_sol, vec_rhs);
      this->watch_total.stop();
      return status;
    }

    virtual bool init_sol(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs) override final
    {
      this->watch_total.start();
      this->watch_init_sol.start();
      bool success = this->cast()._init_sol(vec_sol, vec_rhs);
      this->watch_init_sol.stop();
      this->watch_total.stop();
      return success;
    }

    virtual FEAT::String format_timings() const override final
    {
      return this->cast()._format_timings();
    }

    DefectDataType get_scaling_factor() const
    {
      return this->_scaling_factor;
    }

    void set_scaling_factor(DefectDataType scal_fac)
    {
      this->_scaling_factor = scal_fac;
    }

    virtual double get_assembly_time() const override final
    {
      return this->cast()._get_assembly_time();
    }

    virtual double get_linear_solver_time() const override final
    {
      return this->cast()._get_linear_solver_time();
    }

  }; // CRTP NonLinFlowSolver Baseclass

  template<typename StokesSolver_, typename DefectAssembler_, typename SystemAssembler_, typename DefectFilter_>
  class AlPiNeSteadyFlowSolver : public NonlinearSteadyFlowSolverCRTP<AlPiNeSteadyFlowSolver<StokesSolver_, DefectAssembler_, SystemAssembler_, DefectFilter_>, StokesSolver_, DefectAssembler_, SystemAssembler_, DefectFilter_>
  {
  public:
    typedef NonlinearSteadyFlowSolverCRTP<AlPiNeSteadyFlowSolver, StokesSolver_, DefectAssembler_, SystemAssembler_, DefectFilter_> BaseClass;
    typedef typename BaseClass::DefectVectorType DefectVectorType;
    typedef typename BaseClass::SteadyStokesSolver SteadyStokesSolver;
    typedef typename BaseClass::DefectAssembler DefectAssembler;
    typedef typename BaseClass::SystemAssembler SystemAssembler;
    typedef typename BaseClass::DefectFilter DefectFilter;
    typedef typename BaseClass::SolverVectorType SolverVectorType;
    typedef typename BaseClass::FilterType FilterType;
    typedef typename BaseClass::SolverDataType SolverDataType;
    typedef typename BaseClass::DefectDataType DefectDataType;

    friend BaseClass;

    //use constructors of baseclass
    using BaseClass::BaseClass;

    using BaseClass::no_convert;

    std::size_t min_picard_steps = std::size_t(3);

    // AlPiNeSteadyFlowSolver(const std::shared_ptr<SteadyStokesSolver>& f_solver_, const std::shared_ptr<DefectAssembler>& defect_asm_, const std::shared_ptr<SystemAssembler>& system_asm_,
    //                               DefectFilter& filter_, const Gendie::Logger* logger = nullptr)
    //  : BaseClass(f_solver_, defect_asm_, system_asm_, filter_, logger)
    // {}

  protected:

    bool _parse(const FEAT::PropertyMap* prop)
    {
      bool success = BaseClass::_parse(prop);
      if(!prop)
        return success;
      success &= prop->parse_entry("min-picard-steps", min_picard_steps);
      // {
      // }
      return success;
    }

    FEAT::String _format_string() const
    {
      using String = FEAT::String;
      String s = String("AlPine-Solver:\n") + BaseClass::_format_string();
      return s;
    }

    void _init()
    {
      BaseClass::_init();
    }

    DefectDataType _choose_linear_tolerance(FEAT::Index nl_step, DefectDataType def_nl, DefectDataType def_improve, const std::deque<DefectDataType>& defects) const
    {
      DefectDataType abs_tol = DefectDataType(1e+20);
      if(this->_fixed_ls_tol > DefectDataType(0))
      {
        return abs_tol;
      }
      if(nl_step == 0)
      {
        // shutoff absolute tolerance
      }
      else if(nl_step <= min_picard_steps)
      {
        // linear improvement
        abs_tol = def_nl * def_improve * DefectDataType(0.1);
      }
      else if(nl_step == 1)
      {
        // We furthermore limit this absolute tolerance to ensure that we do not
        // overshoot the mark by overoptimistic quadratic convergence expectations.
        abs_tol = def_nl * def_improve * def_improve * DefectDataType(0.1);
      }
      else if((nl_step % 2 > 0) && (nl_step > 3))
      {
        // We furthermore limit this absolute tolerance to ensure that we do not
        // overshoot the mark by overoptimistic quadratic convergence expectations.
        DefectDataType def_prev_imp = defects.at(defects.size()-2)/defects.at(defects.size()-3);
        abs_tol = def_nl * def_prev_imp * def_prev_imp * def_improve * DefectDataType(0.1);
      }
      else
      {
        DefectDataType def_prev_imp = def_improve;
        if(nl_step > 3u)
          def_prev_imp = defects.at(defects.size()-2)/defects.at(defects.size()-3);
        abs_tol = def_nl * def_prev_imp * DefectDataType(0.1);
      }
      return FEAT::Math::max(abs_tol, this->_tol_abs * DefectDataType(0.01));

    }

    FEAT::Solver::Status _apply(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs)
    {
      this->watch_total.start();
      this->_defects.clear();

      FEAT::Solver::Status nl_status = FEAT::Solver::Status::success;

      //scale down starting solution pressure part
      this->watch_internal_scaling.start();
      vec_sol.local().template at<1>().scale(vec_sol.local().template at<1>(), DefectDataType(1)/this->_scaling_factor);
      this->watch_internal_scaling.stop();


      if(this->_logger)
      {
        this->_logger->print("Using scaling factor: " + stringify_fp_sci(this->_scaling_factor) , info);
        this->_logger->print("\nSolver   #  Defect (abs)   Defect (rel)   Improve   |  LS  fin abs Def   fin rel Def  abs Tol", info);
        this->_logger->print(  "----------------------------------------------------+---------------------------------------------", info);
      }
      for(Index nl_step = 0; nl_step < this->_max_nonlin_iter; ++nl_step)
      {
        // this->defect_asm->set_backend();
        //assemble defect
        this->defect_asm->assemble_defect(this->_vec_def, vec_sol, vec_sol, vec_rhs);

        // and filter
        this->filter.filter_def(this->_vec_def);

        if constexpr(!no_convert)
        {
          this->_vec_solver_def.local().convert(this->_vec_def.local());
        }

        const DefectDataType def_prev = (this->_defects.empty() ? DefectDataType(1) : this->_defects.back());
        const DefectDataType def_nl = this->_vec_def.norm2();
        const DefectDataType def_improve = def_nl / def_prev;
        this->_defects.push_back(def_nl);
        const DefectDataType def_rel = def_nl / this->_defects.front();

        FEAT::String line = "AlPiNe: ";
        line += FEAT::stringify(nl_step).pad_front(2) + ": ";
        line += FEAT::stringify_fp_sci(def_nl, 6) + " / ";
        line += FEAT::stringify_fp_sci(def_nl/this->_defects.front()) + " / ";
        line += FEAT::stringify_fp_sci(def_nl/def_prev, 3);

        if(def_rel > DefectDataType(1E+3))
        {
          if(this->_logger)
          {
            this->_logger->print(line, info);
            this->_logger->print("NONLINEAR SOLVER DIVERGED !!!", error);
          }
          //if(testmode)
            //comm.print("Test-Mode: FAILED");
          return FEAT::Solver::Status::diverged;
        }
        else if(nl_step < this->_min_nonlin_iter)
        {
          // nothing to do here; this else-case exists just to ensure
          // that none of the following cases fires and breaks the loop
        }
        else if((def_nl < this->_tol_abs && def_rel < this->_tol_rel) || def_nl < this->_tol_low_abs)
        {
          if(this->_logger)
            this->_logger->print(line + "\nNonlinear solver converged!\n", info);
          nl_status = FEAT::Solver::Status::success;
          break;
        }
        else if(nl_step >= this->_max_nonlin_iter)
        {
          if(this->_logger)
            this->_logger->print(line + "\nMaximum iterations reached!\n", warning);
          nl_status = FEAT::Solver::Status::max_iter;
          break;
        }
        else if((nl_step >= 3) && (nl_step % 2 == 0u) && (this->_stag_rate < def_improve))
        {
          if(this->_logger)
            this->_logger->print(line + "\nNonlinear solver stagnated!\n", warning);
          nl_status = FEAT::Solver::Status::stagnated;
          break;
        }

        // this->system_asm->set_backend();
        // set newton or picard step
        this->system_asm->set_jacobian(double(nl_step%2u & (nl_step>min_picard_steps)));

        // this call has to work with arbitrary solution vector...
        // also assembles matrix for adjusted system (u, \tilde(p)),  \tilde(p) = p/scaling_factor
        this->system_asm->assemble_matrices(this->flow_solver->get_systems(), vec_sol);

        #ifdef FEAT_DEBUG_MODE
        {
          auto diag_a = this->flow_solver->get_systems().front()->matrix_a.create_vector_l();
          auto diag_m = this->flow_solver->get_systems().front()->velo_mass_matrix.create_vector_l();
          this->flow_solver->get_systems().front()->matrix_a.extract_diag(diag_a, true);
          this->flow_solver->get_systems().front()->velo_mass_matrix.extract_diag(diag_m, true);
          auto sc_factor = diag_a.max_abs_element() / diag_m.max_abs_element();
          if(this->_logger) this->_logger->print("Max abs ele " + stringify(sc_factor), info);
        }
        #endif

        // this->flow_solver->set_backend();

        //todo: init_numeric sets its own backend??
        this->flow_solver->init_numeric();

        this->flow_solver->set_tol_abs(SolverDataType(this->_choose_linear_tolerance(nl_step, def_nl, def_improve, this->_defects)));
        this->flow_solver->set_tol_rel(this->_fixed_ls_tol > DefectDataType(0) ? SolverDataType(this->_fixed_ls_tol) : SolverDataType(1E-2));

        if constexpr(!no_convert)
        {
          this->_vec_solver_def.local().convert(this->_vec_def.local());
        }
        FEAT::Solver::Status status = this->flow_solver->apply(this->_vec_solver_cor, this->_vec_solver_def);
        if constexpr(!no_convert)
        {
          this->_vec_cor.local().convert(this->_vec_solver_cor.local());
        }

        // get solver statistics
        if(this->_logger)
        {
          line += FEAT::String(" | ") + FEAT::stringify(this->flow_solver->get_num_iter()).pad_front(3) + ": "
            + FEAT::stringify_fp_sci(this->flow_solver->get_def_final(), 4) + " / "
            + stringify_fp_sci(this->flow_solver->get_def_final() / this->flow_solver->get_def_initial(), 4);
          if((this->_fixed_ls_tol<DefectDataType(0)) && (nl_step > Index(0)))
            line += FEAT::String(" [") + stringify_fp_sci(this->flow_solver->get_tol_abs(), 4) + "]";
          this->_logger->print(line, info);
        }

        this->flow_solver->done_numeric();

        if(!FEAT::Solver::status_success(status))
        {
          if(this->_logger)
            this->_logger->print("\n LINEAR SOLVER BREAKDOWN\n", error);
          return status;
        }

        // update solution
        vec_sol.axpy(this->_vec_cor, DefectDataType(1));

        if(this->_logger)
          this->_logger->flush_print();

      } // end of inner iteration

      //scale back pressure
      this->watch_internal_scaling.start();
      vec_sol.local().template at<1>().scale(vec_sol.local().template at<1>(), this->_scaling_factor);
      this->watch_internal_scaling.stop();

      return nl_status;
    }

  }; // class AlPiNeSteadyFlowSolver<...>

  template<typename StokesSolver_, typename DefectAssembler_, typename SystemAssembler_, typename DefectFilter_>
  class PseudoUnsteadyFlowSolver : public NonlinearSteadyFlowSolverCRTP<PseudoUnsteadyFlowSolver<StokesSolver_, DefectAssembler_, SystemAssembler_, DefectFilter_>, StokesSolver_, DefectAssembler_, SystemAssembler_, DefectFilter_>
  {
  public:
    typedef NonlinearSteadyFlowSolverCRTP<PseudoUnsteadyFlowSolver, StokesSolver_, DefectAssembler_, SystemAssembler_, DefectFilter_> BaseClass;
    typedef typename BaseClass::DefectVectorType DefectVectorType;
    typedef typename BaseClass::SteadyStokesSolver SteadyStokesSolver;
    typedef typename BaseClass::DefectAssembler DefectAssembler;
    typedef typename BaseClass::SystemAssembler SystemAssembler;
    typedef typename BaseClass::DefectFilter DefectFilter;
    typedef typename BaseClass::SolverVectorType SolverVectorType;
    typedef typename BaseClass::FilterType FilterType;
    typedef typename BaseClass::SolverDataType SolverDataType;
    typedef typename BaseClass::DefectDataType DefectDataType;

    friend BaseClass;

    //use constructors of baseclass
    using BaseClass::BaseClass;

    using BaseClass::no_convert;

    /// The initial time step size
    DefectDataType time_step_size = DefectDataType(1e-3);
    /// Multiplicative factor to increase time-step-size if no linear growth is recognized
    DefectDataType adaptive_factor = DefectDataType(4);
    /// Tolerance to decide on linear convergence
    DefectDataType stable_rel_tolerance = DefectDataType(0.05);
    /// Sample size for stable convergence check
    Index sample_size = Index(4);
    /// Do we use adaptive step size ?
    bool adp_step_size = true;

    // PseudoUnsteadyFlowSolver(const std::shared_ptr<SteadyStokesSolver>& f_solver_, const std::shared_ptr<DefectAssembler>& defect_asm_, const std::shared_ptr<SystemAssembler>& system_asm_,
    //                               DefectFilter& filter_, const Gendie::Logger* logger = nullptr)
    //  : BaseClass(f_solver_, defect_asm_, system_asm_, filter_, logger)
    // {}

    void set_steps_size(DefectDataType _step_size)
    {
      time_step_size = _step_size;
    }

    void set_adaptiv_step_size(bool adaptivity)
    {
      adp_step_size = adaptivity;
    }

  protected:

    bool _parse(const FEAT::PropertyMap* prop)
    {
      bool success = BaseClass::_parse(prop);

      adp_step_size = Gendie::check_for_config_option(prop->query("adaptive-step-size"), true);
      success &= prop->parse_entry("pseudo-ts-size", time_step_size);
      success &= prop->parse_entry("adapt-factor", adaptive_factor);
      success &= prop->parse_entry("stability-tol", stable_rel_tolerance);
      success &= prop->parse_entry("stability-sample-size", sample_size);
      return success;
    }

    FEAT::String _format_string() const
    {
      using String = FEAT::String;
      String s = String("PseudoTimeStepping-Solver:\n") + BaseClass::_format_string();
      s += String("Adapt Timestep ").pad_back(this->padlen, this->pc) + ": " + (adp_step_size ? String("Yes") : String("No")) + "\n";
      s += String("Timestep Size ").pad_back(this->padlen, this->pc) + ": " + stringify_fp_sci(time_step_size, 2, 4) + "\n";
      if(adp_step_size)
      {
        s += String("Stablitiy Tolerance ").pad_back(this->padlen, this->pc) + ": " + stringify_fp_sci(stable_rel_tolerance, 2, 4) + "\n";
        s += String("Stablility Sample Size ").pad_back(this->padlen, this->pc) + ": " + stringify(sample_size) + "\n";
        s += String("Adapt Change Factor ").pad_back(this->padlen, this->pc) + ": " + stringify_fp_sci(adaptive_factor, 2, 4) + "\n";
      }
      return s;
    }

    void _init()
    {
      BaseClass::_init();
    }

    DefectDataType _choose_linear_tolerance(FEAT::Index t_step, DefectDataType def_nl, DefectDataType def_improve, [[maybe_unused]] const std::deque<DefectDataType>& defects) const
    {
      DefectDataType abs_tol = DefectDataType(1e+20);
      if(this->_fixed_ls_tol > DefectDataType(0))
      {
        return abs_tol;
      }
      if(t_step == 0)
      {
        // shutoff absolute tolerance
      }
      else //always use a linear estimate for convergence
      {
        // linear improvement
        abs_tol = def_nl * def_improve * DefectDataType(0.1);
      }
      return FEAT::Math::max(abs_tol, this->_tol_abs * DefectDataType(0.01));

    }

    // TODO: Improve this...
    bool _convergence_stable(FEAT::Index t_step, [[maybe_unused]] DefectDataType def_nl, DefectDataType def_improve, const std::deque<DefectDataType>& defects) const
    {
      Index min_steps = this->_min_nonlin_iter;
      // if we do not have enough data points, simply return true
      if(t_step <= min_steps)
        return true;
      // check if we converge worse than the last four data points
      Index unstable_counter = 0;
      for(Index i = 1; i < sample_size; ++i)
      {
        if((t_step -i) <= min_steps)
          break;
        XASSERT(i <= t_step);
        DefectDataType prev_improve = defects[defects.size()-i]/defects[defects.size()-i-1];
        unstable_counter += Index(def_improve > (1+stable_rel_tolerance)*prev_improve);
      }

      Index actual_sample_size = Math::min(t_step-min_steps, sample_size);
      return unstable_counter <= Index(0.35*double(actual_sample_size));
    }

    FEAT::Solver::Status _apply(DefectVectorType& vec_sol, const DefectVectorType& vec_rhs)
    {
      this->_defects.clear();

      FEAT::Solver::Status nl_status = FEAT::Solver::Status::success;

      // for now, we operate on our internal correction vector
      this->_vec_cor.local().copy(vec_sol.local());

      auto& vec_sol_tmp = this->_vec_cor;

      //scale down starting solution pressure part
      this->watch_internal_scaling.start();
      vec_sol_tmp.local().template at<1>().scale(vec_sol_tmp.local().template at<1>(), DefectDataType(1)/this->_scaling_factor);
      this->watch_internal_scaling.stop();

      if constexpr(!no_convert)
      {
        this->_vec_solver_cor.local().convert(vec_sol_tmp.local());
      }

      if(this->_logger)
      {
        this->_logger->print("Using scaling factor: " + stringify_fp_sci(this->_scaling_factor) , info);
        this->_logger->print("\nSolver   #  Defect (abs)   Defect (rel)   Improve   |  LS  fin abs Def   fin rel Def  abs Tol", info);
        this->_logger->print(  "----------------------------------------------------+---------------------------------------------", info);
      }
      for(Index t_step = 0; t_step < this->_max_nonlin_iter; ++t_step)
      {
        // this->defect_asm->set_backend();

        // assemble defect for breaking criteria
        this->defect_asm->assemble_defect(this->_vec_def, vec_sol_tmp, vec_sol_tmp, vec_rhs);

        // and filter
        this->filter.filter_def(this->_vec_def);


        // calculate defect for our stopping criteria
        const DefectDataType def_prev = (this->_defects.empty() ? DefectDataType(1) : this->_defects.back());
        const DefectDataType def_nl = this->_vec_def.norm2();
        const DefectDataType def_improve = def_nl / def_prev;
        this->_defects.push_back(def_nl);
        const DefectDataType def_rel = def_nl / this->_defects.front();

        FEAT::String line = "Pseudo Time Stepping: ";
        line += FEAT::stringify(t_step).pad_front(2) + ": ";
        line += FEAT::stringify_fp_sci(def_nl, 6) + " / ";
        line += FEAT::stringify_fp_sci(def_nl/this->_defects.front()) + " / ";
        line += FEAT::stringify_fp_sci(def_nl/def_prev, 3);

        // increase time step factor if we do not have linear convergence
        if(!_convergence_stable(t_step, def_nl, def_improve, this->_defects))
        {
          time_step_size *= adaptive_factor;
          if(this->_logger) this->_logger->print("Increase time step size to " + stringify_fp_sci(time_step_size, 2), info);
        }

        if(def_rel > DefectDataType(1E+3))
        {
          if(this->_logger)
          {
            this->_logger->print(line, info);
            this->_logger->print("TIMESTEPPING SCHEMA DIVERGED !!!", error);
          }
          //if(testmode)
            //comm.print("Test-Mode: FAILED");
          return FEAT::Solver::Status::diverged;
        }
        else if(t_step < this->_min_nonlin_iter)
        {
          // nothing to do here; this else-case exists just to ensure
          // that none of the following cases fires and breaks the loop
        }
        else if((def_nl < this->_tol_abs && def_rel < this->_tol_rel) || def_nl < this->_tol_low_abs)
        {
          if(this->_logger)
            this->_logger->print(line + "\nNonlinear solver converged!\n", info);
          nl_status = FEAT::Solver::Status::success;
          break;
        }
        else if(t_step >= this->_max_nonlin_iter)
        {
          if(this->_logger)
            this->_logger->print(line + "\nMaximum iterations reached!\n", warning);
          nl_status = FEAT::Solver::Status::max_iter;
          break;
        }
        else if((t_step >= 3) && (this->_stag_rate < def_improve))
        {
          if(this->_logger)
            this->_logger->print(line + "\nNonlinear solver stagnated!\n", warning);
          nl_status = FEAT::Solver::Status::stagnated;
          break;
        }

        // assembly of actual rhs (implicit euler)
        {
          this->_vec_def.local().copy(vec_rhs.local());
          // this->_vec_def.local().format();
          const auto& velo_gate = this->defect_asm->system_level->gate_velo;
          // we have to create temporary global velocity vectors
          typename DefectAssembler::ConvVectorType loc_def(&velo_gate, this->_vec_def.local().
                                                    template at<0>().clone(FEAT::LAFEM::CloneMode::Shallow));
          typename DefectAssembler::ConvVectorType loc_sol(&velo_gate, vec_sol_tmp.local().
                                                    template at<0>().clone(FEAT::LAFEM::CloneMode::Shallow));

          // we need the velocity mass matrix to be assembled on the finest level
          this->defect_asm->system_level->velo_mass_matrix.apply(loc_def, loc_sol, loc_def, DefectDataType(1)/time_step_size);

          // only works with fbm based assembler
          this->defect_asm->system_level->apply_fbm_filter_to_rhs(loc_def.local());
          // and filter
          this->filter.filter_rhs(this->_vec_def);

          // scale rhs
          this->_vec_def.local().template at<0>().scale(this->_vec_def.local().template at<0>(), DefectDataType(1)/this->_scaling_factor);
        }

        // this->system_asm->set_backend();
        // set newton or picard step
        this->system_asm->set_theta(SolverDataType(1)/SolverDataType(time_step_size));

        // always use simple linearization...
        this->system_asm->set_jacobian(double(0));

        // this call has to work with arbitrary solution vector...
        // also assembles matrix for adjusted system (u, \tilde(p)),  \tilde(p) = p/scaling_factor
        this->system_asm->assemble_matrices(this->flow_solver->get_systems(), vec_sol_tmp);

        // this->flow_solver->set_backend();

        // todo: init numeric sets its own backend?
        this->flow_solver->init_numeric();

        this->flow_solver->set_tol_abs(SolverDataType(this->_choose_linear_tolerance(t_step, def_nl, def_improve, this->_defects)));
        this->flow_solver->set_tol_rel(this->_fixed_ls_tol > DefectDataType(0) ? SolverDataType(this->_fixed_ls_tol) : SolverDataType(1E-2));

        if constexpr(!no_convert)
        {
          this->_vec_solver_def.local().convert(this->_vec_def.local());
        }
        FEAT::Solver::Status status = this->flow_solver->solve(this->_vec_solver_cor, this->_vec_solver_def);
        if constexpr(!no_convert)
        {
          this->_vec_cor.local().convert(this->_vec_solver_cor.local());
        }

        // get solver statistics
        if(this->_logger)
        {
          line += FEAT::String(" | ") + FEAT::stringify(this->flow_solver->get_num_iter()).pad_front(3) + ": "
            + FEAT::stringify_fp_sci(this->flow_solver->get_def_final(), 4) + " / "
            + stringify_fp_sci(this->flow_solver->get_def_final() / this->flow_solver->get_def_initial(), 4);
          if((this->_fixed_ls_tol<DefectDataType(0)) && (t_step > Index(0)))
            line += FEAT::String(" [") + stringify_fp_sci(this->flow_solver->get_tol_abs(), 4) + "]";
          this->_logger->print(line, info);
        }

        this->flow_solver->done_numeric();

        if(!FEAT::Solver::status_success(status))
        {
          if(this->_logger)
            this->_logger->print("\n LINEAR SOLVER BREAKDOWN\n", error);
          return status;
        }

        if(this->_logger)
          this->_logger->flush_print();

      } // end of inner iteration

      // finally, copy over solution from our internal vector
      vec_sol.local().copy(this->_vec_cor.local());
      //scale back pressure

      this->watch_internal_scaling.start();
      vec_sol.local().template at<1>().scale(this->_vec_cor.local().template at<1>(), this->_scaling_factor);
      this->watch_internal_scaling.stop();

      return nl_status;
    }

  }; // class AlPiNeSteadyFlowSolver<...>
}
