#pragma once
#ifndef KERNEL_SOLVER_ITERATIVE_HPP
#define KERNEL_SOLVER_ITERATIVE_HPP 1

// includes, FEAT
#include <kernel/solver/base.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/statistics.hpp>

namespace FEAT
{
  namespace Solver
  {

    /**
     * \brief Solver plot modes enumeration
     *
     */
    enum class PlotMode
    {
      /// No plotting whatsoever
      none = 0,
      /// Plot every iteration (if applicable)
      iter,
      /// Plot a summary after each solver run
      summary,
      /// Plot every iteration (if applicable) and a summary
      all
    };

    /// \cond internal
    inline std::ostream& operator<<(std::ostream& os, PlotMode mode)
    {
      switch(mode)
      {
        case PlotMode::none:
          return os << "none";
        case PlotMode::iter:
          return os << "iter";
        case PlotMode::summary:
          return os << "summary";
        case PlotMode::all:
          return os << "all";
        default:
          return os << "-unknown-";
      }
    }

    inline void operator<<(PlotMode& mode, const String& mode_name)
    {
        if(mode_name == "none")
          mode = PlotMode::none;
        else if(mode_name == "iter")
          mode = PlotMode::iter;
        else if(mode_name == "summary")
          mode = PlotMode::summary;
        else if(mode_name == "all")
          mode = PlotMode::all;
        else
          throw InternalError(__func__, __FILE__, __LINE__, "Unknown PlotMode identifier string "
              +mode_name);
    }
    /// \endcond

    /**
     * \brief Abstract base-class for iterative solvers.
     *
     * This class template acts as an abstract base class for iterative solvers.
     * It also implements various auxiliary features for convergence control.
     *
     * \tparam Vector_
     * The class of the vector that is passed to the solver in the \c solve() method.
     *
     * \author Peter Zajac
     */
    template<typename Vector_>
    class IterativeSolver :
      public SolverBase<Vector_>
    {
    public:
      /// The vector type this solver can be applied to
      typedef Vector_ VectorType;
      /// Floating point type
      typedef typename VectorType::DataType DataType;
      /// The base class
      typedef SolverBase<VectorType> BaseClass;

    protected:
      /// name of the solver in plots
      String _plot_name;
      /// relative tolerance parameter
      DataType _tol_rel;
      /// relative tolerance parameter
      DataType _tol_abs;
      /// relative divergence parameter
      DataType _div_rel;
      /// absolute divergence parameter
      DataType _div_abs;
      /// stagnation rate
      DataType _stag_rate;
      /// minimum number of iterations
      Index _min_iter;
      /// maximum number of iterations
      Index _max_iter;
      /// number of performed iterations
      Index _num_iter;
      /// minimum number of stagnation iterations
      Index _min_stag_iter;
      /// number of consecutive stagnated iterations
      Index _num_stag_iter;
      /// initial defect
      DataType _def_init;
      /// current defect
      DataType _def_cur;
      /// iteration count digits for plotting
      Index _iter_digits;
      /// whether to plot something
      PlotMode _plot_mode;
      /// whether to skip defect computation if possible
      bool _skip_def_calc;

      /**
       * \brief Protected constructor
       *
       * This constructor initialises the following values:
       *
       * - relative tolerance: sqrt(eps) (~1E-8 for double)
       * - absolute tolerance: 1/eps^2 (~1E+32 for double)
       * - relative divergence: 1/eps (~1E+16 for double)
       * - absolute divergence: 1/eps^2 (~1E+32 for double)
       * - stagnation rate: 0.95
       * - minimum iterations: 0
       * - maximum iterations: 100
       * - minimum stagnation iterations: 0
       * - convergence plot: false
       *
       * \param[in] plot_name
       * Specifies the name of the iterative solver. This is used as a prefix for the convergence plot.
       */
      explicit IterativeSolver(const String& plot_name) :
        BaseClass(),
        _plot_name(plot_name),
        _tol_rel(Math::sqrt(Math::eps<DataType>())),
        _tol_abs(DataType(1) / Math::sqr(Math::eps<DataType>())),
        _div_rel(DataType(1) / Math::eps<DataType>()),
        _div_abs(DataType(1) / Math::sqr(Math::eps<DataType>())),
        _stag_rate(DataType(0.95)),
        _min_iter(0),
        _max_iter(100),
        _num_iter(0),
        _min_stag_iter(0),
        _num_stag_iter(0),
        _def_init(0),
        _def_cur(0),
        _iter_digits(Math::ilog10(_max_iter)),
        _plot_mode(PlotMode::none),
#ifdef FEAT_HAVE_MPI
        _skip_def_calc(false) // not allowed as this may cause deadlocks
#else
        _skip_def_calc(true) // no potential problem in non-MPI builds
#endif
      {
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] plot_name
       * The name of the solver.
       *
       */
      explicit IterativeSolver(const String& plot_name, const String& section_name, PropertyMap* section) :
        BaseClass(section_name, section),
        _plot_name(plot_name),
        _tol_rel(Math::sqrt(Math::eps<DataType>())),
        _tol_abs(DataType(1) / Math::sqr(Math::eps<DataType>())),
        _div_rel(DataType(1) / Math::eps<DataType>()),
        _div_abs(DataType(1) / Math::sqr(Math::eps<DataType>())),
        _stag_rate(DataType(0.95)),
        _min_iter(0),
        _max_iter(100),
        _num_iter(0),
        _min_stag_iter(0),
        _num_stag_iter(0),
        _def_init(0),
        _def_cur(0),
        _iter_digits(Math::ilog10(_max_iter)),
        _plot_mode(PlotMode::none),
#ifdef FEAT_HAVE_MPI
        _skip_def_calc(false) // not allowed as this may cause deadlocks
#else
        _skip_def_calc(true) // no potential problem in non-MPI builds
#endif
      {
        Dist::Comm comm(Dist::Comm::world());

        auto plot_mode_p = section->get_entry("plot_mode");
        if (plot_mode_p.second)
        {
          PlotMode plot_mode;
          plot_mode << plot_mode_p.first;
          if(comm.rank() == 0)
          {
            set_plot_mode(plot_mode);
          }
        }

        auto plot_name_p = section->get_entry("plot_name");
        if (plot_name_p.second)
        {
          if(comm.rank() == 0)
          {
            set_plot_name(plot_name_p.first);
          }
        }

        auto tol_abs_p = section->get_entry("tol_abs");
        if (tol_abs_p.second)
          set_tol_abs(DataType(std::stod(tol_abs_p.first)));

        auto tol_rel_p = section->get_entry("tol_rel");
        if (tol_rel_p.second)
          set_tol_rel(DataType(std::stod(tol_rel_p.first)));

        auto div_abs_p = section->get_entry("div_abs");
        if (div_abs_p.second)
          set_div_abs(DataType(std::stod(div_abs_p.first)));

        auto div_rel_p = section->get_entry("div_rel");
        if (div_rel_p.second)
          set_div_rel(DataType(std::stod(div_rel_p.first)));

        auto stag_rate_p = section->get_entry("stag_rate");
        if (stag_rate_p.second)
          set_stag_rate(DataType(std::stod(stag_rate_p.first)));

        auto max_iter_p = section->get_entry("max_iter");
        if (max_iter_p.second)
          set_max_iter(Index(std::stoul(max_iter_p.first)));

        auto min_iter_p = section->get_entry("min_iter");
        if (min_iter_p.second)
          set_min_iter(Index(std::stoul(min_iter_p.first)));

        auto min_stag_iter_p = section->get_entry("min_stag_iter");
        if (min_stag_iter_p.second)
          set_min_stag_iter(Index(std::stoul(min_stag_iter_p.first)));

      }

    public:
      /// Sets the relative tolerance for the solver.
      void set_tol_rel(DataType tol_rel)
      {
        _tol_rel = tol_rel;
      }

      /// Sets the absolute tolerance for the solver.
      void set_tol_abs(DataType tol_abs)
      {
        _tol_abs = tol_abs;
      }

      /// Returns the relative tolerance.
      DataType get_tol_rel() const
      {
        return _tol_rel;
      }

      /// Returns the absolute tolerance.
      DataType get_tol_abs() const
      {
        return _tol_abs;
      }

      /// Sets the relative divergence for the solver.
      void set_div_rel(DataType div_rel)
      {
        _div_rel = div_rel;
      }

      /// Sets the absolute divergence for the solver.
      void set_div_abs(DataType div_abs)
      {
        _div_abs = div_abs;
      }

      /// Returns the relative divergence.
      DataType get_div_rel() const
      {
        return _div_rel;
      }

      /// Returns the absolute divergence.
      DataType get_div_abs() const
      {
        return _div_abs;
      }

      /// Sets the stagnation rate fot the solver.
      void set_stag_rate(DataType rate)
      {
        _stag_rate = rate;
      }

      /// Returns the stagnation rate
      DataType get_stag_rate() const
      {
        return _stag_rate;
      }

      /// Sets the minimum stagnate iteration count for the solver
      void set_min_stag_iter(Index min_iter)
      {
        _min_stag_iter = min_iter;
      }

      /// Returns the minimum stagnation iteration count.
      Index get_min_stag_iter() const
      {
        return _min_stag_iter;
      }

      /// Sets the minimum iteration count for the solver.
      void set_min_iter(Index min_iter)
      {
        _min_iter = min_iter;
      }

      /// Sets the maximum iteration count for the solver.
      void set_max_iter(Index max_iter)
      {
        _max_iter = max_iter;
        _iter_digits = Math::ilog10(_max_iter);
      }

      /// Returns number of performed iterations
      Index get_num_iter() const
      {
        return _num_iter;
      }

      /// Returns the minimal number of iterations
      Index get_min_iter() const
      {
        return _min_iter;
      }

      /// Returns the maximum number of iterations
      Index get_max_iter() const
      {
        return _max_iter;
      }

      /**
       * \brief Specifies whether defect calculation is allowed to be skipped.
       *
       * \warning
       * Skipping defect calculation can lead to deadlocks in parallel solver implementations
       * if one (but not all) processes need to compute the defect for some reason (like plotting) !\n
       * Use this function only when you really know what you are doing!
       *
       * \param[in] skip
       * Specifies whether defect computation is allowed to be skipped.
       */
      void skip_defect_calc(bool skip)
      {
        _skip_def_calc = skip;
      }

      /**
       * \brief Sets the plot mode of the solver.
       *
       * \param[in] plot_mode
       * If set to anything but PlotMode::none, the solver will print information to std::cout
       */
      void set_plot_mode(const PlotMode plot_mode)
      {
        _plot_mode = plot_mode;
      }

      /// Sets the plot name of the solver.
      void set_plot_name(const String& plot_name)
      {
        _plot_name = plot_name;
      }

      /// Returns the plot name of the solver.
      String get_plot_name() const
      {
        return _plot_name;
      }

      /// checks for convergence
      bool is_converged() const
      {
        return is_converged(_def_cur);
      }

      /// checks for convergence
      bool is_converged(const DataType def_cur) const
      {
        return (def_cur <= _tol_abs) && (def_cur <= (_tol_rel * _def_init));
      }

      /// checks for divergence
      bool is_diverged() const
      {
        return is_diverged(_def_cur);
      }

      /// checks for divergence
      bool is_diverged(const DataType def_cur) const
      {
        return (def_cur > _div_abs) || (def_cur > (_div_rel * _def_init));
      }

      /// Returns the initial defect
      DataType get_def_initial() const
      {
        return _def_init;
      }

      /// Returns the final defect
      DataType get_def_final() const
      {
        return _def_cur;
      }

      /// Returns the overall convergence rate.
      DataType get_conv_rate() const
      {
        // no iterations performed?
        if(_num_iter <= Index(0))
          return DataType(0);
        // initial defect zero?
        if(_def_init < Math::eps<DataType>())
          return DataType(0);

        // compute convergence rate: (def_final / def_initial) ^ (1 / #iter)
        return Math::pow(_def_cur / _def_init, DataType(1) / DataType(_num_iter));
      }

      /// \copydoc SolverBase::write_config()
      virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name = "") const override
      {
        XASSERT(parent != nullptr);

        PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

        my_section->add_entry("plot_mode", stringify(_plot_mode));
        my_section->add_entry("tol_rel", stringify_fp_sci(_tol_rel));
        my_section->add_entry("tol_abs", stringify_fp_sci(_tol_abs));
        my_section->add_entry("div_rel", stringify_fp_sci(_div_rel));
        my_section->add_entry("div_abs", stringify_fp_sci(_div_abs));
        my_section->add_entry("stag_rate", stringify_fp_sci(_stag_rate));
        my_section->add_entry("max_iter", stringify(_max_iter));
        my_section->add_entry("min_iter", stringify(_min_iter));
        my_section->add_entry("min_stag_iter", stringify(_min_stag_iter));

        return my_section;

      }

      /**
       * \brief Solver correction method
       *
       * This method applies the solver represented by this object onto a given right-hand-side vector
       * and updates the corresponding solution vector.
       *
       * In contrast to the apply() method of the SolverBase base class, this method uses the
       * vector \p vec_sol as the initial solution vector for the iterative solution process instead of
       * ignoring its contents upon entry and starting with the null vector.
       *
       * \param[in,out] vec_sol
       * The vector that contains the initial solution upon entry and receives the solution
       * of the linear system upon exit.
       *
       * \param[in] vec_rhs
       * The vector that represents the right-hand-side of the linear system to be solved.
       *
       * \returns
       * A Status code that represents the status of the solution step.
       */
      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) = 0;

    protected:
      /**
       * \brief Computes the defect norm.
       *
       * \param[in] vec_def
       * The current defect vector.
       *
       * \param[in] vec_sol
       * The current solution vector approximation.
       */
      virtual DataType _calc_def_norm(const VectorType& vec_def, const VectorType& DOXY(vec_sol))
      {
        return vec_def.norm2();
      }

      /**
       * \brief Plot a summary of the last solver run
       *
       */
      virtual void plot_summary(const Status st) const
      {
        // Print solver summary
        if(!_plot_summary())
          return;

        Dist::Comm comm_world(Dist::Comm::world());

        String msg(this->get_plot_name()+ ": its: "+stringify(this->get_num_iter())+" ("+ stringify(st)+")\n");
        msg += this->get_plot_name()  +": defect norm: "+stringify_fp_sci(this->_def_init)
          + " -> "+stringify_fp_sci(this->_def_cur)
          + ", factor " +stringify_fp_sci(this->_def_cur/this->_def_init);

        comm_world.print(msg);
      }

      /**
       * \brief Plot every iteration?
       *
       * \returns \c true if the plot mode is set to \c iter or \c all.
       */
      bool _plot_iter() const
      {
        return _plot_mode == PlotMode::iter || _plot_mode == PlotMode::all;
      }

      /**
       * \brief Plot summary?
       *
       * \returns \c true if the plot mode is set to \c summary or \c all.
       */
      bool _plot_summary() const
      {
        return _plot_mode == PlotMode::summary || _plot_mode == PlotMode::all;
      }

      /**
       * \brief Internal function: sets the initial defect vector
       *
       * \param[in] vec_def
       * The initial defect vector.
       *
       * \param[in] vec_sol
       * The current solution vector approximation.
       *
       * \returns
       * A Status code.
       */
      virtual Status _set_initial_defect(const VectorType& vec_def, const VectorType& vec_sol)
      {
        // store new defect
        this->_def_init = this->_def_cur = this->_calc_def_norm(vec_def, vec_sol);
        this->_num_iter = Index(0);
        this->_num_stag_iter = Index(0);
        Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_init, this->get_num_iter()));

        // plot?
        if(this->_plot_mode == PlotMode::iter || this->_plot_mode == PlotMode::all)
        {
          std::cout << this->_plot_name
            <<  ": " << stringify(0).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_init) << std::endl;
        }

        // ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(this->_def_init))
          return Status::aborted;

        // check if the initial defect is zero; we test against eps^2 here
        if(this->_def_init <= Math::sqr(Math::eps<DataType>()))
          return Status::success;

        // continue iterating
        return Status::progress;
      }

      /**
       * \brief Internal function: sets the new (next) defect vector
       *
       * This function computes the defect vector's norm, increments the iteration count,
       * plots an output line to std::cout and checks whether any of the stopping criterions is fulfilled.
       *
       * \param[in] vec_def
       * The new defect vector.
       *
       * \param[in] vec_sol
       * The current solution vector approximation.
       *
       * \returns
       * A Status code.
       */
      virtual Status _set_new_defect(const VectorType& vec_def, const VectorType& vec_sol)
      {
        // increase iteration count
        ++this->_num_iter;

        // first, let's see if we have to compute the defect at all
        bool calc_def = !_skip_def_calc;
        calc_def = calc_def || (this->_min_iter < this->_max_iter);
        calc_def = calc_def || (this->_plot_iter());
        calc_def = calc_def || (this->_min_stag_iter > Index(0));

        // save previous defect
        const DataType def_old = this->_def_cur;

        // compute new defect
        if(calc_def)
        {
          this->_def_cur = this->_calc_def_norm(vec_def, vec_sol);
          Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_cur, this->get_num_iter()));
        }

        // plot?
        if(this->_plot_iter())
        {
          std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << " / " << stringify_fp_fix(this->_def_cur / def_old)
            << std::endl;
        }

        // ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(this->_def_cur))
          return Status::aborted;

        // is diverged?
        if(this->is_diverged())
          return Status::diverged;

        // minimum number of iterations performed?
        if(this->_num_iter < this->_min_iter)
          return Status::progress;

        // is converged?
        if(this->is_converged())
          return Status::success;

        // maximum number of iterations performed?
        if(this->_num_iter >= this->_max_iter)
          return Status::max_iter;

        // check for stagnation?
        if(this->_min_stag_iter > Index(0))
        {
          // did this iteration stagnate?
          if(this->_def_cur >= this->_stag_rate * def_old)
          {
            // increment stagnation count
            if(++this->_num_stag_iter >= this->_min_stag_iter)
              return Status::stagnated;
          }
          else
          {
            // this iteration did not stagnate
            this->_num_stag_iter = Index(0);
          }
        }

        // continue iterating
        return Status::progress;
      }

      /**
       * \brief Internal function: sets the new (next) defect norm
       *
       * This function takes a precalculated defect vector's norm, increments the iteration count,
       * plots an output line to std::cout and checks whether any of the stopping criterions is fulfilled.
       *
       * \param[in] def_cur_norm
       * The new defect norm.
       *
       * \returns
       * A Status code.
       *
       * \note This function is preferred over _set_new_defect when using asynchronous mpi operations.
       */
      virtual Status _update_defect(const DataType def_cur_norm)
      {
        // increase iteration count
        ++this->_num_iter;

        // save previous defect
        const DataType def_old = this->_def_cur;

        this->_def_cur = def_cur_norm;

        Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_cur, this->get_num_iter()));

        // plot?
        if(this->_plot_iter())
        {
          std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << " / " << stringify_fp_fix(this->_def_cur / def_old)
            << std::endl;
        }

        // ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(this->_def_cur))
          return Status::aborted;

        // is diverged?
        if(this->is_diverged())
          return Status::diverged;

        // minimum number of iterations performed?
        if(this->_num_iter < this->_min_iter)
          return Status::progress;

        // is converged?
        if(this->is_converged())
          return Status::success;

        // maximum number of iterations performed?
        if(this->_num_iter >= this->_max_iter)
          return Status::max_iter;

        // check for stagnation?
        if(this->_min_stag_iter > Index(0))
        {
          // did this iteration stagnate?
          if(this->_def_cur >= this->_stag_rate * def_old)
          {
            // increment stagnation count
            if(++this->_num_stag_iter >= this->_min_stag_iter)
              return Status::stagnated;
          }
          else
          {
            // this iteration did not stagnate
            this->_num_stag_iter = Index(0);
          }
        }

        // continue iterating
        return Status::progress;
      }
    }; // class IterativeSolver

    /**
     * \brief Override of solve() for IterativeSolver solvers
     */
    template<
      typename Vector_,
      typename Matrix_,
      typename Filter_>
      inline Status solve(
        IterativeSolver<Vector_>& solver,
        Vector_& vec_sol,
        const Vector_& vec_rhs,
        const Matrix_& DOXY(matrix),
        const Filter_& DOXY(filter))
    {
      // simply call the 'correct' method
      return solver.correct(vec_sol, vec_rhs);
    }

    /**
     * \brief Abstract base-class for preconditioned iterative solvers.
     *
     * This class extends the functionality of the IterativeSolver class template by providing overrides
     * for the initialisation and finalisation methods of the SolverBase class template, which take
     * care of forwarding these steps to the preconditioner.
     *
     * \tparam Vector_
     * The vector type this solver can be applied to
     *
     * \author Peter Zajac
     */
    template<typename Vector_>
    class PreconditionedIterativeSolver :
      public IterativeSolver<Vector_>
    {
    public:
      /// The vector type this solver can be applied to
      typedef Vector_ VectorType;
      /// Floating point data type
      typedef typename VectorType::DataType DataType;
      /// Our base class
      typedef IterativeSolver<VectorType> BaseClass;
      /// the interface for the preconditioner
      typedef SolverBase<VectorType> PrecondType;

    protected:
      /// the pointer to the preconditioner
      std::shared_ptr<PrecondType> _precond;

      /**
       * \brief Constructor
       *
       * \param[in] plot_name
       * The name of the solver.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be nullptr.
       */
      explicit PreconditionedIterativeSolver(const String& plot_name,std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass(plot_name),
        _precond(precond)
      {
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] plot_name
       * The name of the solver.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be nullptr.
       */
      explicit PreconditionedIterativeSolver(const String& plot_name, const String& section_name,
      PropertyMap* section, std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass(plot_name, section_name, section),
        _precond(precond)
      {
      }

    public:
      /// virtual destructor
      virtual ~PreconditionedIterativeSolver()
      {
      }

      /// \copydoc SolverBase::write_config()
      virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name) const override
      {
        XASSERT(parent != nullptr);

        PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

        if(_precond == nullptr)
        {
          my_section->add_entry("precon", "none");
        }
        else
        {
          my_section->add_entry("precond", _precond->get_section_name());
          _precond->write_config(parent);
        }

        return my_section;

      }

      /// \copydoc SolverBase::init_symbolic()
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        if(_precond)
          _precond->init_symbolic();
      }

      /// \copydoc SolverBase::init_numeric()
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        if(_precond)
          _precond->init_numeric();
      }

      /// \copydoc SolverBase::done_numeric()
      virtual void done_numeric() override
      {
        if(_precond)
          _precond->done_numeric();
        BaseClass::done_numeric();
      }

      /// \copydoc SolverBase::done_symbolic()
      virtual void done_symbolic() override
      {
        if(_precond)
          _precond->done_symbolic();
        BaseClass::done_symbolic();
      }

    protected:
      /**
       * \brief Applies the preconditioner onto a defect vector.
       *
       * \note
       * If no preconditioner is present, this function will simply copy the input vector's
       * contents into the output vector, therefore emulating an "identity preconditioner".
       *
       * \param[in,out] vec_cor
       * A reference to the vector that shall receive the preconditioned defect.
       *
       * \param[in] vec_def
       * A reference to the vector that is to be preconditioned.
       *
       * \param[in] filter
       * A reference to the system filter. This filter is only used if no preconditioner is present.
       *
       * \returns
       * \c true, if the preconditioner application was successful, otherwise \c false.
       */
      template<typename Filter_>
      bool _apply_precond(VectorType& vec_cor, const VectorType& vec_def, const Filter_& filter)
      {
        if(this->_precond)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), this->_precond->name()));
          return status_success(this->_precond->apply(vec_cor, vec_def));
        }
        else
        {
          vec_cor.copy(vec_def);
          filter.filter_cor(vec_cor);
          return true;
        }
      }
    }; // class PreconditionedIterativeSolver<...>
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_ITERATIVE_HPP
