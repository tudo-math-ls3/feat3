#pragma once
#ifndef KERNEL_SOLVER_ITERATIVE_HPP
#define KERNEL_SOLVER_ITERATIVE_HPP 1

// includes, FEAST
#include <kernel/solver/base.hpp>
#include <kernel/util/statistics.hpp>

namespace FEAST
{
  namespace Solver
  {
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
      typedef Vector_ VectorType;
      typedef typename VectorType::DataType DataType;
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
      bool _plot;

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
      explicit IterativeSolver(String plot_name) :
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
        _plot(false)
      {
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
       * \brief Sets the plot mode of the solver.
       *
       * \param[in] plot
       * If set to \c true, the solver will print a convergence plot to std::cout.
       */
      void set_plot(bool plot)
      {
        _plot = plot;
      }

      /// Sets the plot name of the solver.
      void set_plot_name(const String& name)
      {
        _plot_name = name;
      }

      /// Returns the plot name of the solver.
      String get_plot_name() const
      {
        return _plot_name;
      }

      /// checks for convergence
      bool is_converged() const
      {
        return (_def_cur <= _tol_abs) && (_def_cur <= (_tol_rel * _def_init));
      }

      /// checks for divergence
      bool is_diverged() const
      {
        return (_def_cur > _div_abs) || (_def_cur > (_div_rel * _def_init));
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
        // insert special toe to signal new start of solver
        Statistics::add_solver_toe(this->_branch, double(-1));
        //insert -1 as first defect, to signalize a new starting solver iteration run
        Statistics::add_solver_defect(this->_branch, double(-1));
        Statistics::add_solver_defect(this->_branch, this->_def_init);
        this->_num_iter = Index(0);
        this->_num_stag_iter = Index(0);

        // plot?
        if(this->_plot)
        {
          std::cout << this->_plot_name
            <<  ": " << stringify(0).pad_front(_iter_digits)
            << " : " << scientify(this->_def_init) << std::endl;
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
       * \param[in] vector
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
        bool calc_def = false;
        calc_def = calc_def || (this->_min_iter < this->_max_iter);
        calc_def = calc_def || this->_plot;
        calc_def = calc_def || (this->_min_stag_iter > Index(0));

        // save previous defect
        DataType def_old = this->_def_cur;

        // compute new defect
        if(calc_def)
          this->_def_cur = this->_calc_def_norm(vec_def, vec_sol);

        Statistics::add_solver_defect(this->_branch, this->_def_cur);

        // plot?
        if(this->_plot)
        {
          std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(_iter_digits)
            << " : " << scientify(this->_def_cur)
            << " / " << scientify(this->_def_cur / this->_def_init)
            << std::endl;
        }

        // ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(this->_def_cur))
          return Status::aborted;

        // is diverged?
        if(is_diverged())
          return Status::diverged;

        // minimum number of iterations performed?
        if(this->_num_iter < this->_min_iter)
          return Status::progress;

        // is converged?
        if(is_converged())
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
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Peter Zajac
     */
    template<typename Vector_>
    class PreconditionedIterativeSolver :
      public IterativeSolver<Vector_>
    {
    public:
      typedef Vector_ VectorType;
      typedef typename VectorType::DataType DataType;
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
      explicit PreconditionedIterativeSolver(String plot_name, std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass(plot_name),
        _precond(precond)
      {
      }

    public:
      /// virtual destructor
      virtual ~PreconditionedIterativeSolver()
      {
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        if(_precond)
          _precond->init_symbolic();
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();
        if(_precond)
          _precond->init_numeric();
      }

      virtual void init_branch(String root = "") override
      {
        BaseClass::init_branch(root);
        if(_precond)
          _precond->init_branch(root + "::" + this->name());
      }

      virtual void done_numeric() override
      {
        if(_precond)
          _precond->done_numeric();
        BaseClass::done_numeric();
      }

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
} // namespace FEAST

#endif // KERNEL_SOLVER_ITERATIVE_HPP
