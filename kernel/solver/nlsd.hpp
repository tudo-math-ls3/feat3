#pragma once
#ifndef FEAT_SOLVER_NLSD
#define FEAT_SOLVER_NLSD 1
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlopt_precond.hpp>

#include <deque>
namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Nonlinear Steepest Descent method for finding a minimum of an operator's gradient
     *
     * \tparam Operator_
     * Nonlinear Operator to minimise the gradient of
     *
     * \tparam Filter_
     * Filter to apply to the operator's gradient
     *
     * See \cite NW06 for an overview of optimisation techniques.
     *
     */
    template<typename Operator_, typename Filter_>
    class NLSD : public PreconditionedIterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        /// The nonlinear operator type
        typedef Operator_ OperatorType;
        /// The filter type
        typedef Filter_ FilterType;
        /// The baseclass for all applicable linesearches
        typedef Linesearch<OperatorType, FilterType> LinesearchType;

        /// Type of the operator's gradient has
        typedef typename Operator_::GradientType GradientType;
        /// Input type for the gradient
        typedef typename Operator_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Operator_::DataType DataType;

        /// Our baseclass
        typedef PreconditionedIterativeSolver<VectorType> BaseClass;
        /// Generic preconditioner
        typedef NLOptPrecond<typename Operator_::VectorTypeL, Filter_> PrecondType;

      protected:
        /// Our nonlinear operator
        Operator_& _op;
        /// The filter we apply to the gradient
        Filter_& _filter;
        /// The linesearch used along the descent direction
        std::shared_ptr<LinesearchType> _linesearch;
        /// This will be the preconditioner, or a nullptr. We need to save it ourselves because we cannot access the
        /// prepare() routine through the SolverBase pointer in our BaseClass
        std::shared_ptr<PrecondType> _precond;

        /// defect vector
        VectorType _vec_def;
        /// descend direction vector
        VectorType _vec_dir;

        /// Tolerance for function improvement
        DataType _tol_fval;
        /// Tolerance for the length of the update step
        DataType _tol_step;

        /// Current functional value
        DataType _fval;
        /// Functional value from the previous iteration
        DataType _fval_prev;

      public:
        /// For debugging purposes, all iterates can be logged to here
        std::deque<VectorType>* iterates;

      public:
        /**
         * \brief Standard constructor
         *
         * \param[in, out] op_
         * The (nonlinear) operator. Cannot be const because it saves its own state
         *
         * \param[in] filter_
         * Filter to apply to the operator's gradient
         *
         * \param[in, out] linesearch_
         * The linesearch to be used, cannot be const as internal data changes
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         * \param[in, out] precond
         * Preconditioner, defaults to nullptr. Cannot be const as internal data changes
         *
         */
        explicit NLSD(Operator_& op_, Filter_& filter_, std::shared_ptr<LinesearchType> linesearch_,
        bool keep_iterates = false, std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass("NLSD", precond),
          _op(op_),
          _filter(filter_),
          _linesearch(linesearch_),
          _precond(precond),
          _tol_fval(DataType(0)),
          _tol_step(Math::sqrt(Math::eps<DataType>())),
          iterates(nullptr)
          {
            XASSERT(_linesearch != nullptr);
            if(keep_iterates)
              iterates = new std::deque<VectorType>;
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~NLSD()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();
          // create three temporary vectors
          _vec_def = this->_op.create_vector_r();
          _vec_dir = this->_op.create_vector_r();
          _linesearch->init_symbolic();
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          if(iterates != nullptr)
            iterates->clear();

          //this->_vec_tmp.clear();
          this->_vec_dir.clear();
          this->_vec_def.clear();
          _linesearch->done_symbolic();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "NLSD-"+_linesearch->name();
        }

        /**
         * \brief Sets the tolerance for function value improvement
         *
         * \param[in] tol_fval
         * New tolerance for function value improvement.
         *
         * The convergence check is against the maximum of the absolute and relative function value.
         *
         */
        void set_tol_fval(DataType tol_fval)
        {
          _tol_fval = tol_fval;
        }

        /**
         * \brief Sets the tolerance for the linesearch step size.
         *
         * \param[in] tol_step
         * New tolerance for the linesearch step size.
         *
         * If the linesearch fails to find a new iterate because its relative update is too small, the direction
         * update will fail to produce a new search direction so the NLSD has to be terminated.
         *
         */
        void set_tol_step(DataType tol_step)
        {
          _tol_step = tol_step;
        }

        /// \copydoc BaseClass::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // save defect
          this->_vec_def.copy(vec_def);
          //this->_system_filter.filter_def(this->_vec_def);

          // clear solution vector
          vec_cor.format();

          this->_op.prepare(vec_cor, this->_filter);

          if(this->_precond != nullptr)
            this->_precond->prepare(vec_cor, this->_filter);


          // apply
          return _apply_intern(vec_cor);
        }

        /// \copydoc BaseClass::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& DOXY(vec_rhs)) override
        {
          this->_op.prepare(vec_sol, this->_filter);
          // compute defect
          this->_op.compute_grad(this->_vec_def);
          this->_vec_def.scale(this->_vec_def,DataType(-1));
          this->_filter.filter_def(this->_vec_def);

          if(this->_precond != nullptr)
            this->_precond->prepare(vec_sol, this->_filter);

          // apply
          Status st =_apply_intern(vec_sol);

          return st;
        }

      protected:
        /**
         * \brief Internal function, applies the solver
         *
         * \param[in, out] vec_sol
         * The initial guess, gets overwritten by the solution
         *
         * \returns
         * A solver status code.
         *
         * This does not have a right hand side because that is contained in the gradient of the operator and we
         * always seek grad operator(vec_sol) = 0
         *
         */
        virtual Status _apply_intern(VectorType& vec_sol)
        {
          // Reset member variables in the LineSearch
          _linesearch->reset();

          if(iterates != nullptr)
            iterates->push_back(std::move(vec_sol.clone()));

          // compute initial defect
          Status status = this->_set_initial_defect(this->_vec_def, vec_sol);
          if(status != Status::progress)
            return status;

          this->_fval = this->_op.compute_func();
          // The first direction has to be the steepest descent direction
          this->_vec_dir.clone(this->_vec_def);

          // apply preconditioner to defect vector
          //if(!this->_apply_precond(this->_vec_tmp, this->_vec_def, this->_filter))
          if(!this->_apply_precond(this->_vec_dir, this->_vec_def, this->_filter))
            return Status::aborted;

          // Compute initial eta = <d, r>
          DataType eta(this->_vec_dir.dot(this->_vec_def));

          // If the preconditioned search direction is not a descent direction, reset it to steepest descent
          // TODO: Correct the output of the preconditioner if it turns out to not have been positive definite
          if(eta <= DataType(0))
          {
            this->_vec_dir.clone(this->_vec_def);
          }

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            _fval_prev = _fval;

            // Copy information to the linesearch
            _linesearch->set_initial_fval(this->_fval);
            _linesearch->set_grad_from_defect(this->_vec_def);

            // Call the linesearch to update vec_sol
            status = _linesearch->correct(vec_sol, this->_vec_dir);

            // Copy back information from the linesearch
            this->_fval = _linesearch->get_final_fval();
            _linesearch->get_defect_from_grad(this->_vec_def);

            // Log iterates if necessary
            if(iterates != nullptr)
            {
              iterates->push_back(vec_sol.clone());
            }

            // Compute defect norm. This also performs the convergence/divergence checks.
            status = this->_set_new_defect(this->_vec_def, vec_sol);

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

            // Re-assemble preconditioner if necessary
            if(this->_precond != nullptr)
              this->_precond->prepare(vec_sol, this->_filter);

            // apply preconditioner
            if(!this->_apply_precond(_vec_dir, _vec_def, this->_filter))
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return Status::aborted;
            }

            // Compute new eta
            eta = this->_vec_dir.dot(this->_vec_def);

            // If the preconditioned search direction is not a descent direction, reset it to steepest descent
            // TODO: Correct the output of the preconditioner if it turns out to not have been positive definite
            if(eta <= DataType(0))
            {
              this->_vec_dir.clone(this->_vec_def);
            }
          }

          // We should never come to this point
          return Status::undefined;
        }

        /**
         * \brief Internal function: sets the new defect norm
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
         * A solver status code.
         */
        virtual Status _set_new_defect(const VectorType& vec_def, const VectorType& vec_sol) override
        {
          // increase iteration count
          ++this->_num_iter;

          // first, let's see if we have to compute the defect at all
          bool calc_def = false;
          calc_def = calc_def || (this->_min_iter < this->_max_iter);
          calc_def = calc_def || this->_plot;
          calc_def = calc_def || (this->_min_stag_iter > Index(0));

          // compute new defect
          if(calc_def)
            this->_def_cur = this->_calc_def_norm(vec_def, vec_sol);

          Statistics::add_solver_defect(this->_branch, double(this->_def_cur));

          // plot?
          if(this->_plot)
          {
            std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
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

          // Check for convergence of the gradient norm
          if(this->is_converged())
            return Status::success;

          // Check for convergence wrt. the function value improvement if _tol_fval says so
          if(_tol_fval > DataType(0))
          {
            // This is the factor for the relative funciton value
            DataType scale(Math::max(_fval, _fval_prev));
            // Make sure it is at least 1
            scale = Math::max(scale, DataType(1));
            // Check for success
            if(Math::abs(_fval - _fval_prev)/scale < _tol_fval)
              return Status::success;
          }

          // If the linesearch failed to make progress, the new iterate is too close to the old iterate to compute
          // a new search direction etc. so we have to abort.
          if(_linesearch->get_rel_update() < this->_tol_step)
          {
            return Status::stagnated;
          }

          // If there were too many stagnated iterations, the solver is stagnated
          if(this->_min_stag_iter > 0 && this->_num_stag_iter > this->_min_stag_iter)
            return Status::stagnated;

          // maximum number of iterations performed?
          if(this->_num_iter >= this->_max_iter)
            return Status::max_iter;

          // continue iterating
          return Status::progress;
        }

    }; // class NLSD

    /**
     * \brief Creates a new NLSD solver object
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \param[in] linesearch
     * The linesearch to use.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new NLSD object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Operator_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLSD<Operator_, Filter_>> new_nlsd(
      Operator_& op, Filter_& filter, Linesearch_& linesearch, bool keep_iterates = false)
      {
        return std::make_shared<NLSD<Operator_, Filter_>>(op, filter, linesearch,
        keep_iterates, nullptr);
      }
    template<typename Operator_, typename Filter_, typename Linesearch_, typename Precond_>
    inline std::shared_ptr<NLSD<Operator_, Filter_>> new_nlsd(
      Operator_& op, Filter_& filter, Linesearch_& linesearch, bool keep_iterates,
      std::shared_ptr<Precond_> precond)
      {
        return std::make_shared<NLSD<Operator_, Filter_>>(op, filter, linesearch,
        keep_iterates, precond);
      }
#else
    template<typename Operator_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLSD<Operator_, Filter_>> new_nlsd(
      Operator_& op, Filter_& filter, Linesearch_& linesearch, bool keep_iterates = false,
      std::shared_ptr<NLOptPrecond<typename Operator_::VectorTypeL, Filter_>> precond = nullptr)
      {
        return std::make_shared<NLSD<Operator_, Filter_>>(op, filter, linesearch,
        keep_iterates, precond);
      }
#endif
  } //namespace Solver
} // namespace FEAT

#endif
