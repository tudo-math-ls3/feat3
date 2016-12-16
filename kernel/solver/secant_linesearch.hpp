#pragma once
#ifndef FEAT_KERNEL_SOLVER_SECANT_LINESEARCH
#define FEAT_KERNEL_SOLVER_SECANT_LINESEARCH 1
#include <kernel/base_header.hpp>
#include <kernel/solver/linesearch.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Secant linesearch
     *
     * \tparam Operator_
     * The (nonlinear) operator to be evaluated
     *
     * \tparam Filter_
     * The filter to be applied to the operator's gradient
     *
     * This class implements a linesearch which approximately finds
     * \f[
     *   \alpha^* = \mathrm{argmin} \nabla f(x + \alpha d) \cdot d
     * \f]
     * for a given search direction \f$ d \f$ by approximating the second order derivatives along \f$ d \f$ by a
     * secant.
     *
     */
    template<typename Operator_, typename Filter_>
    class SecantLinesearch : public Linesearch<Operator_, Filter_>
    {
      public:
        /// Filter type to be applied to the gradient of the operator
        typedef Filter_ FilterType;
        /// Input vector type for the operator's gradient
        typedef typename Operator_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Operator_::DataType DataType;
        /// Our base class
        typedef Linesearch<Operator_, Filter_> BaseClass;
        /// Default initial step length
        static constexpr DataType initial_step_default = DataType(1e-2);

      protected:
        /// Step for calculating the "other" secant point in the initial step. Crucial
        DataType _secant_step;
        /// dir^T * preconditioned defect. We want to find the minimum of the functional value along dir
        DataType _eta;

      public:
        /**
         * \brief Standard constructor
         *
         * \param[in, out] op_
         * The (nonlinear) operator. Cannot be const because it saves its own state.
         *
         * \param[in] filter_
         * Filter to apply to the operator's gradient.
         *
         * \param[in] initial_step_
         * Step length for setting the "other" secant point in the first iteration. Crucial.
         *
         * \param[in] keep_iterates_
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit SecantLinesearch(
          Operator_& op_, Filter_& filter_,
          const DataType initial_step_ = initial_step_default,
          const bool keep_iterates_ = false) :
          BaseClass("S-LS", op_, filter_, keep_iterates_),
          _secant_step(initial_step_),
          _eta(DataType(0))
          {
          }

        explicit SecantLinesearch(const String& section_name, PropertyMap* section,
        Operator_& op_, Filter_& filter_) :
          BaseClass("S-LS", section_name, section, op_, filter_),
          _secant_step(initial_step_default),
          _eta(DataType(0))
          {
            auto secant_step_p = section->query("secant_step");
            if(secant_step_p.second)
            {
              set_secant_step(DataType(std::stod(secant_step_p.first)));
            }
          }

        /// \copydoc ~BaseClass()
        virtual ~SecantLinesearch()
        {
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "SecantLinesearch";
        }

        /**
         * \brief Sets the length of the first secant step
         *
         */
        void set_secant_step(DataType secant_step)
        {
          XASSERT(secant_step > DataType(0));

          _secant_step = secant_step;
        }

        /**
         * \brief Applies the solver, setting the initial guess to zero.
         *
         * \param[out] vec_cor
         * Solution, gets zeroed
         *
         * \param[in] vec_dir
         * Search direction
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_dir) override
        {
          // clear solution vector
          vec_cor.format();
          this->_op.prepare(vec_cor, this->_filter);

          // apply
          return _apply_intern(vec_cor, vec_dir);
        }

        /**
         * \brief Applies the solver, making use of an initial guess
         *
         * \param[out] vec_sol
         * Initial guess, gets overwritten by the solution
         *
         * \param[in] vec_dir
         * Search direction
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        virtual Status correct(VectorType& vec_sol, const VectorType& vec_dir) override
        {
          this->_op.prepare(vec_sol, this->_filter);
          // apply
          Status st =_apply_intern(vec_sol, vec_dir);

          return st;
        }

      protected:
        /// \copydoc NewtonRaphsonLinesearch::_apply_intern()
        virtual Status _apply_intern(VectorType& sol, const VectorType& dir)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

          // compute initial defect
          Status status(Status::progress);
          this->_num_iter = Index(0);

          this->_vec_initial_sol.copy(sol);

          // Norm of the search direction
          this->_norm_dir = dir.norm2();
          // Norm of the initial guess
          this->_norm_sol = sol.norm2();

          this->_fval_min = this->_fval_0;
          this->_alpha_min = DataType(0);

          // The second secant point in the first iteration is x + _secant_step * dir
          DataType alpha(_secant_step/this->_norm_dir);
          DataType alpha_hidate(alpha);

          _eta = dir.dot(this->_vec_grad);
          if(_eta > DataType(0))
            throw InternalError(__func__,__FILE__,__LINE__,"Search direction is not a descent direction: "
                +stringify_fp_sci(_eta));

          // The first "other" point for the secant
          sol.axpy(dir, this->_vec_initial_sol, _secant_step/this->_norm_dir);

          _eta = dir.dot(this->_vec_grad);
          this->_def_init = Math::abs(_eta);

          // start iterating
          while(status == Status::progress)
          {
            IterationStats stat(*this);

            // Increase iteration count
            ++this->_num_iter;

            DataType fval(0);
            this->_op.prepare(sol, this->_filter);
            this->_op.eval_fval_grad(fval, this->_vec_grad);
            this->_filter.filter_def(this->_vec_grad);

            if(fval < this->_fval_min)
            {
              this->_fval_min = fval;
              this->_alpha_min = alpha;
            }

            // Set new defect, do convergence checks and compute the new _alpha
            DataType eta_prev = _eta;

            // Compute eta and alpha
            _eta = this->_vec_grad.dot(dir);
            this->_def_cur = Math::abs(_eta);

            // ensure that the defect is neither NaN nor infinity
            if(!Math::isfinite(this->_def_cur))
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            // is diverged?
            if(this->is_diverged())
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::diverged, this->get_num_iter()));
              return Status::diverged;
            }

            // Check if the diffence in etas is too small, thus leading to a huge update of relative size > sqrt(eps)
            if(Math::abs(_eta - eta_prev) < Math::abs(alpha*_eta)*this->_norm_dir*Math::sqrt(Math::eps<DataType>()))
            {
              // If we are not successful, update the solution with the best step found so far
              sol.axpy(dir, this->_vec_initial_sol, this->_alpha_min);
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::stagnated, this->get_num_iter()));
              return Status::stagnated;
            }

            // Update alpha according to secant formula
            alpha_hidate *= _eta/(eta_prev - _eta);
            alpha += alpha_hidate;

            // Update the solution
            sol.axpy(dir, this->_vec_initial_sol, alpha);

            //Statistics::add_solver_defect(this->_branch, double(this->_def_cur));

            // plot?
            if(this->_plot)
            {
              std::cout << this->_plot_name
              <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
              << " : " << stringify_fp_sci(this->_def_cur)
              << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
              << std::endl;
            }

            // minimum number of iterations performed?
            if(this->_num_iter < this->_min_iter)
              continue;

            // is converged?
            if(this->is_converged())
              status = Status::success;

            // maximum number of iterations performed?
            if(this->_num_iter >= this->_max_iter)
              status = Status::max_iter;

            //// check for stagnation?
            //if(this->_min_stag_iter > Index(0))
            //{
            //  // did this iteration stagnate?
            //  if(this->_def_cur >= this->_stag_rate * def_old)
            //  {
            //    // increment stagnation count
            //    ++this->_num_stag_iter;

            //    if(this->_num_stag_iter >= this->_min_stag_iter)
            //      return Status::stagnated;
            //  }
            //  else
            //  {
            //    // this iteration did not stagnate
            //    this->_num_stag_iter = Index(0);
            //  }
            //}

            // If we are not successful, update the solution with the best step found so far
            if(status != Status::progress)
            {
              sol.axpy(dir, this->_vec_initial_sol, this->_alpha_min);
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
              return status;
            }

          }
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return Status::undefined;
        }

    }; // class SecantLinesearch

    /**
     * \brief Creates a new SecantLinesearch object
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] initial_step
     * Length for first secant step.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new SecantLinesearch object.
     */
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<SecantLinesearch<Operator_, Filter_>> new_secant_linesearch(
      Operator_& op, Filter_& filter,
      typename Operator_::DataType initial_step = SecantLinesearch<Operator_, Filter_>::initial_step_default,
      bool keep_iterates = false)
      {
        return std::make_shared<SecantLinesearch<Operator_, Filter_>>(op, filter, initial_step, keep_iterates);
      }

    /**
     * \brief Creates a new SecantLinesearch object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new SecantLinesearch object.
     */
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<SecantLinesearch<Operator_, Filter_>> new_secant_linesearch(
      const String& section_name, PropertyMap* section,
      Operator_& op, Filter_& filter)
      {
        return std::make_shared<SecantLinesearch<Operator_, Filter_>>(section_name, section, op, filter);
      }

  } // namespace Solver
} // namespace FEAT

#endif // FEAT_KERNEL_SOLVER_SECANT_LINESEARCH
