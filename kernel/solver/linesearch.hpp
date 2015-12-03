#pragma once
#ifndef FEAST_KERNEL_SOLVER_LINESEARCH
#define FEAST_KERNEL_SOLVER_LINESEARCH 1
#include <kernel/base_header.hpp>
#include <kernel/solver/iterative.hpp>

#include <deque>

namespace FEAST
{
  namespace Solver
  {
    /**
     * \brief Newton Raphson linesearch
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
     * for a given search direction \f$ d \f$ by applying a Newton Raphson iteration to this.
     *
     */
    template<typename Operator_, typename Filter_>
    class NRLinesearch : public IterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        /// Filter type to be applied to the gradient of the operator
        typedef Filter_ FilterType;
        /// Input vector type for the operator's gradient
        typedef typename Operator_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Operator_::DataType DataType;
        /// Our base class
        typedef IterativeSolver<typename Operator_::VectorTypeR> BaseClass;

      protected:
        /// The (nonlinear) operator
        // Note that this cannot be const, as the operator saves its state and thus changes
        Operator_& _op;
        /// The filter to be applied to the operator's gradient
        const Filter_& _filter;

        /// Gradient vector
        VectorType _vec_grad;
        /// temporary vector
        VectorType _vec_tmp;

        /// Line search parameter
        DataType _alpha;
        /// For keeping track of the total change in _alpha
        DataType _sum_alpha;
        /// The 2-norm of the search direction
        DataType _norm_dir;
        /// The 2-norm of the iterate
        DataType _norm_sol;

      public:
        /// For debugging purposes, it is possible to save all iterates to this
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
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit NRLinesearch(Operator_& op_, const Filter_& filter_, bool keep_iterates = false) :
          BaseClass("NR-LS"),
          _op(op_),
          _filter(filter_),
          _alpha(DataType(0)),
          _sum_alpha(DataType(0)),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;

            this->_min_stag_iter = Index(2);
            this->_min_iter = Index(3);
          }

        /// \copydoc ~BaseClass()
        virtual ~NRLinesearch()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          // Create temporary vector
          _vec_tmp = this->_op.create_vector_r();
          _vec_grad = this->_op.create_vector_r();
          BaseClass::init_symbolic();
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          // Clear temporary vector
          _vec_tmp.clear();
          _vec_grad.clear();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "Newton-Raphson-Linesearch";
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
          this->_op.prepare(vec_cor);

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
          this->_op.prepare(vec_sol);
          // apply
          Status st =_apply_intern(vec_sol, vec_dir);

          return st;
        }

        /**
         * \brief Get the relative update of the solver application
         *
         * The linesearch updates the solution according to
         * \f[
         *   x \mapsto x + \alpha d,
         * \f]
         * so the relative update is
         * \f[
         *   \hat{\alpha} = | \alpha | \frac{\| d \|_2}{ \| x \|_2}.
         * \f]
         *
         * Used for determining if the linesearch stagnated without updating the solution in a significant way.
         *
         * \returns The relative update.
         */
        DataType get_rel_update()
        {
          return Math::abs(_sum_alpha)*_norm_dir/_norm_sol;
        }

      protected:
        /**
         * \brief Internal function: Applies the solver
         *
         * \param[in, out] sol
         * Initial guess, gets overwritten by solution
         *
         * \param[in] dir
         * Search direction
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        virtual Status _apply_intern(VectorType& sol, const VectorType& dir)
        {
          Status status(Status::progress);
          this->_num_iter = Index(0);

          // Norm of the search direction vector
          _norm_dir = dir.norm2();
          // Norm of the initial guess
          _norm_sol = sol.norm2();
          // This will be the sum of all updates
          _sum_alpha = DataType(0);

          // Set state in the operator
          _op.prepare(sol);
          // Compute and filter gradient
          _op.compute_gradient(this->_vec_grad);
          _filter.filter_def(this->_vec_grad);

          // Compute initial defect. We want to minimise d^T * grad(_op)
          this->_def_init = Math::abs(dir.dot(this->_vec_grad));

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            this->_op.prepare(sol);
            _op.compute_gradient(this->_vec_grad);
            _filter.filter_def(this->_vec_grad);

            status = this->_set_new_defect(dir, this->_vec_grad);

            // Update solution: sol <- sol + _alpha*dir
            _sum_alpha += _alpha;
            sol.axpy(dir, sol, _alpha);

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

          }

          return Status::undefined;

        }

        /**
         * \brief Sets new defect norm and does convergence checks
         *
         * \param[in] dir
         * Search direction
         *
         * \param[in] grad
         * Gradient at current iterate
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        Status _set_new_defect(const VectorType& dir, const VectorType& grad) override
        {
          // Increase iteration count
          ++this->_num_iter;

          // First let's see if we have to compute the defect at all
          bool calc_def = false;
          calc_def = calc_def || (this->_min_iter < this->_max_iter);
          calc_def = calc_def || this->_plot;
          calc_def = calc_def || (this->_min_stag_iter > Index(0));

          // Update defect
          this->_def_cur = Math::abs(grad.dot(dir));

          // Compute new _alpha = - grad.dot(dir) / dir.dot(Hess*dir)
          _op.apply_hess(this->_vec_tmp, dir);
          _filter.filter_def(this->_vec_tmp);
          _alpha = -grad.dot(dir)/dir.dot(this->_vec_tmp);

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
            if(Math::sqr(_alpha)*this->_def_init <= Math::eps<DataType>())
            {
              // increment stagnation count
              if(++this->_num_stag_iter >= this->_min_stag_iter)
              {
                return Status::stagnated;
              }
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
    }; // class NRLinesearch

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
    class SecantLinesearch : public IterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        /// Filter type to be applied to the gradient of the operator
        typedef Filter_ FilterType;
        /// Input vector type for the operator's gradient
        typedef typename Operator_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Operator_::DataType DataType;
        /// Our base class
        typedef IterativeSolver<typename Operator_::VectorTypeR> BaseClass;

      protected:
        /// The (nonlinear) operator
        // Note that this cannot be const, as the operator saves its state and thus changes
        Operator_& _op;
        /// The filter to be applied to the operator's gradient
        const Filter_& _filter;

        /// temporary vector
        VectorType _vec_tmp;

        /// Line search parameter
        DataType _alpha;
        /// For keeping track of the total change in _alpha
        DataType _sum_alpha;
        /// Step for calculating the "other" secant point in the initial step. Crucial
        DataType _sigma_0;
        /// dir^T * preconditioned defect. We want to find the minimum of the along dir
        DataType _eta;
        /// The 2-norm of the search direction
        DataType _norm_dir;
        /// The 2-norm of the iterate
        DataType _norm_sol;

      public:
        /// For debugging purposes, it is possible to save all iterates to this
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
         * \param[in] initial_step
         * Step length for setting the "other" secant point in the first iteration. Crucial.
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit SecantLinesearch(Operator_& op_, const Filter_& filter_, DataType initial_step = DataType(1e-2),
        bool keep_iterates = false) :
          BaseClass("S-LS"),
          _op(op_),
          _filter(filter_),
          _alpha(DataType(0)),
          _sum_alpha(DataType(0)),
          _sigma_0(initial_step),
          _eta(DataType(0)),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;

            this->_min_stag_iter = Index(2);
            this->_min_iter = Index(5);

            this->set_plot(true);
          }

        /// \copydoc ~BaseClass()
        virtual ~SecantLinesearch()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          // Create temporary vector
          _vec_tmp = this->_op.create_vector_r();
          BaseClass::init_symbolic();
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          // Clear temporary vector
          _vec_tmp.clear();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "SecantLinesearch";
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
          this->_op.prepare(vec_cor);

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
          this->_op.prepare(vec_sol);
          // apply
          Status st =_apply_intern(vec_sol, vec_dir);

          return st;
        }

        /// \copydoc NRLinesearch::get_rel_update()
        DataType get_rel_update()
        {
          return Math::abs(_sum_alpha)*_norm_dir/_norm_sol;
        }

      protected:
        /// \copydoc NRLinesearch::_apply_intern()
        virtual Status _apply_intern(VectorType& sol, const VectorType& dir)
        {
          // compute initial defect
          Status status(Status::progress);
          this->_num_iter = Index(0);

          // The second secant point in the first iteration is x + _sigma_0 * dir
          _alpha = _sigma_0;

          // Complete update
          _sum_alpha = DataType(0);
          // Norm of the search direction
          _norm_dir = dir.norm2();
          // Norm of the initial guess
          _norm_sol = sol.norm2();

          // The first "other" point for the secant
          this->_vec_tmp.axpy(dir, sol, _sigma_0);
          // Set the state in the operator
          _op.prepare(this->_vec_tmp);
          _op.compute_gradient(this->_vec_tmp);
          _filter.filter_def(this->_vec_tmp);

          _eta = dir.dot(this->_vec_tmp);
          this->_def_init = Math::abs(_eta);

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            this->_op.prepare(sol);
            _op.compute_gradient(this->_vec_tmp);
            _filter.filter_def(this->_vec_tmp);

            // Set new defect, do convergence checks and compute the new _alpha
            status = this->_set_new_defect(dir, this->_vec_tmp);

            // Update the solution
            sol.axpy(dir, sol, _alpha);

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

          }
          return Status::undefined;
        }

        /**
         * \brief Sets new defect norm and does convergence checks
         *
         * \param[in] dir
         * Search direction
         *
         * \param[in] grad
         * Gradient at current iterate
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        Status _set_new_defect(const VectorType& dir, const VectorType& grad) override
        {
          // Increase iteration count
          ++this->_num_iter;

          DataType eta_prev = _eta;

          // Compute eta and alpha
          _eta = grad.dot(dir);
          this->_def_cur = Math::abs(_eta);

          // Check if the diffence in etas is too small, thus leading to a huge update of relative size > sqrt(eps)
          if(Math::abs(_eta - eta_prev) < Math::abs(_alpha*_eta)*_norm_dir*Math::sqrt(Math::eps<DataType>()))
          {
            return Status::stagnated;
          }

          // Update alpha according to secant formula
          _alpha *= _eta/(eta_prev - _eta);
          _sum_alpha += _alpha;

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
            if(Math::sqr(_alpha)*this->_def_init <= Math::eps<DataType>())
            {
              // increment stagnation count
              if(++this->_num_stag_iter >= this->_min_stag_iter)
              {
                return Status::stagnated;
              }
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

    }; // class SecantLinesearch

  } // namespace Solver
} // namespace FEAST
#endif // FEAST_KERNEL_SOLVER_LINESEARCH
