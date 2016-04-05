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
     * \brief Linesearch base class
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
     * for a given search direction \f$ d \f$.
     *
     */
    template<typename Operator_, typename Filter_>
    class Linesearch : public IterativeSolver<typename Operator_::VectorTypeR>
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
        /// Initial solution
        VectorType _vec_initial_sol;
        /// temporary vector
        VectorType _vec_tmp;

        /// Operator functional value
        DataType _fval_min;
        /// Initial functional value
        DataType _fval_0;

        /// Line search parameter
        DataType _alpha_min;
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
         * \param[in] name_
         * String to use in solver plots to console
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
        explicit Linesearch(String name_, Operator_& op_, const Filter_& filter_, bool keep_iterates = false) :
          BaseClass(name_),
          _op(op_),
          _filter(filter_),
          _fval_min(Math::huge<DataType>()),
          _fval_0(Math::huge<DataType>()),
          _alpha_min(DataType(0)),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;
          }

        /// \copydoc ~BaseClass()
        virtual ~Linesearch()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          // Create temporary vectors
          _vec_initial_sol = this->_op.create_vector_r();
          _vec_tmp = this->_op.create_vector_r();
          _vec_grad = this->_op.create_vector_r();
          BaseClass::init_symbolic();
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          // Clear temporary vectors
          _vec_initial_sol.clear();
          _vec_tmp.clear();
          _vec_grad.clear();
          BaseClass::done_symbolic();
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
          return Math::abs(_alpha_min)*_norm_dir/Math::max(_norm_sol, DataType(1));
        }

        /**
         * \brief Gets the functional value of the last iteration
         *
         * \returns The functional value of the last iteration
         */
        DataType get_final_fval() const
        {
          return _fval_min;
        }

        /**
         * \brief Sets the intitial functional value
         *
         * \param[in] f0
         * The intial functional value.
         *
         * This is handy because the linesearch gets called from another solver that in general already evaluated
         * the functional.
         *
         */
        void set_initial_fval(DataType f0)
        {
          _fval_0 = f0;
        }

        /**
         * \brief Sets the initial gradient from a defect vector
         *
         * \param[in] vec_def
         * The intial defect vector.
         *
         * This is handy because the linesearch gets called from another solver that in general already evaluated
         * the functional gradient.
         *
         */
        void get_defect_from_grad(VectorType& vec_def) const
        {
          vec_def.scale(this->_vec_grad, DataType(-1));
        }

        /**
         * \brief Gets a defect vector from the final gradient
         *
         * \param[out] vec_def
         * The final defect vector.
         *
         * This is handy because the other solver can use the gradient of the new solution without re-calculating it.
         *
         */
        void set_grad_from_defect(const VectorType& vec_def)
        {
          this->_vec_grad.scale(vec_def, DataType(-1));
        }

    }; // class Linesearch

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
    class NewtonRaphsonLinesearch : public Linesearch<Operator_, Filter_>
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
        explicit NewtonRaphsonLinesearch(Operator_& op_, const Filter_& filter_, bool keep_iterates = false) :
          BaseClass("NR-LS", op_, filter_, keep_iterates)
          {
          }

        /// \copydoc ~BaseClass()
        virtual ~NewtonRaphsonLinesearch()
        {
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
          return Math::abs(this->_alpha_min)*this->_norm_dir/Math::max(this->_norm_sol, DataType(1));
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

          this->_vec_initial_sol.clone(sol);

          // Norm of the search direction
          this->_norm_dir = dir.norm2();
          // Norm of the initial guess
          this->_norm_sol = sol.norm2();

          this->_fval_min = this->_fval_0;
          this->_alpha_min = DataType(0);

          DataType alpha(0);
          DataType alpha_update(alpha);

          DataType eta = dir.dot(this->_vec_grad);
          if(eta > DataType(0))
            throw InternalError(__func__,__FILE__,__LINE__,"Search direction is not a descent direction: "
                +stringify_fp_sci(eta));

          // Compute initial defect. We want to minimise d^T * grad(_op)
          this->_def_init = Math::abs(dir.dot(this->_vec_grad));

          //sol.axpy(dir, this->_vec_initial_sol, alpha);
          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            // Increase iteration count
            ++this->_num_iter;

            this->_op.prepare(sol, this->_filter);
            this->_op.compute_grad(this->_vec_grad);
            this->_filter.filter_def(this->_vec_grad);
            DataType fval = this->_op.compute_func();

            if(fval < this->_fval_min)
            {
              this->_fval_min = fval;
              this->_alpha_min = alpha;
            }

            // First let's see if we have to compute the defect at all
            bool calc_def = false;
            calc_def = calc_def || (this->_min_iter < this->_max_iter);
            calc_def = calc_def || this->_plot;
            calc_def = calc_def || (this->_min_stag_iter > Index(0));

            // Update defect
            this->_def_cur = Math::abs(this->_vec_grad.dot(dir));

            // Compute new _alpha <- _alpha - grad.dot(dir) / dir.dot(Hess*dir)
            this->_op.apply_hess(this->_vec_tmp, dir);
            this->_filter.filter_def(this->_vec_tmp);

            alpha_update = - this->_vec_grad.dot(dir)/dir.dot(this->_vec_tmp);
            alpha += alpha_update;

            Statistics::add_solver_defect(this->_branch, double(this->_def_cur));
            sol.axpy(dir, this->_vec_initial_sol, alpha);

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

            if(status != Status::progress)
            {
              // If we are not successful, update the solution with the best step found so far
              sol.axpy(dir, this->_vec_initial_sol, this->_alpha_min);

              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

          }

          // We should never come to this point
          return Status::undefined;
        }

    }; // class NewtonRaphsonLinesearch

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
        DataType _sigma_0;
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
          _sigma_0(initial_step_),
          _eta(DataType(0))
          {
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

        /// \copydoc BaseClass::get_rel_update()
        DataType get_rel_update()
        {
          return Math::abs(this->_alpha_min)*this->_norm_dir/this->_norm_sol;
        }

      protected:
        /// \copydoc NewtonRaphsonLinesearch::_apply_intern()
        virtual Status _apply_intern(VectorType& sol, const VectorType& dir)
        {
          // compute initial defect
          Status status(Status::progress);
          this->_num_iter = Index(0);

          this->_vec_initial_sol.clone(sol);

          // Norm of the search direction
          this->_norm_dir = dir.norm2();
          // Norm of the initial guess
          this->_norm_sol = sol.norm2();

          this->_fval_min = this->_fval_0;
          this->_alpha_min = DataType(0);

          // The second secant point in the first iteration is x + _sigma_0 * dir
          DataType alpha(_sigma_0/this->_norm_dir);
          DataType alpha_update(alpha);

          _eta = dir.dot(this->_vec_grad);
          if(_eta > DataType(0))
            throw InternalError(__func__,__FILE__,__LINE__,"Search direction is not a descent direction: "
                +stringify_fp_sci(_eta));

          // The first "other" point for the secant
          sol.axpy(dir, this->_vec_initial_sol, _sigma_0/this->_norm_dir);

          _eta = dir.dot(this->_vec_grad);
          this->_def_init = Math::abs(_eta);

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            // Increase iteration count
            ++this->_num_iter;

            this->_op.prepare(sol, this->_filter);
            DataType fval = this->_op.compute_func();
            this->_op.compute_grad(this->_vec_grad);
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
              return Status::aborted;

            // is diverged?
            if(this->is_diverged())
              return Status::diverged;

            // Check if the diffence in etas is too small, thus leading to a huge update of relative size > sqrt(eps)
            if(Math::abs(_eta - eta_prev) < Math::abs(alpha*_eta)*this->_norm_dir*Math::sqrt(Math::eps<DataType>()))
            {
              // If we are not successful, update the solution with the best step found so far
              sol.axpy(dir, this->_vec_initial_sol, this->_alpha_min);
              return Status::stagnated;
            }

            // Update alpha according to secant formula
            alpha_update *= _eta/(eta_prev - _eta);
            alpha += alpha_update;

            // Update the solution
            sol.axpy(dir, this->_vec_initial_sol, alpha);

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

            if(status != Status::progress)
            {
              // If we are not successful, update the solution with the best step found so far
              sol.axpy(dir, this->_vec_initial_sol, this->_alpha_min);

              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

          }
          return Status::undefined;
        }

    }; // class SecantLinesearch

    /**
     * \brief Strong Wolfe linesearch
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
     *
     */
    template<typename Operator_, typename Filter_>
    class StrongWolfeLinesearch : public Linesearch<Operator_, Filter_>
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
        /// Default tolerance for functional value decrease for the strong Wolfe conditions
        static constexpr DataType tol_decrease_default = DataType(1e-3);
        /// Default tolerance for curvature decrease for the strong Wolfe conditions
        static constexpr DataType tol_curvature_default = DataType(0.3);

      protected:
        /// Last successful line step, needs to start with 1
        DataType _alpha_0;
        /// Hard maximum for the step length
        DataType _alpha_hard_max;
        /// Hard minimum for the step length
        DataType _alpha_hard_min;
        /// Soft maximum for the step length, gets increased if the minimum is not found in the current i.o.u
        DataType _alpha_soft_max;
        /// Soft minimum for the step length, gets increased as better points are found
        DataType _alpha_soft_min;
        /// Initial delta
        DataType _delta_0;
        /// Tolerance for sufficient decrease in the functional value (Wolfe conditions)
        DataType _tol_decrease;
        /// Tolerance for sufficient decrease in the norm of the gradient (Wolfe conditions)
        DataType _tol_curvature;

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
        explicit StrongWolfeLinesearch(Operator_& op_, Filter_& filter_, bool keep_iterates = false,
        DataType tol_decrease_ = tol_decrease_default, DataType tol_curvature_ = tol_curvature_default) :
          BaseClass("SW-LS", op_, filter_, keep_iterates),
          _alpha_0(DataType(1)),
          _alpha_hard_max(DataType(0)),
          _alpha_hard_min(DataType(0)),
          _alpha_soft_max(DataType(0)),
          _alpha_soft_min(DataType(0)),
          _delta_0(Math::huge<DataType>()),
          _tol_decrease(tol_decrease_),
          _tol_curvature(tol_curvature_)
          {
          }

        /// \copydoc ~BaseClass()
        virtual ~StrongWolfeLinesearch()
        {
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "StrongWolfeLinesearch";
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
          auto normed_dir = vec_dir.clone();
          normed_dir.scale(normed_dir, DataType(1)/vec_dir.norm2());

          // apply
          Status st =_apply_intern(vec_sol, normed_dir);

          return st;
        }

        /// \copydoc BaseClass::get_rel_update()
        virtual DataType get_rel_update()
        {
          return this->_alpha_min/this->_norm_dir;
        }

      protected:
        /**
         * \brief Internal function: Applies the solver
         *
         * \param[in, out] vec_sol
         * Initial guess, gets overwritten by solution
         *
         * \param[in] vec_dir
         * Search direction
         *
         * \note This assumes that the initial functional value _fval_0 and the gradient were already set from the
         * calling solver!
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_dir)
        {
          static constexpr DataType extrapolation_width = DataType(4);
          Status status(Status::progress);
          this->_num_iter = Index(0);

          // Save the intial state
          this->_vec_initial_sol.clone(vec_sol);
          // Norm of the search direction vector
          this->_norm_dir = vec_dir.norm2();
          // Norm of the initial guess
          this->_norm_sol = vec_sol.norm2();

          // It is critical that _f_0 was set from the outside!
          DataType fval(this->_fval_0);
          // The functional value from the best step so far
          DataType fval_prev(fval);

          //std::cout << "SWLinesearch start:     sol = " << this->_vec_initial_sol << " dir = " << vec_dir <<
          //  " grad = " << this->_vec_grad << std::endl;
          //std::cout << "SWLinesearch start:     f_0 = " << stringify_fp_sci(this->_fval_0) << " delta_0 = "
          //<< stringify_fp_sci(_delta_0) << std::endl;

          // The step length from the best step so far, so it starts with 0
          DataType alpha_prev(0);
          // The first trial step length. The operator gets prepared and evaluated in the first iteration
          DataType alpha(_alpha_0);

          // Update the solution with the first trial step
          vec_sol.axpy(vec_dir, this->_vec_initial_sol, alpha/this->_norm_dir);

          // Set hard limits to default values if they have not been set
          _alpha_hard_min = DataType(0);
          if(_alpha_hard_max < Math::eps<DataType>())
            _alpha_hard_max = Math::huge<DataType>();

          // Set the soft limits. These are important if we know that the minimum is NOT in the current interval of
          // uncertainity and we need to drive alpha to infinity. See _enforce_step_limits
          _alpha_soft_min = DataType(0);
          _alpha_soft_max = _alpha_0 + extrapolation_width*alpha;
          _alpha_soft_max = Math::min(_alpha_soft_max, _alpha_hard_max);

          // Compute initial defect. We want to minimise <dir^T, grad(_op)>, so everything delta are directional
          // derivatives
          _delta_0 = vec_dir.dot(this->_vec_grad)/this->_norm_dir;
          this->_def_init = Math::abs(_delta_0);
          // The current derivative
          DataType delta(_delta_0);
          // Derivative from the best step so far
          DataType delta_prev(delta);

          //std::cout << "SWLinesearch start: alpha = " << stringify_fp_sci(alpha) << " alpha_prev = " <<
          //  stringify_fp_sci(alpha_prev) << " soft limits = [" <<
          //  stringify_fp_sci(_alpha_soft_min) << ", " << stringify_fp_sci(_alpha_soft_max) << "]" <<std::endl;

          // This is a sanity check
          if(_delta_0 > DataType(0))
            throw InternalError(__func__,__FILE__,__LINE__,"Initial search direction is not a descent direction: "
                +stringify_fp_sci(_delta_0));

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;
            // Increase iteration count
            ++this->_num_iter;

            // Prepare and evaluate
            this->_op.prepare(vec_sol, this->_filter);
            fval = this->_op.compute_func();

            // Compute and filter the gradient
            this->_op.compute_grad(this->_vec_grad);
            this->_filter.filter_def(this->_vec_grad);

            // New directional derivative and new defect. Note that we did not normalise vec_dir, so we have to take
            // care of that here.
            delta = this->_vec_grad.dot(vec_dir)/this->_norm_dir;
            this->_def_cur = Math::abs(delta);

            Statistics::add_solver_defect(this->_branch, double(this->_def_cur));

            //std::cout << "SWLinesearch:   sol = " << vec_sol << " grad = " << this->_vec_grad << std::endl;
            //std::cout << "SWLinesearch: alpha = " << stringify_fp_sci(alpha) << " f   = " <<
            //  stringify_fp_sci(fval) << " delta = " << stringify_fp_sci(delta) << std::endl;
            //std::cout << "       _prev: alpha = " << stringify_fp_sci(alpha_prev) << " f   = " <<
            //  stringify_fp_sci(fval_prev) << " delta = " << stringify_fp_sci(delta_prev) << std::endl;
            //std::cout << " soft limits = [" << stringify_fp_sci(_alpha_soft_min) << ", " <<
            //  stringify_fp_sci(_alpha_soft_max) << "]" <<std::endl;

            // Check if the sufficient decrease condition is violated. If so, one call to _bracket gets the solution
            if( (fval > this->_fval_0 + _tol_decrease*alpha*_delta_0) || (fval > fval_prev && this->_num_iter > 1) )
            {
              //std::cout << "SWLinesearch: sufficient decrease violated" << std::endl;
              //if(fval > this->_fval_0 + _tol_decrease*alpha*_delta_0)
              //{
              //  std::cout << "  f = " << stringify_fp_sci(fval) << " > " <<
              //    stringify_fp_sci(this->_fval_0 + _tol_decrease*alpha*_delta_0) <<
              //    " = _f_0 + _tol_decrease*_alpha*_delta_0)" << std::endl;
              //}

              //if(fval > fval_prev)
              //{
              //  std::cout << "  f = " << stringify_fp_sci(fval) << " > " << stringify_fp_sci(fval_prev ) << " = f_prev" << std::endl;
              //}

              // Call _bracket with alpha_lo = alpha_prev and alpha_hi = _alpha because the old step was better
              status = _bracket(vec_sol, vec_dir, alpha_prev, alpha, delta_prev, delta, fval_prev, fval);

              // _bracket writes the new values for alpha etc. to the first entries of the corresponding set
              alpha = alpha_prev;
              delta = delta_prev;
              fval = fval_prev;

              this->_alpha_min = alpha;
              this->_fval_min = fval;

              // Save the last successful step length for the next call
              _alpha_0 = this->_alpha_min;

            }

            // plot?
            if(this->_plot)
            {
              std::cout << this->_plot_name
              <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
              << " : " << stringify_fp_sci(this->_def_cur)
              << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
              << std::endl;
            }

            // Stop if _bracket was successful or encountered an error
            if(status != Status::progress)
              return status;

            // Ensure that the defect is neither NaN nor infinity
            if(!Math::isfinite(this->_def_cur))
              return Status::aborted;

            // Is diverged?
            if(this->is_diverged())
              return Status::diverged;

            // If we come to here, the sufficient decrease condition is satisfied and the solver did not encounter
            // any errors

            // If the curvature condition is satisfied, the current alpha is good and we are successful
            if( this->_def_cur < -_tol_curvature*_delta_0)
            {
              //std::cout << "SWLinesearch: Sufficient decrease condition fulfilled" << std::endl;
              //std::cout << "   f = " << stringify_fp_sci(fval) << " < " <<
              //  stringify_fp_sci(this->_fval_0 + _tol_decrease*alpha*_delta_0) << std:: endl;
              //std::cout << "   and curvature condition satisfied: " << stringify_fp_sci(this->_def_cur)
              //<< " < " << stringify_fp_sci(-_tol_curvature*_delta_0) << std::endl;

              // Save minimum parameters
              this->_alpha_min = alpha;
              this->_fval_min = fval;
              // Save the last successful step length for the next call
              _alpha_0 = this->_alpha_min;

              return Status::success;
            }

            // If the derivative changes sign, we are successful as well
            if(delta >= DataType(0))
            {
              //std::cout << "SWLinesearch: derivative is positive! " << stringify_fp_sci(delta) << std::endl;
              status = _bracket(vec_sol, vec_dir, alpha, alpha_prev, delta, delta_prev, fval, fval_prev);

              // Save minimum parameters
              this->_alpha_min = alpha;
              this->_fval_min = fval;
              // Save the last successful step length for the next call
              _alpha_0 = this->_alpha_min;

              return status;
            }

            // If we come to here, we do not know if the minimum is in the current interval of uncertainity, so we
            // need to update it. _polynomial_fit does that for us.

            DataType alpha_new(0);

            bool min_in_interval(false);

            // We discard the status returned by _polynomial_fit because it means nothing at this point
            _polynomial_fit(alpha_new, alpha, alpha_prev, fval, fval_prev, delta, delta_prev, min_in_interval,
            (delta < 0));
            _enforce_step_limits(alpha_new, alpha, alpha_prev);

            // Since the sufficient decrease condition was satisfied, we update alpha_prev
            alpha_prev = alpha;
            fval_prev = fval;
            delta_prev = delta;

            // Update alpha. fval and delta get updated at the start of the next iteration.
            alpha = alpha_new;

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

            // If the maximum number of iterations was performed, return the iterate for the best step so far
            if(this->_num_iter > this->_max_iter)
            {
              vec_sol.axpy(vec_dir, this->_vec_initial_sol, alpha_prev/this->_norm_dir);
              return Status::max_iter;
            }

            // Update solution: sol <- sol + _alpha*dir
            vec_sol.axpy(vec_dir, this->_vec_initial_sol, alpha/this->_norm_dir);

          }

          // We should never come to this point
          return Status::undefined;
        }

        /**
         * \brief Bracketing iteration that finds a minimum of a 1d function
         *
         * \param[in, out] vec_sol
         * The current solution, gets updated in the iteration
         *
         * \param[in] vec_dir
         * The search direction
         *
         * \param[in,out] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in,out] alpha_hi
         * The step size for the higher functional value
         *
         * \param[in,out] delta_lo
         * Derivative at alpha_lo
         *
         * \param[in,out] delta_hi
         * Derivative at alpha_hi
         *
         * \param[in,out] f_lo
         * Functional value at alpha_lo
         *
         * \param[in,out] f_hi
         * Functional value at alpha_hi
         *
         * \returns
         * Status::success if no abnormal termination (i.e. triggered by _polynomial_fit) occured
         *
         */
        Status _bracket(VectorType& vec_sol, const VectorType& vec_dir, DataType& alpha_lo, DataType& alpha_hi,
        DataType& delta_lo, DataType& delta_hi, DataType& f_lo, DataType& f_hi)
        {
          // Check if the starting point is already good enough. This should have been checked before, anyway.
          if(f_lo < this->_fval_0 +_tol_decrease*alpha_lo*_delta_0 && Math::abs(delta_lo) < -_tol_curvature*_delta_0)
            return Status::success;

          // If the derivative at the inital step length with the lower functional value is positive, we are
          // searching "backwards" and need to pass this information to _polynomial_fit
          bool df_was_negative(delta_lo < DataType(0));
          // We need to remember if the minimum is already known to be in the interval of uncertainity for subsequent
          // calls to _polynomial_fit
          bool min_in_interval(false);

          // Trial step and functional value/gradient there
          DataType alpha_new(0);
          DataType f_new(0);
          DataType delta_new(0);

          Status st(Status::progress);
          int iter(0);
          while(st == Status::progress && iter < 100)
          {
            ++iter;

            // Find new trial step
            st = _polynomial_fit(alpha_new, alpha_lo, alpha_hi, f_lo, f_hi, delta_lo, delta_hi, min_in_interval,
            df_was_negative);
            _enforce_step_limits(alpha_new, alpha_lo, alpha_hi);

            // Update the solution and evaluate functional/gradient there
            vec_sol.axpy(vec_dir, this->_vec_initial_sol, alpha_new/this->_norm_dir);
            this->_op.prepare(vec_sol, this->_filter);

            f_new = this->_op.compute_func();
            this->_op.compute_grad(this->_vec_grad);
            delta_new = this->_vec_grad.dot(vec_dir);

            //std::cout << "Bracket iter " << iter << " sol = " << vec_sol << " dir = " << vec_dir
            //<< " grad = " << this->_vec_grad << std::endl;
            //std::cout << "  alpha_new " << stringify_fp_sci(alpha_new) << " f_new " << stringify_fp_sci(f_new) << " delta_new " << stringify_fp_sci(delta_new) << std::endl;
            //std::cout << "  alpha_lo  " << stringify_fp_sci(alpha_lo) << " f_lo  " << stringify_fp_sci(f_lo) << " delta_lo  " << stringify_fp_sci(delta_lo) << std::endl;
            //std::cout << "  alpha_hi  " << stringify_fp_sci(alpha_hi) << " f_hi  " << stringify_fp_sci(f_hi) << " delta_hi  " << stringify_fp_sci(delta_hi) << std::endl;

            // Update of interval of uncertainity. We might already be successful.

            // Check if the sufficient decrease condition is violated
            if(f_new > this->_fval_0 +_tol_decrease*alpha_new*_delta_0 || ( (f_new > f_lo) && iter > 1) )
            {
              //std::cout << "  sufficient decrease violated" << std::endl;
              //std::cout << "  f_0 + tol_decrease*alpha_new*_delta_0 " <<
              //  stringify_fp_sci(this->_fval_0 +_tol_decrease*alpha_new*_delta_0)
              //  << " f_lo = " << stringify_fp_sci(f_lo) << "f_new = " << stringify_fp_sci(f_new) << std::endl;

              alpha_hi = alpha_new;
              f_hi = f_new;
              delta_hi = delta_new;
            }
            // If we have a sufficient decrease, we need to check the curvature condition
            else
            {
              //std::cout << "  bracket: sufficient decrease satisfied:  f = " << stringify_fp_sci(f_new) << " <= "
              //<< stringify_fp_sci(this->_fval_0 +_tol_decrease*alpha_new*_delta_0) << std::endl;

              // If the curvature condition is satisfied too, we are successful
              if(Math::abs(delta_new) < -_tol_curvature*_delta_0)
              {
                //std::cout << "           curvature condition satisfied " << stringify_fp_sci(Math::abs(delta_new))
                //<< " < "  << stringify_fp_sci(-_tol_curvature*_delta_0) << std::endl;

                // Ok, alpha_new contains a better step than alpha_lo
                alpha_lo = alpha_new;
                f_lo = f_new;
                delta_lo = delta_new;

                st = Status::success;
              }
              // If the curvature condition is violated, we replace either alpha_lo or alpha_hi by alpha_new
              else
              {
                //std::cout << "  curvature condition violated: " << stringify_fp_sci(Math::abs(delta_new)) << " > "
                //<< stringify_fp_sci(-_tol_curvature*_delta_0) << std::endl;

                // This is the original condition from GN06
                //if(delta_new*(alpha_hi - alpha_lo) >= DataType(0))
                // This seems to work better and is uses in ALGLIB
                if(delta_new < DataType(0))
                {
                  alpha_hi = alpha_lo;
                  f_hi = f_lo;
                  delta_hi = delta_lo;
                }
              }

              // Since the sufficient decrease condition was satisfied, we have a new alpha_lo
              alpha_lo = alpha_new;
              f_lo = f_new;
              delta_lo = delta_new;

            }
          }

          // The bracketing iterations ends

          // If we have not been successful, use the best step size obtained so far.
          if(st != Status::success)
          {
            vec_sol.axpy(vec_dir, this->_vec_initial_sol, alpha_lo/this->_norm_dir);

            this->_op.prepare(vec_sol, this->_filter);
            this->_op.compute_grad(this->_vec_grad);
          }
          else
          {
            _alpha_0 = alpha_lo;
          }

          //std::cout << "Bracket stops at sol = " << vec_sol << " grad = " << this->_vec_grad << std::endl;
          //std::cout << "  alpha_lo = " << stringify_fp_sci(alpha_lo) << " f_lo = "<< stringify_fp_sci(f_lo)
          //<< " delta_lo = " << stringify_fp_sci(delta_lo) << std::endl;
          //std::cout << "  alpha_hi = " << stringify_fp_sci(alpha_hi) << " f_hi = "<< stringify_fp_sci(f_hi)
          //<< " delta_hi = " << stringify_fp_sci(delta_hi) << std::endl;

          return st;
        } // _bracket

        /**
         * \brief The great magick trick to find a minimum of a 1d function
         *
         * \param[out] alpha_new
         * The step corresponding to the minumum
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] alpha
         * The step size for the higher functional value
         *
         * \param[in] df_lo
         * Derivative at alpha_lo
         *
         * \param[in] df
         * Derivative at alpha_hi
         *
         * \param[in] f_lo
         * Functional value at alpha_lo
         *
         * \param[in] f
         * Functional value at alpha_hi
         *
         * \param[in,out] min_in_interval
         * Do we know the minimum has to be in the interval (or any superset we worked on before)? Because this gets
         * called several times, this knowledge has to be passed back and forth.
         *
         * \param[in] df_was_negative
         * If the first df_lo was POSITIVE, we have some sort of backwards search (decreasing step lengths) and this
         * changes how the cases get treated (see Case 1).
         *
         * \returns
         * Status::success as there is nothing that can fail.
         *
         * This routine does a quadratic interpolation using f_lo, df_lo and df, and a cubic interpolation using
         * f_lo, df_lo, f_hi and df_hi. The great magick is to determine which is better, and this depends on a lot
         * of factors, which leads to 4 distinct cases.
         *
         * This is also used if we do not know in advance that the minimum is in the given interval, so it performs
         * extrapolation as well.
         *
         */
        Status _polynomial_fit(DataType& alpha_new, DataType alpha_lo, DataType alpha, DataType f_lo, DataType f, DataType df_lo, DataType df, bool& min_in_interval, bool df_was_negative) const
        {
          //std::cout << "Polynomial_fit: min_in_interval = " << min_in_interval << std::endl;
          //std::cout << "  alpha_lo = " << stringify_fp_sci(alpha_lo) << " f_lo = "<< stringify_fp_sci(f_lo)
          //<< " df_lo = " << stringify_fp_sci(df_lo) << std::endl;
          //std::cout << "  alpha    = " << stringify_fp_sci(alpha) << " f    = "<< stringify_fp_sci(f)
          //<< " df    = " << stringify_fp_sci(df) << std::endl;

          // Default: Quadratic interpolant using f_lo, df_lo and df
          DataType alpha_q = _argmin_quadratic(alpha_lo, alpha, f_lo, f, df_lo, df, true);
          DataType alpha_c = _argmin_cubic(alpha_lo, alpha, f_lo, f, df_lo, df);

          // Case 1: This basically means that alpha was the new step and the function value increased.
          if( ((alpha_lo-alpha) < DataType(0)) && df_was_negative )
          {
            // The first derivative was negative, but the function value increases, so the minimum has to be in this
            // interval (or have been in the original interval, which this interval is a subset of)
            min_in_interval = true;
            //std::cout << "  Case 1 : ";
            // use diffent 2nd order approximation to mimick ALGLIB
            alpha_q = _argmin_quadratic(alpha_lo, alpha, f_lo, f, df_lo, df, false);
            // If the cubic minimum is closer to the old best step length, take it. Otherwise, take the mean
            if(Math::abs(alpha_c - alpha_lo) < Math::abs(alpha_q - alpha_lo))
              alpha_new = alpha_c;
            else
              alpha_new = (alpha_q + alpha_c)/DataType(2);

          }
          // Case 2: The first derivative changes sign
          // We also know that the function value increased in the last trial step
          else if( df_lo*df < DataType(0) )
          {
            // Because the derivative changed sign, it has to be zero somewhere in between, so the minimum is in the
            // current interval
            min_in_interval = true;
            //std::cout << "  Case 2 : ";
            // Take the step closer to the new step alpha
            Math::abs(alpha - alpha_c) > Math::abs(alpha - alpha_q) ? alpha_new = alpha_q : alpha_new = alpha_c;
          }
          // Case 3: The absolute value of the derivative increases
          // We also know that the function value increased in the last trial step and that the derivative did not
          // change sign.
          // This means we have to further search in the direction of alpha.
          else if(Math::abs(df) > Math::abs(df_lo))
          {
            //std::cout << "  Case 3";
            // If the cubic step is closer to the trial step
            if(Math::abs(alpha - alpha_c) < Math::abs(alpha - alpha_q))
            {
              //std::cout << "a: ";
              min_in_interval ? alpha_new = alpha_c : alpha_new = alpha_q;
            }
            else
            {
              //std::cout << "b: ";
              min_in_interval ? alpha_new = alpha_q : alpha_new = alpha_c;
            }
          }
          // Case 4: The absolute value of the derivative did not increase
          // We also know that the function value increased in the last trial step and that the derivative did not
          // change sign.
          else
          {
            //std::cout << "  Case 4: ";
            // Here we can just take the cubic step because the cubic interpolation sets it to the alpha_soft_max
            // if it recognises that the minimum is outside the current interval.
            alpha_new = alpha_c;
          }

          //std::cout << "alpha_new = " << stringify_fp_sci(alpha_new) <<
          //" alpha_c = " << stringify_fp_sci(alpha_c) << " alpha_q = " << stringify_fp_sci(alpha_q) << " min_in_interval " << min_in_interval << std::endl;
          //
          return Status::success;
        }

        /**
         * \brief Enforces hard and soft step limits, adjusting the soft limits if neccessary
         *
         * \param[out] alpha_new
         * New step that satisfies \f$ \alpha_{\mathrm{new}} \in [\alpha_{\mathrm{min}}, \alpha_{\mathrm{max}}] \f$
         * and \f$ \alpha_{\mathrm{new}} \in [\alpha_{\mathrm{soft_min}}, \alpha_{\mathrm{soft_max}}] \f$
         *
         * \param[in] alpha
         * Step size to enforce the limits for
         *
         * \param[in] alpha_prev
         * Best step found so far
         *
         * If alpha > alpha_soft_max, alpha_soft_min and alpha_soft_max get increased according to a heuristic to
         * make sure that we can increase alpha fast enought that the linesearch does not spend too much time on
         * tiny trial steps.
         */
        void _enforce_step_limits(DataType& alpha_new, DataType alpha, DataType alpha_prev)
        {
          alpha_new = Math::max(alpha_new, _alpha_hard_min);
          alpha_new = Math::min(alpha_new, _alpha_hard_max);

          if(alpha_new >= _alpha_soft_max)
          {
            alpha_new = _alpha_soft_max;
            //std::cout << "SWLinesearch: increase soft min from " << stringify_fp_sci(_alpha_soft_min) << " to ";
            _alpha_soft_min = alpha;
            //std::cout << stringify_fp_sci(_alpha_soft_min) << std::endl;
            //std::cout << "SWLinesearch: increase soft max from  " << stringify_fp_sci(_alpha_soft_max) << " to ";
            _alpha_soft_max = alpha + DataType(4)*(alpha_new - alpha_prev);
            //std::cout << stringify_fp_sci(_alpha_soft_max) << std::endl;
          }

        }

        /**
         * \brief Computes the minimum of a quadratic interpolation polynomial
         *
         * This has two variants:
         *  a) Interpolate f_lo, df_lo and f_hi
         *  b) Interpolate f_lo, df_lo and df_hi
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] alpha_hi
         * The step size for the higher functional value
         *
         * \param[in] f_lo
         * Functional value at alpha_lo
         *
         * \param[in] f_hi
         * Functional value at alpha_hi
         *
         * \param[in] df_lo
         * Derivative at alpha_lo
         *
         * \param[in] df_hi
         * Derivative at alpha_hi
         *
         * \param[in] interpolate_derivative
         * If this is true, variant a) is chosen
         *
         * \returns
         * The step corresponding to the minumum
         *
         */
        DataType _argmin_quadratic(DataType alpha_lo, DataType alpha_hi, DataType f_lo, DataType f_hi, DataType df_lo, DataType df_hi, bool interpolate_derivative) const
        {
          DataType alpha(alpha_lo);

          // Quadratic interpolation using f_lo, df_lo, df_hi
          if(interpolate_derivative)
            alpha -= df_lo * (alpha_hi - alpha_lo)/(df_hi - df_lo);
          // Quadratic interpolation using f_lo, f_hi, df_lo
          else
            alpha -= df_lo*(alpha_hi - alpha_lo)/(DataType(2)*( (f_hi - f_lo)/(alpha_hi - alpha_lo) - df_lo ) );

          return alpha;
        }

        /**
         * \brief Computes the minimum of a cubic interpolation polynomial
         *
         * If the minimum of the cubic interpolation polynomial is in the interiour, it is computed here. If it is
         * not in the interiour, it lies across from the point with the lower function value. In this case, we set
         * the return value to _alpha_soft_{min, max} accordingly.
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] alpha_hi
         * The step size for the higher functional value
         *
         * \param[in] f_lo
         * Functional value at alpha_lo
         *
         * \param[in] f_hi
         * Functional value at alpha_hi
         *
         * \param[in] df_lo
         * Derivative at alpha_lo
         *
         * \param[in] df_hi
         * Derivative at alpha_hi
         *
         * \returns
         * The step corresponding to the minumum or _alpha_soft_{min, max} if there is no minimum in the interiour.
         *
         */
        DataType _argmin_cubic(DataType alpha_lo, DataType alpha_hi, DataType f_lo, DataType f_hi, DataType df_lo, DataType df_hi) const
        {
          DataType alpha(alpha_lo);

          DataType d1 = - DataType(3)*(f_hi - f_lo)/(alpha_hi - alpha_lo) + (df_hi + df_lo);

          // Scale the computation of r for better numerical stability
          DataType scale = Math::max(Math::abs(df_lo), Math::abs(df_hi));
          scale = Math::max(scale, Math::abs(d1));

          DataType r = Math::sqr(d1/scale) - (df_lo/scale) * (df_hi/scale);

          DataType d2(0);

          if(r > DataType(0))
          {
            d2 = Math::signum(alpha_lo - alpha_hi) * scale *
              Math::sqrt( Math::sqr(d1/scale) - (df_lo/scale) * (df_hi/scale));

            alpha -= (alpha_lo - alpha_hi)*(-d1 + df_lo + d2)/(d2 - df_hi + df_lo + d2);
          }
          // r <= 0 means that the minimum is not in the interiour, so it has to lie across the endpoint with the
          // lower function value
          else
          {
            alpha_lo < alpha_hi ? alpha = _alpha_soft_min : alpha = _alpha_soft_max;
          }

          return alpha;
        }
    }; // class StrongWolfeLinesearch

    /**
     * \brief Creates a new NewtonRaphsonLinesearch object
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new NewtonRaphsonLinesearch object.
     */
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<NewtonRaphsonLinesearch<Operator_, Filter_>> new_newton_raphson_linesearch(
      Operator_& op, Filter_& filter, bool keep_iterates = false)
      {
        return std::make_shared<NewtonRaphsonLinesearch<Operator_, Filter_>>(op, filter, keep_iterates);
      }

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
     * \brief Creates a new StrongWolfeLinesearch object
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new StrongWolfeLinesearch object.
     */
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<StrongWolfeLinesearch<Operator_, Filter_>> new_strong_wolfe_linesearch(
      Operator_& op, Filter_& filter, bool keep_iterates = false,
      typename Operator_::DataType tol_decrease = StrongWolfeLinesearch<Operator_, Filter_>::tol_decrease_default,
      typename Operator_::DataType tol_curvature = StrongWolfeLinesearch<Operator_, Filter_>::tol_curvature_default)
      {
        return std::make_shared<StrongWolfeLinesearch<Operator_, Filter_>>
          (op, filter, keep_iterates, tol_decrease, tol_curvature);
      }

  } // namespace Solver
} // namespace FEAST
#endif // FEAST_KERNEL_SOLVER_LINESEARCH
