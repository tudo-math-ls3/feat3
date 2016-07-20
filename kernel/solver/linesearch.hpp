#pragma once
#ifndef FEAT_KERNEL_SOLVER_LINESEARCH
#define FEAT_KERNEL_SOLVER_LINESEARCH 1
#include <kernel/base_header.hpp>
#include <kernel/solver/iterative.hpp>
#include <bitset>

#include <deque>

namespace FEAT
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
        Filter_& _filter;

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
        /// Threshold for trimming function value and gradient
        DataType _trim_threshold;

        /// Line search parameter
        DataType _alpha_min;
        /// The 2-norm of the search direction
        DataType _norm_dir;
        /// The 2-norm of the iterate
        DataType _norm_sol;

        /// Tolerance for the update step
        DataType _tol_step;

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
        explicit Linesearch(String name_, Operator_& op_, Filter_& filter_, bool keep_iterates = false) :
          BaseClass(name_),
          _op(op_),
          _filter(filter_),
          _fval_min(Math::huge<DataType>()),
          _fval_0(Math::huge<DataType>()),
          _trim_threshold(Math::huge<DataType>()),
          _alpha_min(DataType(0)),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          _tol_step(DataType(Math::pow(Math::eps<DataType>(), DataType(0.85)))),
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
         * \brief Resets various member variables in case the solver is reused
         */
        virtual void reset()
        {
          _fval_min = Math::huge<DataType>();
          _fval_0 = Math::huge<DataType>();
          _trim_threshold = Math::huge<DataType>();
          _alpha_min = DataType(0);
          _norm_dir = DataType(0);
          _norm_sol = DataType(0);

          if(iterates != nullptr)
            iterates->clear();
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
        virtual DataType get_rel_update()
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
          if(_trim_threshold == Math::huge<DataType>())
            _trim_threshold = DataType(10)*(Math::abs(_fval_0) + DataType(1));
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

        /**
         * \brief Gets the initial solution the linesearch started with
         *
         * \returns A const reference to the initial solution vector.
         *
         * This is useful for rejecting steps that are too small in the calling solver, like NLCG or NLSD.
         */
        const VectorType& get_initial_sol() const
        {
          return _vec_initial_sol;
        }

        /**
         * \brief Trims the function value and gradient according to some threshhold
         *
         * \param[in,out] func
         * Function value.
         *
         * This sets
         * \f[
         *   trim(f, \nabla f) =
         *     \begin{cases}
         *       (f, \nabla f), & f \leq f_t \\
         *       (f_t, 0 ), & f > f_t
         *     \end{cases},
         * \f]
         * where \f$ f_t : = 10*(|f_0|+1) \f$.
         *
         * This is useful and even mandatory when \f$ f \f$ has singularities, so that if the linesearch steps into
         * a singularity the very high function and gradient values do not pollute the solution process, since they
         * get rejected by the linesearch anyway but might be used in computations.
         *
         */
        virtual void trim_func_grad(DataType& func)
        {
          if(func > _trim_threshold)
          {
            func = _trim_threshold;
            this->_vec_grad.format(DataType(0));
          }
        }

    }; // class Linesearch

    /**
     * \brief Fixed step line search
     *
     * This does not perform a line search, it just returns a constant step size
     *
     */
    template<typename Operator_, typename Filter_>
    class FixedStepLinesearch : public Linesearch<Operator_, Filter_>
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

      private:
        DataType _steplength;

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
        explicit FixedStepLinesearch(Operator_& op_, Filter_& filter_, const DataType steplength_,
        bool keep_iterates = false) :
          BaseClass("FS-LS", op_, filter_, keep_iterates),
          _steplength(steplength_)
          {
            this->set_max_iter(0);
            this->_alpha_min = _steplength;
          }

        /// \copydoc ~BaseClass()
        virtual ~FixedStepLinesearch()
        {
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "Fixed-Step-Linesearch";
        }

        /// \copydoc BaseClass::get_rel_update()
        virtual DataType get_rel_update() override
        {
          return _steplength;
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
          vec_cor.axpy(vec_dir, vec_cor, _steplength);

          this->_op.prepare(vec_cor, this->_filter);

          return Status::success;
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
          vec_sol.axpy(vec_dir, vec_sol, _steplength);

          this->_op.prepare(vec_sol, this->_filter);

          this->_fval_min = this->_op.compute_func();
          this->_op.compute_grad(this->_vec_grad);
          this->_filter.filter_def(this->_vec_grad);

          return Status::success;
        }
    }; // class FixedStepLinesearch

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
        explicit NewtonRaphsonLinesearch(Operator_& op_, Filter_& filter_, bool keep_iterates = false) :
          BaseClass("NR-LS", op_, filter_, keep_iterates)
          {
            this->set_max_iter(20);
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
          DataType alpha_hidate(alpha);

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
            IterationStats stat(*this);

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

            // Update defect
            this->_def_cur = Math::abs(this->_vec_grad.dot(dir));

            // Compute new _alpha <- _alpha - grad.dot(dir) / dir.dot(Hess*dir)
            this->_op.apply_hess(this->_vec_tmp, dir);
            this->_filter.filter_def(this->_vec_tmp);

            alpha_hidate = - this->_vec_grad.dot(dir)/dir.dot(this->_vec_tmp);
            alpha += alpha_hidate;

            //Statistics::add_solver_defect(this->_branch, double(this->_def_cur));
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

            // If we are not successful, update the solution with the best step found so far
            if(status != Status::progress)
            {
              sol.axpy(dir, this->_vec_initial_sol, this->_alpha_min);
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
            this->set_max_iter(20);
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
          DataType alpha_hidate(alpha);

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
            IterationStats stat(*this);

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
     * This is quite complicated code in some places, and very sensitive to annihilation as lots of information is
     * condensed into very few floating point values (i.e. the dot products of the gradient and the search direction
     * for a highly nonlinear and high dimensional problem).
     *
     * The line search works with two potentially overlapping intervals: The interval of uncertainty
     * \f$ [\alpha_{\mathrm{soft~min}}, \alpha_{\mathrm{soft~max}}] \f$ and the current search interval
     * \f$ [\alpha_{\mathrm{lo}}, \alpha_{\mathrm{hi}}] \f$.
     *
     * Intitially, we have an initial step length \f$ \alpha_0 \f$ which defaults to 1 if it was not set explicitly
     * before and set \f$ \alpha_{\mathrm{soft~min}} = 0, \alpha_{\mathrm{soft~max}} = \alpha_0 + c(\alpha_0 -
     * \alpha_{\mathrm{soft~min}}) \f$, where the constant is \f$ c = 4\f$ which is pretty arbitrary.
     * The search interval is the degenerate interval \f$ [\alpha_0, \alpha_0]\f$. Not that we know that the
     * derivative \f$ \left< d, grad(\alpha_0) \right> < 0 \f$ because \f$ d \f$ is a descent direction.
     *
     * At this point we do not know if the minimum we are looking for is in either interval (it is not, in general),
     * so we first need to expand the interval of uncertainty until we are sure that it contains the minimum.
     * This is the case if i.e. the derivative changes sign for some step size \f$ \alpha > \alpha_0 \f$, or if the
     * function value decreases so the sufficient decrease condition is satisfied along with the curvature condition.
     *
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
        /// Lower bound of the interval of uncertainty
        DataType _alpha_soft_max;
        /// Upper bound of the interval of uncertainty
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
         * \param[in] tol_decrease_
         * Tolerance for sufficient decrease in function value.
         *
         * \param[in] tol_curvature_
         * Tolerance for the curvature condition.
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
            this->set_max_iter(20);
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

        /// \copydoc BaseClass::reset()
        virtual void reset() override
        {
          BaseClass::reset();
          _alpha_0 = DataType(1);
          _alpha_hard_max = DataType(0);
          _alpha_hard_min = DataType(0);
          _alpha_soft_max = DataType(0);
          _alpha_soft_min = DataType(0);
          _delta_0 = Math::huge<DataType>();
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
        virtual DataType get_rel_update() override
        {
          return this->_alpha_min*this->_norm_dir;
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
          this->_vec_initial_sol = vec_sol.clone(LAFEM::CloneMode::Deep);
          // Norm of the search direction vector. It was normalised befor, but there might be rounding errors
          this->_norm_dir = vec_dir.norm2();
          // Norm of the initial guess
          this->_norm_sol = vec_sol.norm2();

          // Compute initial defect. We want to minimise <dir^T, grad(_op)>, so everything delta are directional
          // derivatives
          _delta_0 = vec_dir.dot(this->_vec_grad);
          this->_def_init = Math::abs(_delta_0);
          XASSERTM(_delta_0 < DataType(0),"Initial search direction is not a descent direction!");

          // Set hard limits to default values if they have not been set
          _alpha_hard_min = DataType(0);
          // This bogus value should be around 1e50 for double precision. It is chosen make comparing results with
          // ALGLIB easier
          if(_alpha_hard_max < Math::eps<DataType>())
            _alpha_hard_max = Math::pow(Math::huge<DataType>(), DataType(0.1622));

          // Set the intervall of uncertainity
          _alpha_soft_min = DataType(0);
          _alpha_soft_max = _alpha_0 + extrapolation_width*_alpha_0;
          _alpha_soft_max = Math::min(_alpha_soft_max, _alpha_hard_max);

          //std::cout << "Linesearch initial alpha " << _alpha_0 << std::endl;

          DataType alpha_lo(0);
          // It is critical that _f_0 was set from the outside!
          DataType fval_lo(this->_fval_0);
          DataType df_lo(_delta_0);

          DataType alpha_hi(0);
          DataType fval_hi(this->_fval_0); // = this->_op.compute_func();
          DataType df_hi(_delta_0); // = this->_vec_grad.dot(vec_dir);

          // Set the first step
          DataType alpha(_alpha_0);
          DataType fval(0);
          DataType df(0);

          // This is the width of the search interval
          DataType width(Math::abs(_alpha_hard_max - _alpha_hard_min));
          DataType width_old(DataType(2)*width);

          // Do we know the interval of uncertainty?
          bool interval_known(false);
          // Is the minimum in the interval of uncertainty?
          bool min_in_interval(false);
          // Does the minimum lie outside the current search interval, requiring us to drive the step to its
          // boundary?
          bool drive_to_bndry(false);

          // start iterating
          while(status == Status::progress)
          {

            IterationStats stat(*this);
            ++(this->_num_iter);

            // If we know the minimum is in the search interval, the interval of uncertainty is the search inverval
            if(min_in_interval)
            {
              _alpha_soft_min = Math::min(alpha_lo, alpha_hi);
              _alpha_soft_max = Math::max(alpha_lo, alpha_hi);
            }
            // Enlarge the interval of uncertainty
            else
            {
              _alpha_soft_min = alpha_lo;
              _alpha_soft_max = alpha + extrapolation_width*Math::abs(alpha - alpha_lo);
            }
            //std::cout << "Set alpha_smin " << _alpha_soft_min << " alpha_smax " << _alpha_soft_max << std::endl;

            // Update solution: sol <- initial_sol + _alpha*dir
            vec_sol.axpy(vec_dir, this->_vec_initial_sol, alpha);
            //std::cout << "Linesearch alpha " << alpha << std::endl;
            //std::cout << "initial_sol " << *(this->_vec_initial_sol) << std::endl;
            //std::cout << "dir " << *vec_dir << std::endl;
            //std::cout << "sol " << *vec_sol << std::endl;

            // Prepare and evaluate
            this->_op.prepare(vec_sol, this->_filter);
            fval = this->_op.compute_func();

            // Compute and filter the gradient
            this->_op.compute_grad(this->_vec_grad);
            this->trim_func_grad(fval);
            this->_filter.filter_def(this->_vec_grad);

            // New directional derivative and new defect
            df = this->_vec_grad.dot(vec_dir);
            this->_def_cur = Math::abs(df);

            //Statistics::add_solver_defect(this->_branch, double(this->_def_cur));

            // plot?
            if(this->_plot)
            {
              std::cout << this->_plot_name
              <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
              << " : " << stringify_fp_sci(this->_def_cur)
              << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
              << " : " << stringify_fp_sci(fval) << " : " << stringify_fp_sci(alpha)
              << std::endl;
            }

            status = Status::progress;

            // Stagnation due to rounding errors
            //if(min_in_interval && (alpha <= _alpha_soft_min || alpha >= _alpha_soft_max))
            //{
            //  std::cout << "Rounding errors" << std::endl;
            //  status = Status::stagnated;
            //}
            //if( alpha == _alpha_hard_max)
            //if( alpha == _alpha_hard_min)

            // If the maximum number of iterations was performed, return the iterate for the best step so far
            if(this->_num_iter >= this->_max_iter)
            {
              //std::cout << "nfevmax " << this->_num_iter-1 << " >= " << this->_max_iter-1 << std::endl;
              status = Status::max_iter;
            }

            if(min_in_interval && (_alpha_soft_max - _alpha_soft_min) <= this->_tol_step*_alpha_soft_max)
            {
              //std::cout << "interval width " << _alpha_soft_max - _alpha_soft_min << " : " << this->_tol_step*_alpha_soft_max << std::endl;
              status = Status::stagnated;
            }

            // If the strong Wolfe conditions hold, we are successful
            if(fval < this->_fval_0 +_tol_decrease*alpha*_delta_0
                && Math::abs(df) < -_tol_curvature*_delta_0)
                {
                  status = Status::success;
                }

            // Stop if _bracket was successful or encountered an error
            if(status == Status::success)
            {
              this->_vec_tmp.axpy(vec_sol, this->_vec_initial_sol, -DataType(1));
              if(fval >= this->_fval_0 || this->_vec_tmp.norm2() == DataType(0))
                status = Status::stagnated;
            }

            if(status != Status::progress)
              break;

            if(!interval_known
                && (fval < this->_fval_0 + this->_tol_decrease*alpha*_delta_0)
                && (df >= Math::min(this->_tol_decrease, this->_tol_curvature)))
                {
                  interval_known = true;
                }

            //std::cout << "interval known " << interval_known << std::endl;
            // If we do not know that the minimum was in the previous interval of uncertainty, we need to compute a
            // new step size to expand the interval of uncertainty at the start of the next iteration.
            if(!interval_known && (fval <= fval_lo) && fval > this->_fval_0 + alpha*this->_tol_decrease*_delta_0)
            {
              DataType fval_m(fval - alpha*this->_tol_decrease*_delta_0);
              DataType df_m(df - this->_tol_decrease*_delta_0);
              DataType fval_lo_m(fval_lo - alpha_lo*this->_tol_decrease*_delta_0);
              DataType df_lo_m(df_lo - this->_tol_decrease*_delta_0);
              DataType fval_hi_m(fval_hi - alpha_hi*this->_tol_decrease*_delta_0);
              DataType df_hi_m(df_hi - this->_tol_decrease*_delta_0);

              // Note that the expansion step might already give us the information that the minimum is the the
              // new search interval
              _polynomial_fit(
                alpha, fval_m, df_m, alpha_lo, fval_lo_m, df_lo_m,
                alpha_hi, fval_hi_m, df_hi_m, min_in_interval, drive_to_bndry);

              fval_lo = fval_lo_m + alpha_lo*this->_tol_decrease*_delta_0;
              df_lo = df_lo_m + this->_tol_decrease*_delta_0;
              fval_hi = fval_hi_m * alpha_hi*this->_tol_decrease*_delta_0;
              df_hi = df_hi_m + this->_tol_decrease*_delta_0;

            }
            else
            {
              _polynomial_fit(alpha, fval, df, alpha_lo, fval_lo, df_lo, alpha_hi, fval_hi, df_hi, min_in_interval,
              drive_to_bndry);

              //std::cout << "width " << width << " width_old " << width_old << " min_in_interval " << min_in_interval <<std::endl;

              if(min_in_interval)
              {
                if(Math::abs(alpha_hi - alpha_lo) >= DataType(0.66)*width_old)
                {
                  //std::cout << "Forcing " << alpha << " to ";
                  alpha = alpha_lo + DataType(0.5)*(alpha_hi - alpha_lo);
                  //std::cout << alpha << std::endl;
                }
                width_old = width;
                width = Math::abs(alpha_lo - alpha_hi);
              }
            }

          }

          this->_alpha_min = alpha;//*this->_norm_dir;
          this->_fval_min = fval;

          if(status == Status::success)
            this->_alpha_0 = alpha*this->_norm_dir;
          else
          {
            this->_alpha_min = alpha_lo;//*this->_norm_dir;
            this->_fval_min = fval_lo;
            //std::cout << "Unusual termination alpha_lo " << alpha_lo << "fval_lo " << fval_lo << std::endl;
            vec_sol.axpy(vec_dir, this->_vec_initial_sol, alpha_lo);

            // Prepare and evaluate
            this->_op.prepare(vec_sol, this->_filter);

            // Compute and filter the gradient
            this->_op.compute_grad(this->_vec_grad);
            this->trim_func_grad(fval);
            this->_filter.filter_def(this->_vec_grad);
          }

          return status;
        }

        /**
         * \brief The great magick trick to find a minimum of a 1d function
         *
         * \param[out] alpha
         * The trial step.
         *
         * \param[in] fval
         * Functional value at alpha
         *
         * \param[in] df
         * Derivative value at alpha
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] fval_lo
         * Functional value at alpha_lo
         *
         * \param[in] df_lo
         * Derivative at alpha_lo
         *
         * \param[in] alpha_hi
         * The step size for the lower functional value
         *
         * \param[in] fval_hi
         * Functional value at alpha_hi
         *
         * \param[in] df_hi
         * Derivative at alpha_hi
         *
         * \param[in,out] min_in_interval
         * Do we know the minimum has to be in the interval (or any superset we worked on before)? Because this gets
         * called several times, this knowledge has to be passed back and forth.
         *
         * \param[in,out] drive_to_bndry
         * In the cases 1 and 3 this is set to true, in cases 2 and 4 to false. This determines if convergence can be
         * improved by driving the step size to the appropriate boundary more quickly.
         *
         * \returns
         * Status::success as there is nothing that can fail.
         *
         * This routine does a quadratic interpolation using fval_lo, df_lo and df, and a cubic interpolation using
         * fval_lo, df_lo, fval_hi and df_hi. The great magick is to determine which is better, and this depends on a
         * lot of factors, which leads to 4 distinct cases.
         *
         */
        Status _polynomial_fit(
          DataType& alpha, DataType& fval, DataType& df,
          DataType& alpha_lo, DataType& fval_lo, DataType& df_lo,
          DataType& alpha_hi, DataType& fval_hi, DataType& df_hi,
          bool& min_in_interval, bool& drive_to_bndry)
        {
          DataType alpha_c(0);
          DataType alpha_q(0);
          DataType alpha_new(0);

          // Case 1: The function value increased.
          if( fval > fval_lo)
          {
            // Computation of cubic step
            alpha_c  = _argmin_cubic(alpha_lo, fval_lo, df_lo, alpha, fval, df);
            // The first derivative was negative, but the function value increases, so the minimum has to be in this
            // interval (or have been in the original interval, which this interval is a subset of)
            min_in_interval = true;
            drive_to_bndry = true;
            // Use 2nd order approximation NOT using the derivative df because it might be garbage (i.e. if we hit a
            // singularity) because the function value increased
            alpha_q = _argmin_quadratic(alpha_lo, fval_lo, df_lo, alpha, fval, df, false);

            // If the cubic minimum is closer to the old best step length, take it. Otherwise, take the mean
            if(Math::abs(alpha_c - alpha_lo) < Math::abs(alpha_q - alpha_lo))
              alpha_new = alpha_c;
            else
              alpha_new = (alpha_q + alpha_c)/DataType(2);
            //std::cout << "Case 1: alpha " << alpha_new << " q " << alpha_q << " c " << alpha_c << std::endl;

          }
          // Case 2: The first derivative changes sign
          // We also know that the function value increased in the last trial step
          else if( df_lo*df < DataType(0) )
          {
            alpha_c  = _argmin_cubic(alpha, fval, df, alpha_lo, fval_lo, df_lo);
            // Because the derivative changed sign, it has to be zero somewhere in between, so the minimum is in the
            // current interval
            min_in_interval = true;
            drive_to_bndry = false;
            // Default: Quadratic interpolant using fval_lo, df_lo and df
            alpha_q = _argmin_quadratic(alpha_lo, fval_lo, df_lo, alpha, fval, df, true);
            // Take the step closer to the new step alpha
            Math::abs(alpha - alpha_c) > Math::abs(alpha - alpha_q) ? alpha_new = alpha_c : alpha_new = alpha_q;
            //std::cout << "Case 2: alpha " << alpha_new << " q " << alpha_q << " c " << alpha_c << std::endl;
          }
          // Case 3: The absolute value of the derivative increases
          // We also know that the function value increased in the last trial step and that the derivative did not
          // change sign.
          // This means we have to further search in the direction of alpha.
          else if(Math::abs(df) < Math::abs(df_lo))
          {
            alpha_c  = _argmin_cubic_case_3(alpha, fval, df, alpha_lo, fval_lo, df_lo);
            drive_to_bndry = true;

            // Quadratic interpolant using fval_lo, df_lo and df
            alpha_q = _argmin_quadratic(alpha, fval, df, alpha_lo, fval_lo, df_lo, true);

            // If the cubic step is closer to the trial step
            if(Math::abs(alpha - alpha_c) < Math::abs(alpha - alpha_q))
            {
              min_in_interval ? alpha_new = alpha_c : alpha_new = alpha_q;
            }
            else
            {
              min_in_interval ? alpha_new = alpha_q : alpha_new = alpha_c;
            }
            //std::cout << "Case 3: alpha " << alpha_new << " q " << alpha_q << " c " << alpha_c << std::endl;
          }
          // Case 4: The absolute value of the derivative did not increase
          // We also know that the function value increased in the last trial step and that the derivative did not
          // change sign.
          else
          {
            drive_to_bndry = false;
            // Note that the arguments for argmin_cubic are different in this case
            // Here we can just take the cubic step because the cubic interpolation sets it to the alpha_soft_max
            // if it recognises that the minimum is outside the current interval.
            if(min_in_interval)
            {
              alpha_c  = _argmin_cubic(alpha, fval , df, alpha_hi, fval_hi, df_hi);
              alpha_new = alpha_c;
            }
            else
            {
              alpha > alpha_lo ? alpha_new = _alpha_soft_max : alpha_new = _alpha_soft_min;
            }
            //std::cout << "Case 4: alpha " << alpha_new << std::endl;
          }

            // Update the inverval of uncertainity. Has to happen befor we clamp the step to the admissable interval.
            // Check which step(s) we need to replace
            if(fval > fval_lo)
            {
              alpha_hi = alpha;
              fval_hi = fval;
              df_hi = df;
            }
            else
            {
              if(df * df_lo < DataType(0))
              {
                alpha_hi = alpha_lo;
                fval_hi = fval_lo;
                df_hi = df_lo;
              }
              alpha_lo = alpha;
              fval_lo = fval;
              df_lo = df;
            }

            _clamp_step(alpha_new);

            // If we know the minimum is in the interval, we can drive alpha to the corresponding end more quickly
            if(min_in_interval && drive_to_bndry)
            {
              //std::cout << "Bracketing adjust " << alpha_lo << " " << alpha_hi << " " << alpha_new << std::endl;
              //std::cout << "Bracketing adjust alpha_new from " << alpha_new;
              if(alpha_lo < alpha_hi)
                alpha_new = Math::min(alpha_lo + DataType(0.66)*(alpha_hi - alpha_lo), alpha_new);
              else
                alpha_new = Math::max(alpha_lo + DataType(0.66)*(alpha_hi - alpha_lo), alpha_new);
              //std::cout << " to " << alpha_new << std::endl;
            }

            alpha = alpha_new;

            //std::cout << "Polynomial fit: new " << alpha << " " << fval << " " << df <<std::endl;
            //std::cout << "Polynomial fit: lo  " << alpha_lo << " " << fval_lo << " " << df_lo <<std::endl;
            //std::cout << "Polynomial fit: hi  " << alpha_hi << " " << fval_hi << " " << df_hi <<std::endl;

          return Status::success;
        }

        /**
         * \brief Enforces hard and soft step limits, adjusting the soft limits if neccessary
         *
         * \param[in,out] alpha_new
         * The step size to be clamped.
         *
         */
        void _clamp_step(DataType& alpha_new) const
        {
          alpha_new = Math::max(alpha_new, _alpha_hard_min);
          alpha_new = Math::min(alpha_new, _alpha_hard_max);

          alpha_new = Math::max(alpha_new, _alpha_soft_min);
          alpha_new = Math::min(alpha_new, _alpha_soft_max);
        }

        /**
         * \brief Computes the minimum of a quadratic interpolation polynomial
         *
         * This has two variants:
         *  a) Interpolate fval_lo, df_lo and fval_hi
         *  b) Interpolate fval_lo, df_lo and df_hi
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] alpha_hi
         * The step size for the higher functional value
         *
         * \param[in] fval_lo
         * Functional value at alpha_lo
         *
         * \param[in] fval_hi
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
        DataType _argmin_quadratic(const DataType alpha_lo, const DataType fval_lo, const DataType df_lo, const DataType alpha_hi, const DataType fval_hi, const DataType df_hi, const bool interpolate_derivative) const
        {
          DataType alpha(alpha_lo);

          // Quadratic interpolation using fval_lo, df_lo, df_hi
          if(interpolate_derivative)
            alpha += df_lo /(df_lo - df_hi) * (alpha_hi - alpha_lo);
          // Quadratic interpolation using fval_lo, fval_hi, df_lo
          else
            alpha += df_lo/( (fval_hi - fval_lo)/( alpha_hi - alpha_lo) - df_lo)/DataType(2)*(alpha_lo - alpha_hi);

          return alpha;
        }

        /**
         * \brief Computes the minimum of a cubic interpolation polynomial
         *
         * If the minimum of the cubic interpolation polynomial is in the interiour, it is computed here. If it is
         * not in the interiour, it lies across one of the endpoints depending on the functional and derivative
         * values. In this case, we set the return value to _alpha_soft_{min, max} accordingly.
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] alpha_hi
         * The step size for the higher functional value
         *
         * \param[in] fval_lo
         * Functional value at alpha_lo
         *
         * \param[in] fval_hi
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
         * \note The code is extremely ugly and sensitive to annihilation.
         *
         */
        DataType _argmin_cubic(DataType alpha_lo, DataType fval_lo, DataType df_lo,
        DataType alpha_hi, DataType fval_hi, DataType df_hi) const
        {
          DataType alpha(alpha_lo);

          DataType d1 = DataType(3)*(fval_lo - fval_hi)/(alpha_hi - alpha_lo) + df_lo + df_hi;
          // Scale the computation of r for better numerical stability
          DataType scale = Math::max(Math::abs(df_lo), Math::abs(df_hi));
          scale = Math::max(scale, Math::abs(d1));

          DataType r = Math::sqr(d1/scale) - (df_lo/scale) * (df_hi/scale);

          //std::cout << "fval_hi " << fval_hi << " fval_lo "<< fval_lo << std::endl;
          //std::cout << "alpha_hi " << alpha_hi << " alpha_lo "<< alpha_lo << std::endl;
          //std::cout << "df_hi " << df_hi << " df_lo "<< df_lo << std::endl;
          //std::cout << "d1 " << d1 << std::endl;
          //std::cout << "scale " << scale << std::endl;

          DataType d2(0);

          if(r > DataType(0))
          {
            d2 = Math::signum(alpha_hi - alpha_lo) * scale * Math::sqrt(r);

            DataType p(d2 - df_lo + d1);

            DataType q(0);
            // This sorting is done to avoid annihilation
            df_hi*df_lo > DataType(0) ?  q = d2 +(df_hi-df_lo) + d2 : q = d2 - df_lo + d2 + df_hi;

            //std::cout << "d2 " << d2 << std::endl;
            //std::cout << "p " << p << " q " << q << std::endl;

            alpha += (alpha_hi - alpha_lo)*(p/q);
          }
          // r <= 0 means that the minimum is not in the interior, so it has to lie across the endpoint with the
          // lower function value
          else
          {
            //d2 = Math::signum(alpha_hi - alpha_lo) * scale *
            //  Math::sqrt(Math::max(r, DataType(0)));
            //DataType q(0);
            //if(df_hi*df_lo > 0)
            //  q = d2 +(df_hi-df_lo) + d2;
            //else
            //  q = d2 - df_lo + d2 +df_hi;
            //std::cout << "d2 " << d2 << std::endl;
            //std::cout << "p " << (d2 - df_lo + d1) << " q " << q << std::endl;

            (alpha_lo < alpha_hi && df_lo > DataType(0) )? alpha = _alpha_soft_min: alpha = _alpha_soft_max;
          }

          return alpha;
        }

        /**
         * \copydoc _argmin_cubic
         * This is the weird check if the function value does not tend to infinity in the direction of the cubic
         * step.
         */
        DataType _argmin_cubic_case_3(
          const DataType alpha_lo, const DataType fval_lo, const DataType df_lo,
          const DataType alpha_hi, const DataType fval_hi, DataType df_hi) const
        {
          DataType alpha_c(alpha_lo);

          DataType d1 = DataType(3)*(fval_hi - fval_lo)/(alpha_lo - alpha_hi) + df_hi + df_lo;
          DataType scale = Math::max(Math::abs(df_lo), Math::abs(df_hi));
          scale = Math::max(scale, Math::abs(d1));

          DataType r = Math::sqr(d1/scale) - (df_lo/scale) * (df_hi/scale);

          DataType d2 = Math::signum(alpha_hi - alpha_lo) * scale * Math::sqrt(Math::max(r, DataType(0)));

          DataType p(d2 - df_lo + d1);
          DataType q(d2 +(df_hi - df_lo) + d2);

          //std::cout << "fval_hi " << fval_hi << " fval_lo "<< fval_lo << std::endl;
          //std::cout << "alpha_hi " << alpha_hi << " alpha_lo "<< alpha_lo << std::endl;
          //std::cout << "df_hi " << df_hi << " df_lo "<< df_lo << std::endl;
          //std::cout << "d1 " << d1 << std::endl;
          //std::cout << "scale " << scale << std::endl;
          //std::cout << "d2 " << d2 << std::endl;
          //std::cout << "p " << (d2 - df_lo + d1) << " q " << q << std::endl;

          if( p/q < DataType(0) && d2 != DataType(0))
            alpha_c += p/q*(alpha_hi - alpha_lo);
          else
          {
            //std::cout << "Weird ";
            alpha_lo > alpha_hi ? alpha_c = _alpha_soft_max : alpha_c = _alpha_soft_min;
          }

          return alpha_c;
        }

    }; // class StrongWolfeLinesearch

    /**
     * \brief Creates a new FixedStepLinesearch object
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] steplength
     * The step length the "line search" is to take
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new FixedStepLinesearch object.
     */
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<FixedStepLinesearch<Operator_, Filter_>> new_fixed_step_linesearch(
      Operator_& op, Filter_& filter, typename Operator_::DataType steplength, bool keep_iterates = false)
      {
        return std::make_shared<FixedStepLinesearch<Operator_, Filter_>>(op, filter, steplength, keep_iterates);
      }

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
     * \param[in] tol_decrease
     * Tolerance for sufficient decrease in function value.
     *
     * \param[in] tol_curvature
     * Tolerance for the curvature condition.
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
} // namespace FEAT
#endif // FEAT_KERNEL_SOLVER_LINESEARCH
