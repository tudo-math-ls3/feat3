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
          _op.compute_grad(this->_vec_grad);
          _filter.filter_def(this->_vec_grad);

          // Compute initial defect. We want to minimise d^T * grad(_op)
          this->_def_init = Math::abs(dir.dot(this->_vec_grad));

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            this->_op.prepare(sol);
            _op.compute_grad(this->_vec_grad);
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

          // We should never come to this point
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
            this->_min_iter = Index(1);

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

          _op.prepare(sol);
          _op.compute_grad(this->_vec_tmp);
          _filter.filter_def(this->_vec_tmp);
          _eta = dir.dot(this->_vec_tmp);
          if(_eta > DataType(0))
            throw InternalError("Search direction is not a descent direction: "+stringify_fp_sci(_eta));

          // The first "other" point for the secant
          this->_vec_tmp.axpy(dir, sol, _sigma_0);
          // Set the state in the operator
          _op.prepare(this->_vec_tmp);
          _op.compute_grad(this->_vec_tmp);
          _filter.filter_def(this->_vec_tmp);

          _eta = dir.dot(this->_vec_tmp);
          this->_def_init = Math::abs(_eta);

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            this->_op.prepare(sol);
            _op.compute_grad(this->_vec_tmp);
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
    class SWLinesearch : public IterativeSolver<typename Operator_::VectorTypeR>
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

        VectorType _vec_initial_sol;
        /// Gradient vector
        VectorType _vec_grad;
        /// temporary vector
        VectorType _vec_tmp;

        /// Line search parameter
        DataType _alpha;
        DataType _alpha_0;
        DataType _alpha_max;
        DataType _alpha_min;
        DataType _alpha_prev;
        DataType _delta;
        DataType _delta_0;
        DataType _delta_prev;
        DataType _f;
        DataType _f_0;
        DataType _f_prev;
        /// The 2-norm of the search direction
        DataType _norm_dir;
        /// The 2-norm of the iterate
        DataType _norm_sol;
        DataType _tol_decrease;
        DataType _tol_curvature;

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
        explicit SWLinesearch(Operator_& op_, const Filter_& filter_, bool keep_iterates = false) :
          BaseClass("SW-LS"),
          _op(op_),
          _filter(filter_),
          _alpha(DataType(0)),
          _alpha_0(DataType(1)),
          _alpha_max(DataType(0)),
          _alpha_min(DataType(0)),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          _tol_decrease(DataType(1e-4)),
          _tol_curvature(DataType(0.9)),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;

            this->_min_stag_iter = Index(2);
            this->_min_iter = Index(1);
            //this->_max_iter = Index(4);
          }

        /// \copydoc ~BaseClass()
        virtual ~SWLinesearch()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          _vec_initial_sol = this->_op.create_vector_r();
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
          return "Strong Wolfe Linesearch";
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
          auto normed_dir = vec_dir.clone();
          normed_dir.scale(normed_dir, DataType(1)/vec_dir.norm2());
          std::cout << "normed_dir = " << normed_dir << std::endl;

          // apply
          Status st =_apply_intern(vec_sol, normed_dir);

          return st;
        }

        virtual DataType get_rel_update()
        {
          return DataType(1);
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
        virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_dir)
        {
          static constexpr DataType extrapolation_width = DataType(4);
          Status status(Status::progress);
          this->_num_iter = Index(0);

          _vec_initial_sol.clone(vec_sol);

          // Norm of the search direction vector
          _norm_dir = vec_dir.norm2();
          // Norm of the initial guess
          _norm_sol = vec_sol.norm2();

          //DataType alpha_max(DataType(1)/_norm_dir);

          // Set state in the operator
          _op.prepare(vec_sol);

          _f_0 = _op.compute_func();
          _f = _f_0;
          _f_prev = _f_0;

          // Compute and filter gradient
          _op.compute_grad(this->_vec_grad);
          _filter.filter_def(this->_vec_grad);

          // Compute initial defect. We want to minimise d^T * grad(_op)
          _delta_0 = vec_dir.dot(this->_vec_grad)/_norm_dir;
          this->_def_init = Math::abs(_delta_0);

          std::cout << "SWLinesearch start:     sol = " << _vec_initial_sol << " dir = " << vec_dir << " grad = "
          << this->_vec_grad << std::endl;
          std::cout << "SWLinesearch start:     f_0 = " << stringify_fp_sci(_f_0) << " delta_0 = "
          << stringify_fp_sci(_delta_0) << std::endl;

          _alpha_prev = DataType(0);

          _alpha = _alpha_0;
          _alpha = DataType(0.464349);
          vec_sol.axpy(vec_dir, _vec_initial_sol, _alpha/_norm_dir);

          _alpha_min = DataType(0);
          if(_alpha_max < Math::eps<DataType>())
            _alpha_max = Math::huge<DataType>();

          DataType alpha_soft_max = _alpha_min + extrapolation_width*_alpha;
          alpha_soft_max = Math::min(alpha_soft_max, _alpha_max);

          _delta = _delta_0;
          _delta_prev = _delta_0;

          std::cout << "SWLinesearch start: alpha = " << stringify_fp_sci(_alpha) << " soft max = "
          << stringify_fp_sci(alpha_soft_max) << std::endl;

          if(_delta_0 > DataType(0))
          {
            throw InternalError("Initial search direction is not a descent direction: "+stringify_fp_sci(_delta_0));
          }


          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            // Save old data
            //DataType def_prev(this->_def_cur);
            _delta_prev = _delta;
            _f_prev = _f;

            // Prepare and evaluate
            this->_op.prepare(vec_sol);
            _f = _op.compute_func();

            // Increase iteration count
            ++this->_num_iter;

            // First let's see if we have to unconditionally compute the defect
            bool calc_def = false;
            calc_def = calc_def || (this->_min_iter < this->_max_iter);
            calc_def = calc_def || this->_plot;
            calc_def = calc_def || (this->_min_stag_iter > Index(0));


            // Check if the sufficient decrease condition is violated. If so, one call to _backet gets the solution
            if( (_f > _f_0 + _tol_decrease*_alpha*_delta_0) || (_f > _f_prev && this->_num_iter > 1) )
            {
              std::cout << "SWLinesearch: sufficient decrease violated" << std::endl;
              if(_f > _f_0 + _tol_decrease*_alpha*_delta_0)
                std::cout << "  f = " << stringify_fp_sci(_f) << " > " << stringify_fp_sci(_f_0 + _tol_decrease*_alpha*_delta_0) << " = _f_0 + _tol_decrease*_alpha*_delta_0)" << std::endl;
              if(_f > _f_prev)
                std::cout << "  f = " << stringify_fp_sci(_f) << " > " << stringify_fp_sci(_f_prev ) << " = f_prev" << std::endl;
                status = _bracket(vec_sol, vec_dir, _alpha_prev, _alpha, _delta_prev, _delta, _f_prev, _f);

                _alpha = _alpha_prev;
                _delta = _delta_prev;
                _f = _f_prev;

                _alpha_0 = _alpha;

                // Bracket already updated the solution and the defect etc.
                calc_def = false;

                this->_def_cur = Math::abs(_delta);
                Statistics::add_solver_defect(this->_branch, double(this->_def_cur));
            }

            if(calc_def)
            {
              _op.compute_grad(this->_vec_grad);
              _filter.filter_def(this->_vec_grad);
              _delta = this->_vec_grad.dot(vec_dir)/_norm_dir;

              this->_def_cur = Math::abs(_delta);
              Statistics::add_solver_defect(this->_branch, double(this->_def_cur));
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

            if(status != Status::progress)
              return status;

            // ensure that the defect is neither NaN nor infinity
            if(!Math::isfinite(this->_def_cur))
              return Status::aborted;

            // is diverged?
            if(this->is_diverged())
              return Status::diverged;

            //// minimum number of iterations performed?
            //if(this->_num_iter < this->_min_iter)
            //  return Status::progress;

            // If we did not compute the defect before, we have to do it now
            if(!calc_def)
            {
              _op.compute_grad(this->_vec_grad);
              _filter.filter_def(this->_vec_grad);
              _delta = this->_vec_grad.dot(vec_dir)/_norm_dir;

              this->_def_cur = Math::abs(_delta);
              Statistics::add_solver_defect(this->_branch, double(this->_def_cur));
            }

            // If we come to here, the sufficient decrease condition is satisfied. If the curvature condition is
            // satisfied, the current _alpha is good and we are successfull
            if( this->_def_cur < -_tol_curvature*_delta_0)
            {
              std::cout << "SWLinesearch: Sufficient decrease condition fulfilled" << std::endl;
              std::cout << "   f = " << stringify_fp_sci(_f) << " < " << stringify_fp_sci(_f_0 + _tol_decrease*_alpha*_delta_0) << std:: endl;
              std::cout << "   and curvature condition satisfied: " << stringify_fp_sci(this->_def_cur)
              << " < " << stringify_fp_sci(-_tol_curvature*_delta_0) << std::endl;

              //status = _bracket(vec_sol, vec_dir, _alpha, _alpha_prev, _delta, _delta_prev, _f, _f_prev);
              _alpha_0 = _alpha;
              return Status::success;
            }
            if(_delta >= DataType(0))
            {
              std::cout << "SWLinesearch: derivative is positive! " << stringify_fp_sci(_delta) << std::endl;
              status = _bracket(vec_sol, vec_dir, _alpha, _alpha_prev, _delta, _delta_prev, _f, _f_prev);
              _alpha_0 = _alpha;
              return status;
            }

            //_alpha_max = _alpha_prev + DataType(4)*(_alpha - _alpha_prev);
            _alpha_prev = _alpha;
            _delta_prev = _delta;
            //_alpha = _argmin_quadratic(_alpha, _alpha_max, _f_prev, _f, _delta_prev);
            _alpha = _alpha + DataType(0.33)*(alpha_soft_max - _alpha_min);

            _alpha = Math::max(_alpha, _alpha_min);
            _alpha = Math::min(_alpha, _alpha_max);

            if(_alpha > alpha_soft_max)
            {
              std::cout << "SWLinesearch: increase alpha_min from " << stringify_fp_sci(_alpha_min) << " to "
                << stringify_fp_sci(_alpha) << std::endl;
              std::cout << "SWLinesearch: increase soft max from  " << stringify_fp_sci(alpha_soft_max) << " to ";
              alpha_soft_max = _alpha + extrapolation_width*(_alpha - _alpha_min);
              _alpha_min = _alpha;
              std::cout << stringify_fp_sci(alpha_soft_max) << std::endl;
            }

            std::cout << "SWLinesearch: alpha = " << stringify_fp_sci(_alpha) << " soft max "
            << stringify_fp_sci(alpha_soft_max) << std::endl;

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

            // Update solution: sol <- sol + _alpha*dir
            vec_sol.axpy(vec_dir, _vec_initial_sol, _alpha/_norm_dir);

            if(this->_num_iter > this->_max_iter)
              return Status::max_iter;

          }

          // We should never come to this point
          return Status::undefined;

        }


        Status _bracket(VectorType& vec_sol, const VectorType& vec_dir, DataType& alpha_lo, DataType& alpha_hi, DataType& delta_lo, DataType& delta_hi, DataType& f_lo, DataType& f_hi)
        {
          Status st(Status::progress);

          DataType f_new(Math::Limits<DataType>::max());
          int iter(0);
          while(st == Status::progress && iter < 100)
          {
            ++iter;
          std::cout << "Bracket iter " << iter << " sol = " << vec_sol << " dir = " << vec_dir
            << " grad = " << this->_vec_grad << std::endl;
          std::cout << "  alpha_lo " << stringify_fp_sci(alpha_lo) << " f_lo " << stringify_fp_sci(f_lo) << " delta_lo " << stringify_fp_sci(delta_lo) << std::endl;
          std::cout << "  alpha_hi " << stringify_fp_sci(alpha_hi) << " f_hi " << stringify_fp_sci(f_hi) << " delta_hi " << stringify_fp_sci(delta_hi) << std::endl;

            DataType alpha_new(DataType(0));

            _find_in_bracket(alpha_new, alpha_lo, alpha_hi, f_lo, f_hi, delta_lo);

            vec_sol.axpy(vec_dir, _vec_initial_sol, alpha_new/_norm_dir);

            this->_op.prepare(vec_sol);
            f_new = _op.compute_func();

            // Update of intervall of uncertainity
            if(f_new > this->_f_0 +_tol_decrease*alpha_new*_delta_0 || f_new > f_lo)
            {
              std::cout << "  sufficient decrease violated" << std::endl;
              std::cout << "  f_0 + tol_decrease*alpha_new*_delta_0 " << stringify_fp_sci(this->_f_0 +_tol_decrease*alpha_new*_delta_0)
                << " f_lo " << stringify_fp_sci(f_lo) << "f_new " << stringify_fp_sci(f_new) << std::endl;
              alpha_hi = alpha_new;
              f_hi = f_new;
            }
            else
            {
              std::cout << "  bracket: sufficient decrease satisfied:  f = " << stringify_fp_sci(f_new) << " <= "
              << stringify_fp_sci(this->_f_0 +_tol_decrease*alpha_new*_delta_0) << std::endl;

              _op.compute_grad(this->_vec_grad);
              DataType delta_new = this->_vec_grad.dot(vec_dir);

              if(Math::abs(delta_new) < -_tol_curvature*_delta_0)
              {
                std::cout << "           curvature condition satisfied " << stringify_fp_sci(Math::abs(delta_new))
                << " < "  << stringify_fp_sci(-_tol_curvature*_delta_0) << std::endl;

                alpha_lo = alpha_new;
                f_lo = f_new;
                delta_lo = delta_new;

                st = Status::success;
              }
              else
              {
                std::cout << "  curvature condition violated: " << stringify_fp_sci(Math::abs(delta_new)) << " > "
                << stringify_fp_sci(-_tol_curvature*_delta_0) << std::endl;

                if(delta_new*(alpha_hi - alpha_lo) >= DataType(0))
                {
                  alpha_hi = alpha_lo;
                  f_hi = f_lo;
                }
              }

              alpha_lo = alpha_new;
              f_lo = f_new;
              delta_lo = delta_new;

            }

          }
          _alpha_0 = alpha_lo;

          if(st != Status::success)
          {
            vec_sol.axpy(vec_dir, _vec_initial_sol, alpha_lo/_norm_dir);

            this->_op.prepare(vec_sol);
            _op.compute_grad(this->_vec_grad);
          }

          std::cout << "Bracket stops at sol = " << vec_sol << " grad = " << this->_vec_grad << std::endl;
          std::cout << "  alpha_lo = " << stringify_fp_sci(alpha_lo) << " f_lo = "<< stringify_fp_sci(f_lo)
          << " delta_lo = " << stringify_fp_sci(delta_lo) << std::endl;
          std::cout << "  alpha_hi = " << stringify_fp_sci(alpha_hi) << " f_hi = "<< stringify_fp_sci(f_hi)
          << " delta_hi = " << stringify_fp_sci(delta_hi) << std::endl;

          return st;
        }

        Status _find_in_bracket(DataType& alpha_new, DataType alpha_lo, DataType alpha_hi, DataType f_lo, DataType f_hi, DataType delta_lo)
        {


          DataType alpha_q = _argmin_quadratic(alpha_lo, alpha_hi, f_lo, f_hi, delta_lo);
          //DataType alpha_c = _argmin_cubic(alpha_lo, alpha_hi, f_lo, f_hi, delta_lo, delta_hi);
          //
          alpha_new = alpha_q;
          return Status::success;
        }

        DataType _argmin_quadratic(DataType alpha_lo, DataType alpha_hi, DataType f_0, DataType f_1, DataType df_0)
        {
          DataType alpha(alpha_lo);

          alpha -= df_0*(alpha_hi - alpha_lo)/(DataType(2)*( (f_1 - f_0)/(alpha_hi - alpha_lo) - df_0 ) );
          std::cout << "Argmin_quadratic: alpha_new = " <<stringify_fp_sci(alpha) << ", input alpha = [" << stringify_fp_sci(alpha_lo) << " / " << stringify_fp_sci(alpha_hi) << "]" << std::endl;
          std::cout << "  f = " << stringify_fp_sci(f_0) << " / " << stringify_fp_sci(f_1) << ", df = " << stringify_fp_sci(df_0) << std::endl;

          return alpha;
        }

        DataType _argmin_cubic(DataType alpha_lo, DataType alpha_hi, DataType f_0, DataType f_1, DataType df_lo, DataType df_hi)
        {
          DataType alpha(alpha_lo);

          DataType d1 = (df_hi + df_lo) + DataType(3)*(f_hi - f_lo)/(alpha_hi - alpha_lo);

          DataType scale = Math::max(Math::abs(df_lo), Math::abs(df_hi));
          scale = Math::max(scale, Math::abs(d1));

          DataType d2 = Math::sign(alpha_lo - alpha_hi) * scale Math::sqrt( Math::sqr(d1/scale) - df_lo/scale * df_hi/scale);

          alpha += (alpha_hi - alpha_lo)*(-d2 - df_lo + d1)/(d2 + df_hi - df_lo + d2);
          return alpha;
        }
    }; // class SWLinesearch

  } // namespace Solver
} // namespace FEAST
#endif // FEAST_KERNEL_SOLVER_LINESEARCH
