#pragma once
#ifndef FEAT_KERNEL_SOLVER_NEWTON_RAPHSON_LINESEARCH
#define FEAT_KERNEL_SOLVER_NEWTON_RAPHSON_LINESEARCH 1
#include <kernel/base_header.hpp>
#include <kernel/solver/linesearch.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Newton Raphson linesearch
     *
     * \tparam Functional_
     * The (nonlinear) functional to be evaluated
     *
     * \tparam Filter_
     * The filter to be applied to the functional's gradient
     *
     * This class implements a linesearch which approximately finds
     * \f[
     *   \alpha^* = \mathrm{argmin} \nabla f(x + \alpha d) \cdot d
     * \f]
     * for a given search direction \f$ d \f$ by applying a Newton Raphson iteration to this.
     *
     */
    template<typename Functional_, typename Filter_>
    class NewtonRaphsonLinesearch : public Linesearch<Functional_, Filter_>
    {
      public:
        /// Filter type to be applied to the gradient of the functional
        typedef Filter_ FilterType;
        /// Input vector type for the functional's gradient
        typedef typename Functional_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Functional_::DataType DataType;
        /// Our base class
        typedef Linesearch<Functional_, Filter_> BaseClass;

      public:
        /**
         * \brief Standard constructor
         *
         * \param[in, out] functional
         * The (nonlinear) functional. Cannot be const because it saves its own state
         *
         * \param[in] filter
         * Filter to apply to the functional's gradient
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit NewtonRaphsonLinesearch(Functional_& functional, Filter_& filter, bool keep_iterates = false) :
          BaseClass("NR-LS", functional, filter, keep_iterates)
          {
          }

        /**
         * \brief Constructor using a PropertyMap
         *
         * \param[in] section_name
         * The name of the config section, which it does not know by itself.
         *
         * \param[in] section
         * A pointer to the PropertyMap section configuring this solver.
         *
         * \param[in] functional
         * The functional.
         *
         * \param[in] filter
         * The system filter.
         *
         * \returns
         * A shared pointer to a new NewtonRaphsonLinesearch object.
         */
        explicit NewtonRaphsonLinesearch(const String& section_name, PropertyMap* section,
        Functional_& functional, Filter_& filter) :
          BaseClass("NR-LS", section_name, section, functional, filter)
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
          this->_functional.prepare(vec_cor, this->_filter);

          // apply
          Status st(_apply_intern(vec_cor, vec_dir));
          this->plot_summary(st);
          return st;
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
          this->_functional.prepare(vec_sol, this->_filter);
          // apply
          Status st(_apply_intern(vec_sol, vec_dir));
          this->plot_summary(st);
          return st;
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
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

          Status status(Status::progress);

          DataType alpha(0);
          // It is critical that this and _vec_grad were set from the outside!
          DataType fval(this->_fval_0);
          DataType df(vec_dir.dot(this->_vec_grad));

          // Setup
          status = this->_startup(fval, df, vec_sol, vec_dir);

          // Compute new _alpha <- _alpha - grad.dot(vec_dir) / vec_dir.dot(Hess*vec_dir)
          this->_functional.apply_hess(this->_vec_tmp, vec_dir);
          this->_filter.filter_def(this->_vec_tmp);

          alpha -= this->_vec_grad.dot(vec_dir)/vec_dir.dot(this->_vec_tmp);

          // start iterating
          while(status == Status::progress)
          {
            IterationStats stat(*this);

            // Increase iteration count
            ++this->_num_iter;

            this->_functional.prepare(vec_sol, this->_filter);
            this->_functional.eval_fval_grad(fval, this->_vec_grad);
            // DO NOT call trim() here, as this simple linesearch accepts bad steps too, which would then be trimmed
            // and used.
            //this->trim_func_grad(fval);
            this->_filter.filter_def(this->_vec_grad);

            df = vec_dir.dot(this->_vec_grad);

            if(fval < this->_fval_min)
            {
              this->_fval_min = fval;
              this->_alpha_min = alpha;
            }

            status = this->_check_convergence(fval, df, alpha);

            if(status != Status::progress)
            {
              break;
            }

            // Compute new _alpha <- _alpha - grad.dot(vec_dir) / vec_dir.dot(Hess*vec_dir)
            this->_functional.apply_hess(this->_vec_tmp, vec_dir);
            this->_filter.filter_def(this->_vec_tmp);

            alpha -= this->_vec_grad.dot(vec_dir)/vec_dir.dot(this->_vec_tmp);

            // Update solution
            vec_sol.axpy(vec_dir, this->_vec_initial_sol, alpha);

          }

          // If we are not successful, we update the best step length and need to re-evaluate everything for that
          // step
          if(status != Status::success)
          {
            vec_sol.axpy(vec_dir, this->_vec_initial_sol, this->_alpha_min);

            // Prepare and evaluate
            this->_functional.prepare(vec_sol, this->_filter);
            this->_functional.eval_fval_grad(fval, this->_vec_grad);
            // DO NOT call trim() here, as this simple linesearch accepts bad steps too, which would then be trimmed
            // and used.
            this->trim_func_grad(fval);
            this->_filter.filter_def(this->_vec_grad);
          }

          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return status;
        }

    }; // class NewtonRaphsonLinesearch

    /**
     * \brief Creates a new NewtonRaphsonLinesearch object
     *
     * \param[in] functional
     * The functional.
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
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<NewtonRaphsonLinesearch<Functional_, Filter_>> new_newton_raphson_linesearch(
      Functional_& functional, Filter_& filter, bool keep_iterates = false)
      {
        return std::make_shared<NewtonRaphsonLinesearch<Functional_, Filter_>>(functional, filter, keep_iterates);
      }

    /**
     * \brief Creates a new NewtonRaphsonLinesearch object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] functional
     * The functional
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new NewtonRaphsonLinesearch object.
     */
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<NewtonRaphsonLinesearch<Functional_, Filter_>> new_newton_raphson_linesearch(
      const String& section_name, PropertyMap* section,
      Functional_& functional, Filter_& filter)
      {
        return std::make_shared<NewtonRaphsonLinesearch<Functional_, Filter_>>(section_name, section, functional, filter);
      }
  } // namespace Solver
} // namespace FEAT

#endif // FEAT_KERNEL_SOLVER_NEWTON_RAPHSON_LINESEARCH
