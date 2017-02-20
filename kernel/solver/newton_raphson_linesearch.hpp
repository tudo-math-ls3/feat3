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
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

          Status status(Status::progress);
          this->_num_iter = Index(0);

          this->_vec_initial_sol.copy(sol);

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

          // Compute initial defect. We want to minimise d^T * grad(_functional)
          this->_def_init = Math::abs(dir.dot(this->_vec_grad));

          //sol.axpy(dir, this->_vec_initial_sol, alpha);
          // start iterating
          while(status == Status::progress)
          {
            IterationStats stat(*this);

            // Increase iteration count
            ++this->_num_iter;

            DataType fval(0);
            this->_functional.prepare(sol, this->_filter);
            this->_functional.eval_fval_grad(fval, this->_vec_grad);
            this->_filter.filter_def(this->_vec_grad);

            if(fval < this->_fval_min)
            {
              this->_fval_min = fval;
              this->_alpha_min = alpha;
            }

            // Update defect
            this->_def_cur = Math::abs(this->_vec_grad.dot(dir));

            // Compute new _alpha <- _alpha - grad.dot(dir) / dir.dot(Hess*dir)
            this->_functional.apply_hess(this->_vec_tmp, dir);
            this->_filter.filter_def(this->_vec_tmp);

            alpha_hidate = - this->_vec_grad.dot(dir)/dir.dot(this->_vec_tmp);
            alpha += alpha_hidate;

            //Statistics::add_solver_defect(this->_branch, double(this->_def_cur));
            sol.axpy(dir, this->_vec_initial_sol, alpha);

            // plot?
            if(this->_plot_iter())
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
              {
                Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
                return status;
              }
            }

          }

          // We should never come to this point
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return Status::undefined;
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
