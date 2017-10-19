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
     * for a given search direction \f$ d \f$ by approximating the second order derivatives along \f$ d \f$ by a
     * secant.
     *
     */
    template<typename Functional_, typename Filter_>
    class SecantLinesearch : public Linesearch<Functional_, Filter_>
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
        /// Default initial step length
        static constexpr DataType secant_step_default = DataType(1e-2);

      protected:
        /// Default step for calculating the "other" secant point in the initial step. Crucial
        DataType _secant_step;

      public:
        /**
         * \brief Standard constructor
         *
         * \param[in, out] functional
         * The (nonlinear) functional. Cannot be const because it saves its own state.
         *
         * \param[in] filter
         * Filter to apply to the functional's gradient.
         *
         * \param[in] secant_step
         * Step length for setting the "other" secant point in the first iteration. Crucial.
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit SecantLinesearch(
          Functional_& functional, Filter_& filter,
          const DataType secant_step = secant_step_default,
          const bool keep_iterates = false) :
          BaseClass("S-LS", functional, filter, keep_iterates),
          _secant_step(secant_step)
          {
            this->_alpha_0 = _secant_step;
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
         * \param[in] functional
         * The functional.
         *
         * \param[in] filter
         * The system filter.
         *
         */
        explicit SecantLinesearch(const String& section_name, PropertyMap* section,
        Functional_& functional, Filter_& filter) :
          BaseClass("S-LS", section_name, section, functional, filter),
          _secant_step(secant_step_default)
          {
            auto secant_step_p = section->query("secant_step");
            if(secant_step_p.second)
            {
              set_secant_step(DataType(std::stod(secant_step_p.first)));
            }

            this->_alpha_0 = _secant_step;
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

        /// \copydoc BaseClass::reset()
        virtual void reset() override
        {
          BaseClass::reset();
          this->_alpha_0 = _secant_step;
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

        /// \copydoc SolverBase::write_config()
        virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name) const override
        {

          PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

          my_section->add_entry("secant_step", stringify_fp_sci(_secant_step));

          return my_section;
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

          // The step length wrt. to the NORMALISED search direction
          DataType alpha(0);
          // The functional value
          DataType fval(0);
          /// <vec_pn, vec_grad>. We want to find the minimum of the functional value along vec_pn
          DataType df(0);

          // Perform initialisations and checks
          Status st = this->_startup(alpha, fval, df, vec_sol, vec_dir);
          // Use the additional information about the preconditioned search direction's length?
          if(this->_dir_scaling)
          {
            alpha *= this->_norm_dir;
          }

          DataType alpha_update(alpha);

          // The second secant point in the first iteration is x + _secant_step * vec_dir
          // The first "other" point for the secant
          vec_sol.axpy(this->_vec_pn, this->_vec_initial_sol, alpha);

          // start iterating
          while(st == Status::progress)
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

            if(fval < this->_fval_min)
            {
              this->_fval_min = fval;
              this->_alpha_min = alpha;
            }

            // Set new defect, do convergence checks and compute the new _alpha
            DataType df_prev = df;
            df = this->_vec_pn.dot(this->_vec_grad);

            st = this->_check_convergence(fval, df, alpha);

            // Check if the diffence in etas is too small, thus leading to a huge update of relative size > sqrt(eps)
            if(Math::abs(df - df_prev) < Math::abs(alpha_update*df)*Math::sqrt(Math::eps<DataType>()))
            {
              st = Status::stagnated;
            }

            if(st != Status::progress)
            {
              break;
            }

            // Update alpha according to secant formula
            alpha_update *= df/(df_prev - df);
            alpha += alpha_update;

            // Update the solution
            vec_sol.axpy(this->_vec_pn, this->_vec_initial_sol, alpha);
          }

          // If we are successful, we could save the last step length as the new initial step length. This is
          // disabled by default, as it can lead to stagnation e.g. for minimising the Rosenbrock function.
          //if(st == Status::success)
          //{
          //  this->_alpha_0 = this->_alpha_min/this->_norm_dir;
          //}
          // If we are not successful, we update the best step length and need to re-evaluate everything for that
          // step
          if(st != Status::success)
          {
            vec_sol.axpy(this->_vec_pn, this->_vec_initial_sol, this->_alpha_min);

            // Prepare and evaluate
            this->_functional.prepare(vec_sol, this->_filter);
            this->_functional.eval_fval_grad(fval, this->_vec_grad);
            // DO NOT call trim() here, as this simple linesearch accepts bad steps too, which would then be trimmed
            // and used.
            //this->trim_func_grad(fval);
            this->_filter.filter_def(this->_vec_grad);
          }
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return st;
        }

    }; // class SecantLinesearch

    /**
     * \brief Creates a new SecantLinesearch object
     *
     * \param[in] functional
     * The functional.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] secant_step
     * Length for first secant step.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new SecantLinesearch object.
     */
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<SecantLinesearch<Functional_, Filter_>> new_secant_linesearch(
      Functional_& functional, Filter_& filter,
      typename Functional_::DataType secant_step = SecantLinesearch<Functional_, Filter_>::secant_step_default,
      bool keep_iterates = false)
      {
        return std::make_shared<SecantLinesearch<Functional_, Filter_>>(functional, filter, secant_step,
        keep_iterates);
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
     * \param[in] functional
     * The functional
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new SecantLinesearch object.
     */
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<SecantLinesearch<Functional_, Filter_>> new_secant_linesearch(
      const String& section_name, PropertyMap* section,
      Functional_& functional, Filter_& filter)
      {
        return std::make_shared<SecantLinesearch<Functional_, Filter_>>(section_name, section, functional, filter);
      }

  } // namespace Solver
} // namespace FEAT

#endif // FEAT_KERNEL_SOLVER_SECANT_LINESEARCH
