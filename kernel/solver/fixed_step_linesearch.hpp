#pragma once
#ifndef FEAT_KERNEL_SOLVER_FIXED_STEP_LINESEARCH
#define FEAT_KERNEL_SOLVER_FIXED_STEP_LINESEARCH 1
#include <kernel/base_header.hpp>
#include <kernel/solver/linesearch.hpp>

namespace FEAT
{
  namespace Solver
  {

    /**
     * \brief Fixed step line search
     *
     * This does not perform a line search, it just returns a constant step size
     *
     */
    template<typename Functional_, typename Filter_>
    class FixedStepLinesearch : public Linesearch<Functional_, Filter_>
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

      private:
        /// The length of the step
        DataType _step_length;

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
         * \param[in] step_length
         * The length of the step returned by this "line search".
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit FixedStepLinesearch(Functional_& functional, Filter_& filter, const DataType step_length,
        bool keep_iterates = false) :
          BaseClass("FS-LS", functional, filter, keep_iterates),
          _step_length(step_length)
          {
            this->_alpha_min = _step_length;
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
         * A shared pointer to a new FixedStepLinesearch object.
         */
        explicit FixedStepLinesearch(const String& section_name, PropertyMap* section,
        Functional_& functional, Filter_& filter) :
          BaseClass("FS-LS", section_name, section, functional, filter),
          _step_length(-1)
          {
            this->set_max_iter(0);
            this->_alpha_min = _step_length;

            // Check if we need to set the step length
            auto step_length_p = section->query("step_length");
            if(step_length_p.second)
            {
              set_step_length(std::stod(step_length_p.first));
            }
            else
            {
              throw InternalError(__func__,__FILE__,__LINE__,
              name()+" config section is missing the mandatory step_length key!");
            }

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
        virtual DataType get_rel_update() const override
        {
          return _step_length*this->_norm_dir;
        }

        /**
         * \brief Sets the step length
         */
        virtual void set_step_length(DataType step_length)
        {
          XASSERT(step_length > DataType(0));

          _step_length = step_length;
        }

        /// \copydoc SolverBase::write_config()
        virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name) const override
        {

          PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

          my_section->add_entry("step_length", stringify_fp_sci(_step_length));

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
          this->_norm_dir = vec_dir.norm2();
          // clear solution vector
          vec_cor.format();
          vec_cor.axpy(vec_dir, vec_cor, _step_length);

          this->_functional.prepare(vec_cor, this->_filter);

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
          this->_norm_dir = vec_dir.norm2();
          vec_sol.axpy(vec_dir, vec_sol, _step_length);

          this->_functional.prepare(vec_sol, this->_filter);
          this->_functional.eval_fval_grad(this->_fval_min, this->_vec_grad);
          this->_filter.filter_def(this->_vec_grad);

          return Status::success;
        }
    }; // class FixedStepLinesearch

    /**
     * \brief Creates a new FixedStepLinesearch object
     *
     * \param[in] functional
     * The functional.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] step_length
     * The step length the "line search" is to take.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false.
     *
     * \returns
     * A shared pointer to a new FixedStepLinesearch object.
     */
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<FixedStepLinesearch<Functional_, Filter_>> new_fixed_step_linesearch(
      Functional_& functional, Filter_& filter, typename Functional_::DataType step_length, bool keep_iterates = false)
      {
        return
          std::make_shared<FixedStepLinesearch<Functional_, Filter_>>(functional, filter, step_length, keep_iterates);
      }

    /**
     * \brief Creates a new FixedStepLinesearch object using a PropertyMap
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
     * A shared pointer to a new FixedStepLinesearch object.
     */
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<FixedStepLinesearch<Functional_, Filter_>> new_fixed_step_linesearch(
      const String& section_name, PropertyMap* section, Functional_& functional, Filter_& filter)
      {
        return
          std::make_shared<FixedStepLinesearch<Functional_, Filter_>>(section_name, section, functional, filter);
      }
  } // namespace Solver
} // namespace FEAT
#endif // FEAT_KERNEL_SOLVER_FIXED_STEP_LINESEARCH
