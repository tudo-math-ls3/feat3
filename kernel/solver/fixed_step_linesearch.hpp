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
        /// The length of the step
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
            this->_alpha_min = _steplength;
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
         * \param[in] op
         * The operator
         *
         * \param[in] filter
         * The system filter.
         *
         * \returns
         * A shared pointer to a new FixedStepLinesearch object.
         */
        explicit FixedStepLinesearch(const String& section_name, PropertyMap* section,
        Operator_& op_, Filter_& filter_) :
          BaseClass("FS-LS", section_name, section, op_, filter_),
          _steplength(-1)
          {
            this->set_max_iter(0);
            this->_alpha_min = _steplength;

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
        virtual DataType get_rel_update() override
        {
          return _steplength;
        }

        /**
         * \brief Sets the step length
         */
        virtual void set_step_length(DataType step_length)
        {
          XASSERT(step_length > DataType(0));

          _steplength = step_length;
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
          this->_op.eval_fval_grad(this->_fval_min, this->_vec_grad);
          this->_filter.filter_def(this->_vec_grad);

          return Status::success;
        }
    }; // class FixedStepLinesearch

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
     * \brief Creates a new FixedStepLinesearch object using a PropertyMap
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
     * A shared pointer to a new FixedStepLinesearch object.
     */
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<FixedStepLinesearch<Operator_, Filter_>> new_fixed_step_linesearch(
      const String& section_name, PropertyMap* section, Operator_& op, Filter_& filter)
      {
        return std::make_shared<FixedStepLinesearch<Operator_, Filter_>>(section_name, section, op, filter);
      }
  } // namespace Solver
} // namespace FEAT
#endif // FEAT_KERNEL_SOLVER_FIXED_STEP_LINESEARCH
