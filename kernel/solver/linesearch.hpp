#pragma once
#ifndef FEAT_KERNEL_SOLVER_LINESEARCH
#define FEAT_KERNEL_SOLVER_LINESEARCH 1
#include <kernel/base_header.hpp>
#include <kernel/solver/iterative.hpp>

#include <deque>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Linesearch base class
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
     * for a given search direction \f$ d \f$.
     *
     */
    template<typename Functional_, typename Filter_>
    class Linesearch : public IterativeSolver<typename Functional_::VectorTypeR>
    {
      public:
        /// Filter type to be applied to the gradient of the functional
        typedef Filter_ FilterType;
        /// Input vector type for the functional's gradient
        typedef typename Functional_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Functional_::DataType DataType;
        /// Our base class
        typedef IterativeSolver<typename Functional_::VectorTypeR> BaseClass;

      protected:
        /// The (nonlinear) functional
        // Note that this cannot be const, as the functional saves its state and thus changes
        Functional_& _functional;
        /// The filter to be applied to the functional's gradient
        Filter_& _filter;

        /// Gradient vector
        VectorType _vec_grad;
        /// Initial solution
        VectorType _vec_initial_sol;
        /// temporary vector
        VectorType _vec_tmp;

        /// Functional functional value
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
         * \param[in] plot_name
         * String to use in solver plots to console
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
        explicit Linesearch(const String& plot_name, Functional_& functional, Filter_& filter,
        bool keep_iterates = false) :
          BaseClass(plot_name),
          _functional(functional),
          _filter(filter),
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
            {
              iterates = new std::deque<VectorType>;
            }
          }

        /**
         * \brief Constructor using a PropertyMap
         *
         * \param[in] plot_name
         * String to use in solver plots to console
         *
         * \param[in] section_name
         * The name of the config section, which it does not know by itself
         *
         * \param[in] section
         * A pointer to the PropertyMap section configuring this solver
         *
         * \param[in, out] functional
         * The (nonlinear) functional. Cannot be const because it saves its own state
         *
         * \param[in] filter
         * Filter to apply to the functional's gradient
         *
         */
        explicit Linesearch(const String& plot_name, const String& section_name, PropertyMap* section,
        Functional_& functional, Filter_& filter) :
          BaseClass(plot_name, section_name, section),
          _functional(functional),
          _filter(filter),
          _fval_min(Math::huge<DataType>()),
          _fval_0(Math::huge<DataType>()),
          _trim_threshold(Math::huge<DataType>()),
          _alpha_min(DataType(0)),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          _tol_step(DataType(Math::pow(Math::eps<DataType>(), DataType(0.85)))),
          iterates(nullptr)
          {
            // Check if we have to keep the iterates
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second)
            {
              iterates = new std::deque<VectorType>;
            }

            // Check if we have to keep the iterates
            auto tol_step_p = section->query("tol_step");
            if(tol_step_p.second)
            {
              set_tol_step(DataType(std::stod(tol_step_p.first)));
            }
          }

        /// \copydoc ~BaseClass()
        virtual ~Linesearch()
        {
          if(iterates != nullptr)
          {
            delete iterates;
          }
        }

        /// \copydoc SolverBase::write_config()
        virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name = "") const override
        {
          XASSERT(parent != nullptr);

          Dist::Comm comm(Dist::Comm::world());

          PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

          my_section->add_entry("keep_iterates", stringify(iterates == nullptr ? 0 : 1));
          my_section->add_entry("tol_step", stringify_fp_sci(_tol_step));

          return my_section;

        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          // Create temporary vectors
          _vec_initial_sol = this->_functional.create_vector_r();
          _vec_tmp = this->_functional.create_vector_r();
          _vec_grad = this->_functional.create_vector_r();
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
          {
            iterates->clear();
          }
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
         * The initial functional value.
         *
         * This is handy because the linesearch gets called from another solver that in general already evaluated
         * the functional.
         *
         */
        void set_initial_fval(DataType f0)
        {
          _fval_0 = f0;
          if(_trim_threshold == Math::huge<DataType>())
          {
            _trim_threshold = DataType(10)*(Math::abs(_fval_0) + DataType(1));
          }
        }

        /**
         * \brief Sets the initial gradient from a defect vector
         *
         * \param[in] vec_def
         * The initial defect vector.
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
         * \brief Trims the function value and gradient according to some threshold
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

        /**
         * \brief Sets the step length tolerance
         *
         */
        void set_tol_step(DataType tol_step)
        {
          XASSERT(tol_step > DataType(0));

          _tol_step = tol_step;
        }
    }; // class Linesearch
  } // namespace Solver
} // namespace FEAT
#endif // FEAT_KERNEL_SOLVER_LINESEARCH
