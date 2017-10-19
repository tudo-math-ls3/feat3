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
     *   \alpha^* = \mathrm{argmin} \left< \nabla f(x + \alpha d), d \right>
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
        /// descend direction vector, normalised for better numerical stability
        VectorType _vec_pn;

        /// Functional functional value
        DataType _fval_min;
        /// Initial functional value
        DataType _fval_0;
        /// Threshold for trimming function value and gradient
        DataType _trim_threshold;

        /// Initial line search parameter
        DataType _alpha_0;
        /// Line search parameter
        DataType _alpha_min;
        /// Initial <vec_dir, vec_grad>
        DataType _delta_0;
        /// The 2-norm of the search direction
        DataType _norm_dir;
        /// The 2-norm of the iterate
        DataType _norm_sol;

        /// Tolerance for sufficient decrease in the norm of the gradient (Wolfe conditions)
        DataType _tol_curvature;
        /// Tolerance for sufficient decrease in the functional value (Wolfe conditions)
        DataType _tol_decrease;
        /// Tolerance for the update step
        DataType _tol_step;
        /// Use the search direction norm for step scaling? This can be important if the search direction was
        /// preconditioned (e.g. the Newton direction)
        bool _dir_scaling;

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
          _alpha_0(1),
          _alpha_min(0),
          _delta_0(Math::huge<DataType>()),
          _norm_dir(0),
          _norm_sol(0),
          _tol_curvature(DataType(0.3)),
          _tol_decrease(DataType(1e-3)),
          _tol_step(DataType(Math::pow(Math::eps<DataType>(), DataType(0.85)))),
          _dir_scaling(false),
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
          _alpha_0(1),
          _alpha_min(DataType(0)),
          _delta_0(Math::huge<DataType>()),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          _tol_curvature(0.3),
          _tol_decrease(1e-3),
          _tol_step(DataType(Math::pow(Math::eps<DataType>(), DataType(0.85)))),
          _dir_scaling(false),
          iterates(nullptr)
          {
            // Check if we have to keep the iterates
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second)
            {
              iterates = new std::deque<VectorType>;
            }

            auto tol_curvature_p = section->query("tol_curvature");
            if(tol_curvature_p.second)
            {
              set_tol_curvature(DataType(std::stod(tol_curvature_p.first)));
            }

            auto tol_decrease_p = section->query("tol_decrease");
            if(tol_decrease_p.second)
            {
              set_tol_decrease(DataType(std::stod(tol_decrease_p.first)));
            }
            // Check if we have to keep the iterates
            auto tol_step_p = section->query("tol_step");
            if(tol_step_p.second)
            {
              set_tol_step(DataType(std::stod(tol_step_p.first)));
            }
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~Linesearch()
        {
          if(iterates != nullptr)
          {
            delete iterates;
          }
        }

        /**
         * \brief Plot a summary of the last solver run
         */
        virtual void plot_summary(const Status st) const override
        {
          // Print solver summary
          if(this->_plot_summary())
          {
            Dist::Comm comm_world(Dist::Comm::world());

            String msg(this->get_plot_name()+ ": its: "+stringify(this->get_num_iter())+" ("+ stringify(st)+")"
                +", evals: "+stringify(_functional.get_num_func_evals())+" (func) "
                + stringify(_functional.get_num_grad_evals()) + " (grad) "
                + stringify(_functional.get_num_hess_evals()) + " (hess)"
                +" last step: "+stringify_fp_sci(_alpha_min)+"\n");
            msg +=this->get_plot_name()+": fval: "+stringify_fp_sci(_fval_0)
              + " -> "+stringify_fp_sci(_fval_min)
              + ", factor "+stringify_fp_sci(_fval_min/_fval_0)+"\n";
            msg += this->get_plot_name()  +": <dir, grad>: "+stringify_fp_sci(this->_def_init)
              + " -> "+stringify_fp_sci(this->_def_cur)
              + ", factor " +stringify_fp_sci(this->_def_cur/this->_def_init);
            comm_world.print(msg);
          }
        }

        /// \copydoc SolverBase::write_config()
        virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name = "") const override
        {
          XASSERT(parent != nullptr);

          Dist::Comm comm(Dist::Comm::world());

          PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

          my_section->add_entry("keep_iterates", stringify(iterates == nullptr ? 0 : 1));
          my_section->add_entry("tol_decrease", stringify_fp_sci(_tol_decrease));
          my_section->add_entry("tol_curvature", stringify_fp_sci(_tol_curvature));
          my_section->add_entry("tol_step", stringify_fp_sci(_tol_step));

          return my_section;

        }

        /**
         * \brief Sets the tolerance for the sufficient decrease in curvature
         */
        void set_tol_curvature(DataType tol_curvature)
        {
          XASSERT(tol_curvature > DataType(0));

          _tol_curvature = tol_curvature;
        }

        /**
         * \brief Sets the tolerance for the sufficient decrease in functional value
         */
        void set_tol_decrease(DataType tol_decrease)
        {
          XASSERT(tol_decrease > DataType(0));

          _tol_decrease = tol_decrease;
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

        /**
         * \brief Determines if search direction scaling is to be used
         *
         * Internally, the search direction vector gets normalised ( \f$ \| p \|_2 = 1 \f$) for numerical stablility.
         * Some line searches have a heuristic for choosing the initial step length \f$ \alpha_0 \f$, which can use
         * \f$ \| p \|_2 \f$ if \c _dir_scaling is set to \c true. This is useful when using a preconditioner
         * (e.g. Newton).
         */
        void set_dir_scaling(const bool b)
        {
          _dir_scaling = b;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          // Create temporary vectors
          _vec_initial_sol = this->_functional.create_vector_r();
          _vec_tmp = this->_functional.create_vector_r();
          _vec_pn = this->_functional.create_vector_r();
          _vec_grad = this->_functional.create_vector_r();
          BaseClass::init_symbolic();
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          // Clear temporary vectors
          _vec_initial_sol.clear();
          _vec_pn.clear();
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
          _alpha_0 = DataType(1),
          _alpha_min = DataType(0);
          _delta_0 = Math::huge<DataType>();
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
         * As this base class normalises \f$ d \f$ in _startup(), it just return \f$ \alpha \f$. Derived classes
         * NOT performing this normalisation must overwrite this.
         *
         * \returns The relative update.
         */
        virtual DataType get_rel_update() const
        {
          return Math::abs(this->_alpha_min);
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

      protected:

        /**
         * \brief Performs the startup of the iteration
         *
         * \param[out] alpha
         * The initial step size value.
         *
         * \param[out] fval
         * The initial functional value.
         *
         * \param[out] delta
         * The initial direction derivative value: \f$ \delta = \left< dir, grad \right>\f$.
         *
         * \param[in] vec_sol
         * The initial guess.
         *
         * \param[in] vec_dir
         * The search direction.
         *
         * The routine sets the initial data from "iteration 0" (needed for checking the strong Wolfe conditions
         * later) and does some error checking.
         *
         * \returns Status::progress if no error occurred.
         */
        virtual Status _startup(DataType& alpha, DataType& fval, DataType& delta, const VectorType& vec_sol, const VectorType& vec_dir)
        {
          Status status(Status::progress);

          this->_num_iter = Index(0);
          this->_vec_initial_sol.copy(vec_sol);

          this->_vec_pn.copy(vec_dir);
          this->_norm_dir = this->_vec_pn.norm2();

          // First scale so that all entries are |.| < 1
          this->_vec_pn.scale(this->_vec_pn, DataType(1)/this->_vec_pn.max_element());
          // Now scale so that vec_pn.norm2() == 1. Note that this will only be approximately true
          this->_vec_pn.scale(this->_vec_pn, DataType(1)/this->_vec_pn.norm2());

          // It is critical that fval_0 was set from outside using set_initial_fval!
          this->_fval_min = this->_fval_0;
          this->_delta_0 = this->_vec_pn.dot(this->_vec_grad);

          if(this->_delta_0 > DataType(0))
          {
            throw InternalError(__func__,__FILE__,__LINE__,"Search direction is not a descent direction: "
                +stringify_fp_sci(this->_delta_0));
          }

          // Compute initial defect. We want to minimise d^T * grad(_functional)
          this->_def_init = Math::abs(this->_delta_0);

          // Norm of the initial guess
          this->_norm_sol = vec_sol.norm2();

          this->_alpha_min = DataType(0);

          // plot?
          if(this->_plot_iter())
          {
            std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_init)
            << " : " << stringify_fp_sci(this->_fval_0) << ", ||dir|| = " << stringify_fp_sci(this->_norm_dir)
            << std::endl;
          }

          if(!Math::isfinite(this->_fval_0) || !Math::isfinite(this->_delta_0) || !Math::isfinite(this->_norm_dir))
          {
            status = Status::aborted;
          }

          // Set initial step size
          alpha = _alpha_0;
          // Other intitial values
          fval = _fval_0;
          delta = _delta_0;

          return status;
        }

        /**
         * \brief Performs the line search convergence checks using the strong Wolfe conditions
         *
         * \param[in] fval
         * The current functional value
         *
         * \param[in] df
         * The current directional derivative value: \f$ df = \left< \nabla f(x_0), d \right>\f$.
         *
         * \param[in] alpha
         * The step size of the last update.
         *
         * Also sets the current defect.
         * The strong Wolfe conditions are
         * \f[
         *    f(x) < f(x_0) + \alpha c_1 \left< \nabla f(x_0), d \right> \wedge
         *    |\left< \nabla f(x), d \right>| < - c_2 \left< \nabla f(x_0), d \right>,
         * \f]
         * where \f$c_1\f$ is known as \c tol_decrease and \f$ c_2 \f$ is known as \c tol_curvature.
         *
         * \returns Status::success if the strong Wolfe conditions hold.
         */
        virtual Status _check_convergence(const DataType fval, const DataType df, const DataType alpha)
        {
          Status status(Status::progress);

          this->_def_cur = Math::abs(df);

          Statistics::add_solver_expression(
            std::make_shared<ExpressionDefect>(this->name(), this->_def_cur, this->get_num_iter()));

          // plot?
          if(this->_plot_iter())
          {
            std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << " : " << stringify_fp_sci(fval) << " : " << stringify_fp_sci(alpha)
            << std::endl;
          }

          // ensure that the defect is neither NaN nor infinity
          if(!Math::isfinite(this->_def_cur))
          {
            status = Status::aborted;
          }

          // is diverged?
          if(this->is_diverged())
          {
            status = Status::diverged;
          }

          // If the maximum number of iterations was performed, return the iterate for the best step so far
          if(this->_num_iter >= this->_max_iter)
          {
            status = Status::max_iter;
          }

          // If the strong Wolfe conditions hold, we are successful
          if((fval < this->_fval_0 +_tol_decrease*alpha*_delta_0) && (Math::abs(df) < -_tol_curvature*_delta_0))
          {
            status = Status::success;
          }

          return status;
        }

    }; // class Linesearch
  } // namespace Solver
} // namespace FEAT
#endif // FEAT_KERNEL_SOLVER_LINESEARCH
