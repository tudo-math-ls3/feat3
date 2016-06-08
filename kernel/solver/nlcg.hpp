#pragma once
#ifndef FEAT_SOLVER_NLCG
#define FEAT_SOLVER_NLCG 1
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlopt_precond.hpp>

#include <deque>
namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Enum for NLCG search direction updates
     */
    enum class NLCGDirectionUpdate
    {
      undefined = 0,
      DaiYuan,
      DYHSHybrid,
      FletcherReeves,
      HagerZhang, // warning: Experimental, inefficient
      HestenesStiefel,
      PolakRibiere
    };

    /// \cond internal
    /**
     * \brief Streaming operator for NLCGDirectionUpdates
     */
    inline std::ostream& operator<<(std::ostream& os, NLCGDirectionUpdate update)
    {
      switch(update)
      {
        case NLCGDirectionUpdate::undefined:
          return os << "undefined";
        case NLCGDirectionUpdate::DaiYuan:
          return os << "DaiYuan";
        case NLCGDirectionUpdate::DYHSHybrid:
          return os << "DYHSHybrid";
        case NLCGDirectionUpdate::FletcherReeves:
          return os << "FletcherReeves";
        case NLCGDirectionUpdate::HagerZhang:
          return os << "HagerZhang";
        case NLCGDirectionUpdate::HestenesStiefel:
          return os << "HestenesStiefel";
        case NLCGDirectionUpdate::PolakRibiere:
          return os << "PolakRibiere";
        default:
          return os << "-unknown-";
      }
    }

    inline void operator<<(NLCGDirectionUpdate& update, const String& update_name)
    {
        if(update_name == "undefined")
          update = NLCGDirectionUpdate::undefined;
        else if(update_name == "DaiYuan")
          update = NLCGDirectionUpdate::DaiYuan;
        else if(update_name == "DYHSHybrid")
          update = NLCGDirectionUpdate::DYHSHybrid;
        else if(update_name == "FletcherReeves")
          update = NLCGDirectionUpdate::FletcherReeves;
        else if(update_name == "HagerZhang")
          update = NLCGDirectionUpdate::HagerZhang;
        else if(update_name == "HestenesStiefel")
          update = NLCGDirectionUpdate::HestenesStiefel;
        else if(update_name == "PolakRibiere")
          update = NLCGDirectionUpdate::PolakRibiere;
        else
          throw InternalError(__func__, __FILE__, __LINE__, "Unknown NLCGDirectionUpdate identifier string "
              +update_name);
    }

    /// \endcond

    /**
     * \brief Nonlinear Conjugate Gradient method for finding a minimum of an operator's gradient
     *
     * \tparam Operator_
     * Nonlinear Operator to minimise the gradient of
     *
     * \tparam Filter_
     * Filter to apply to the operator's gradient
     *
     * \tparam Linesearch_
     * Type of linesearch to use along the descent direction
     *
     * See \cite NW06 for an overview of optimisation techniques.
     *
     * Possible update strategies for the search direction are Dai-Yuan \cite DY99, Fletcher-Reeves \cite FR64,
     * Hager-Zhang \cite HZ05, Hestenes-Stiefel \cite HS52 and Polak-Ribiere \cite PR64.
     */
    template<typename Operator_, typename Filter_>
    class NLCG : public PreconditionedIterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        /// The nonlinear operator type
        typedef Operator_ OperatorType;
        /// The filter type
        typedef Filter_ FilterType;
        /// Our type of linesearch
        typedef Solver::Linesearch<Operator_, Filter_> LinesearchType;

        /// Type of the operator's gradient has
        typedef typename Operator_::GradientType GradientType;
        /// Input type for the gradient
        typedef typename Operator_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Operator_::DataType DataType;

        /// Our baseclass
        typedef PreconditionedIterativeSolver<VectorType> BaseClass;
        /// Generic preconditioner
        typedef NLOptPrecond<typename Operator_::VectorTypeL, Filter_> PrecondType;
        /// Default search direction update
        static constexpr NLCGDirectionUpdate direction_update_default = NLCGDirectionUpdate::DYHSHybrid;

      protected:
        /// Our nonlinear operator
        Operator_& _op;
        /// The filter we apply to the gradient
        Filter_& _filter;
        /// The linesearch used along the descent direction
        std::shared_ptr<LinesearchType> _linesearch;
        /// This will be the preconditioner, or a nullptr. We need to save it ourselves because we cannot access the
        /// prepare() routine through the SolverBase pointer in our BaseClass
        std::shared_ptr<PrecondType> _precond;

        /// Method to update the search direction
        NLCGDirectionUpdate _direction_update;

        /// defect vector
        VectorType _vec_r;
        /// descend direction vector
        VectorType _vec_p;
        /// temporary vector
        VectorType _vec_z;

        /// Tolerance for function improvement
        DataType _tol_fval;
        /// Tolerance for the length of the update step
        DataType _tol_step;

        /// vec_p <- vec_r + _beta * vec_p
        DataType _beta;
        /// _eta = <vec_r, vec_p>
        DataType _eta;
        /// Current functional value
        DataType _fval;
        /// Functional value from the previous iteration
        DataType _fval_prev;

        /// Number of subsequent steepest descent steps
        Index _num_restarts;
        /// Maximum number of subsequent restarts (meaning steepest descent steps) before aborting
        Index _max_num_restarts;

      public:
        /// Restart frequency, defaults to problemsize+1
        Index restart_freq;
        /// For debugging purposes, all iterates can be logged to here
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
         * \param[in, out] linesearch_
         * The linesearch to be used, cannot be const as internal data changes
         *
         * \param[in] du_
         * Which direction update to use
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         * \param[in, out] precond
         * Preconditioner, defaults to nullptr. Cannot be const as internal data changes
         *
         */
        explicit NLCG(Operator_& op_, Filter_& filter_, std::shared_ptr<LinesearchType> linesearch_,
        const NLCGDirectionUpdate du_ = direction_update_default,
        bool keep_iterates = false, std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass("NLCG", precond),
          _op(op_),
          _filter(filter_),
          _linesearch(linesearch_),
          _precond(precond),
          _direction_update(du_),
          _tol_fval(DataType(0)),
          _tol_step(Math::sqrt(Math::eps<DataType>())),
          _beta(0),
          _num_restarts(0),
          _max_num_restarts(0),
          restart_freq(_op.columns() + Index(4)),
          iterates(nullptr)
          {
            XASSERT(_linesearch != nullptr);

            this->_min_stag_iter = 0;

            if(keep_iterates)
              iterates = new std::deque<VectorType>;
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~NLCG()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();
          // create three temporary vectors
          _vec_r = this->_op.create_vector_r();
          _vec_p = this->_op.create_vector_r();
          _vec_z = this->_op.create_vector_r();
          _linesearch->init_symbolic();
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          if(iterates != nullptr)
            iterates->clear();

          this->_vec_z.clear();
          this->_vec_p.clear();
          this->_vec_r.clear();
          restart_freq = Index(0);
          _linesearch->done_symbolic();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "NLCG-"+_linesearch->name()+"-"+stringify(_direction_update);
        }

        /**
         * \brief Sets the tolerance for function value improvement
         *
         * \param[in] tol_fval
         * New tolerance for function value improvement.
         *
         * The convergence check is against the maximum of the absolute and relative function value.
         *
         */
        void set_tol_fval(DataType tol_fval)
        {
          _tol_fval = tol_fval;
        }

        /**
         * \brief Sets the tolerance for the linesearch step size.
         *
         * \param[in] tol_step
         * New tolerance for the linesearch step size.
         *
         * If the linesearch fails to find a new iterate because its relative update is too small, the direction
         * update will fail to produce a new search direction so the NLCG has to be terminated.
         *
         */
        void set_tol_step(DataType tol_step)
        {
          _tol_step = tol_step;
        }

        /// \copydoc BaseClass::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_r) override
        {
          // save defect
          this->_vec_r.copy(vec_r);
          //this->_system_filter.filter_def(this->_vec_r);

          // clear solution vector
          vec_cor.format();

          // Evaluate the operator at the new state
          this->_op.prepare(vec_cor, this->_filter);

          // Prepare the preconditioner (if any)
          if(this->_precond != nullptr)
            this->_precond->prepare(vec_cor, this->_filter);

          // apply
          return _apply_intern(vec_cor);
        }

        /// \copydoc BaseClass::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& DOXY(vec_rhs)) override
        {
          // Evaluate the operator at the new state
          this->_op.prepare(vec_sol, this->_filter);

          // Compute defect
          this->_op.compute_grad(this->_vec_r);
          this->_vec_r.scale(this->_vec_r,DataType(-1));
          this->_filter.filter_def(this->_vec_r);

          // Prepare the preconditioner (if any)
          if(this->_precond != nullptr)
            this->_precond->prepare(vec_sol, this->_filter);

          // apply
          Status st =_apply_intern(vec_sol);

          return st;
        }

        /**
         * \brief Sets the direction update method
         *
         * \param[in] update_
         * NLCGDirectionUpdate to use.
         */
        void set_direction_update(NLCGDirectionUpdate update_)
        {
          _direction_update = update_;
        }

      protected:
        /**
         * \brief Internal function, applies the solver
         *
         * \param[in, out] vec_sol
         * The initial guess, gets overwritten by the solution
         *
         * \returns
         * A solver status code.
         *
         * This does not have a right hand side because that is contained in the gradient of the operator and we
         * always seek grad operator(vec_sol) = 0
         *
         */
        virtual Status _apply_intern(VectorType& vec_sol)
        {
          // Reset member variables in the LineSearch
          _linesearch->reset();

          if(iterates != nullptr)
            iterates->push_back(std::move(vec_sol.clone()));

          // Compute intitial function value
          this->_fval = this->_op.compute_func();
          // Compute initial defect
          Status status = this->_set_initial_defect(this->_vec_r, vec_sol);
          if(status != Status::progress)
            return status;

          // apply preconditioner to defect vector
          if(!this->_apply_precond(this->_vec_z, this->_vec_r, this->_filter))
            return Status::aborted;

          // The first direction has to be the (preconditioned) steepest descent direction
          this->_vec_p.clone(this->_vec_z);

          // Compute initial eta = < p, r>
          this->_eta = this->_vec_p.dot(this->_vec_r);

          // If the preconditioner was not pd in the first step, reset the search direction to steepest descent
          if(_eta <= DataType(0))
          {
            this->_vec_p.clone(this->_vec_r);
            _eta = this->_vec_p.dot(this->_vec_r);
          }

          // compute initial gamma = < z, r >
          DataType gamma = this->_vec_z.dot(this->_vec_r);
          if(this->_def_init <= this->_tol_rel)
            return Status::success;

          Index its_since_restart(0);
          _num_restarts = Index(0);
          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            ++its_since_restart;
            _fval_prev = _fval;
            _beta = DataType(0);

            // Copy information to the linesearch
            _linesearch->set_initial_fval(this->_fval);
            _linesearch->set_grad_from_defect(this->_vec_r);

            // Call the linesearch to update vec_sol
            Status linesearch_status = _linesearch->correct(vec_sol, this->_vec_p);

            // If the linesearch failed to make progress, the new iterate is too close to the old iterate to compute
            // a new search direction etc. so we have to abort without updating the solution
            if( (_linesearch->get_rel_update() < this->_tol_step) && (Math::abs(_beta) < Math::eps<DataType>())
                && (linesearch_status != Status::success ))
            {
              vec_sol.clone(_linesearch->get_initial_sol());
              return Status::stagnated;
            }

            // Copy back information from the linesearch
            this->_fval = _linesearch->get_final_fval();
            _linesearch->get_defect_from_grad(this->_vec_r);

            // Log iterates if necessary
            if(iterates != nullptr)
              iterates->push_back(vec_sol.clone());

            // Compute defect norm. This also performs the convergence/divergence checks.
            status = this->_set_new_defect(this->_vec_r, vec_sol);

            // Something might have gone wrong in applying the preconditioner, so check status again.
            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

            // Update preconditioner if necessary
            if(this->_precond != nullptr)
              this->_precond->prepare(vec_sol, this->_filter);

            status = compute_beta(_beta, gamma);

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

            // If a restart is scheduled, reset beta to 0
            if( (restart_freq > 0 && its_since_restart%restart_freq == 0) ||
                linesearch_status != Status::success)
            {
              _beta = DataType(0);
              its_since_restart++;
            }

            // Restarting means discarding the new search direction and setting the new search direction to the
            // (preconditioned) steepest descent direction
            if(_beta == DataType(0))
            {
              // Uncomment the line below to deviate from the ALGLIBMinCG behaviour
              // its_since_restart = Index(1);
              _num_restarts++;
              // Discard the old search direction and perform (preconditioned) steepest descent
              this->_vec_p.clone(this->_vec_z);
            }
            else
            {
              _num_restarts = Index(0);
              this->_vec_p.axpy(this->_vec_p, this->_vec_z, _beta);
            }


            // Now that we know the new search direction, we can update _eta
            _eta = this->_vec_p.dot(this->_vec_r);

            // Safeguard as this should not happen
            // TODO: Correct the output of the preconditioner if it turns out to not have been positive definite
            if(_eta <= DataType(0))
            {
              this->_vec_p.clone(this->_vec_r);
              _eta = this->_vec_p.dot(this->_vec_r);
              //return Status::aborted;
            }

          }

          // We should never come to this point
          return Status::undefined;
        }

        /**
         * \brief Computes the parameter beta for the search direction update
         *
         * \param[out] beta
         * Parameter for updating the search direction: \f$ d_{k+1} = z_{k+1} + \beta_{k+1} p_k \f$
         *
         * \param[in, out] gamma
         * This is the scalar product of the defect and the preconditioned defect
         *
         * If the returned \f$ \beta = 0,\f$ the next search direction is the steepest descent direction. Many of the
         * update formulas need the old (preconditioned) defect we keep in _vec_temp, so the preconditioner is always
         * called in these routines.
         *
         * \returns A solver status code.
         */
        Status compute_beta(DataType& beta, DataType& gamma)
        {
          switch(_direction_update)
          {
            case NLCGDirectionUpdate::DaiYuan:
              return dai_yuan(beta, gamma);
            case NLCGDirectionUpdate::DYHSHybrid:
              return dy_hs_hybrid(beta, gamma);
            case NLCGDirectionUpdate::FletcherReeves:
              return fletcher_reeves(beta, gamma);
            case NLCGDirectionUpdate::HagerZhang:
              // Warning: This update is not very efficient in its current implementation.
              return hager_zhang(beta, gamma);
            case NLCGDirectionUpdate::HestenesStiefel:
              return hestenes_stiefel(beta, gamma);
            case NLCGDirectionUpdate::PolakRibiere:
              return polak_ribiere(beta, gamma);
            default:
              return Status::undefined;
          }
        }

        /**
         * \brief Internal function: sets the new defect norm
         *
         * This function computes the defect vector's norm, increments the iteration count,
         * plots an output line to std::cout and checks whether any of the stopping criterions is fulfilled.
         *
         * \param[in] vec_r
         * The new defect vector.
         *
         * \param[in] vec_sol
         * The current solution vector approximation.
         *
         * \returns
         * A solver status code.
         */
        virtual Status _set_new_defect(const VectorType& vec_r, const VectorType& vec_sol) override
        {
          // increase iteration count
          ++this->_num_iter;

          // Check for convergence wrt. the function value improvement if _tol_fval says so
          if(_tol_fval > DataType(0))
          {
            // This is the factor for the relative function value
            DataType scale(Math::max(this->_fval, _fval_prev));
            // Make sure it is at least 1
            scale = Math::max(scale, DataType(1));
            // Check for success
            if(Math::abs(_fval_prev - this->_fval) <= _tol_fval*scale)
              return Status::success;
          }

          // first, let's see if we have to compute the defect at all
          bool calc_def = false;
          calc_def = calc_def || (this->_min_iter < this->_max_iter);
          calc_def = calc_def || this->_plot;
          calc_def = calc_def || (this->_min_stag_iter > Index(0));

          // compute new defect
          if(calc_def)
            this->_def_cur = this->_calc_def_norm(vec_r, vec_sol);

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

          // Check for convergence of the gradient norm
          if(this->is_converged())
            return Status::success;

          // If there were too many subsequent restarts, the solver is stagnated
          if(_max_num_restarts > Index(0) && _num_restarts > _max_num_restarts)
            return Status::stagnated;

          // If there were too many stagnated iterations, the solver is stagnated
          if(this->_min_stag_iter > 0 && this->_num_stag_iter > this->_min_stag_iter)
            return Status::stagnated;

          // maximum number of iterations performed?
          if(this->_num_iter >= this->_max_iter)
            return Status::max_iter;

          // continue iterating
          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The Dai-Yuan update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad z_k := M^{-1} r_k, \quad  \gamma_{k+1} := \left< r_{k+1}, z_{k+1} \right>,
         *   \quad \eta_{k+\frac{1}{2}} := \left<r_{k+1}, p_k \right>, \quad \eta_k := \left<r_k, d_k\right>
         * \f]
         * \f[
         *   \beta_{k+1} := \frac{\gamma_{k+1}}{ - \eta_{k+\frac{1}{2}} - \eta_k}
         * \f]
         *
         * Note the sign changes due to this being formulated in terms of the defect rather than the function's
         * gradient.
         */
        Status dai_yuan(DataType& beta, DataType& gamma)
        {

          // _vec_z still contais the old (preconditioned) defect
          DataType eta_old(this->_eta);
          // eta_mid = < r[k+1], d[k] >
          DataType eta_mid(this->_vec_r.dot(this->_vec_p));

          // apply preconditioner
          if(!this->_apply_precond(_vec_z, _vec_r, this->_filter))
            return Status::aborted;

          gamma = this->_vec_r.dot(this->_vec_z);

          beta = gamma/( - eta_mid + eta_old);

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The hybrid Dai-Yuan/Hestenes-Stiefel update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad z_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, z_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} := \left< r_k, z_{k+1} \right>, \quad
         *   \eta_{k+\frac{1}{2}} := \left<r_{k+1}, p_k \right>, \quad \eta_k := \left<r_k, d_k\right>
         * \f]
         * \f[
         *   \beta_{k+1} := \max \left\{0, \min \left\{ \frac{\gamma_{k+1} - \gamma_{k+\frac{1}{2}}}{ - \eta_{k+\frac{1}{2}}+ \eta_k}, \frac{\gamma_{k+1}}{ - \eta_{k+\frac{1}{2}}+ \eta_k} \right\} \right\}
         * \f]
         *
         * Note the sign changes due to this being formulated in terms of the defect rather than the function's
         * gradient.
         */
        Status dy_hs_hybrid(DataType& beta, DataType& gamma)
        {
          DataType eta_old = this->_eta;
          // _vec_p still contains the old search direction
          DataType eta_mid = this->_vec_r.dot(this->_vec_p);
          // _vec_z still contains the old (preconditioned) defect
          DataType gamma_mid = this->_vec_z.dot(this->_vec_r);

          // apply preconditioner
          if(!this->_apply_precond(_vec_z, _vec_r, this->_filter))
            return Status::aborted;

          gamma = this->_vec_r.dot(this->_vec_z);

          beta = Math::max(DataType(0), Math::min(gamma, gamma - gamma_mid)/( -eta_mid + eta_old));

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         * The Fletcher-Reeves update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad z_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, z_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} = \left<r_{k+1}, p_k \right>
         * \f]
         * \f[
         *   \beta_{k+1} := \frac{\gamma_{k+1}}{\gamma_k}
         * \f]
         */
        Status fletcher_reeves(DataType& beta, DataType& gamma)
        {
          DataType gamma_old(gamma);

          // apply preconditioner
          if(!this->_apply_precond(_vec_z, _vec_r, this->_filter))
            return Status::aborted;

          gamma = this->_vec_r.dot(this->_vec_z);
          beta = gamma/gamma_old;

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The Hager-Zhang update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad z_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, z_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} := \left< r_k, z_{k+1} \right>, \quad
         *   \eta_{k+\frac{1}{2}} := \left<r_{k+1}, p_k \right>, \quad \eta_k := \left<r_k, d_k\right>
         * \f]
         * \f[
         *   \bar{\beta}_{k+1} := \left< r_{k+1} -r_k - \frac{2 \|r_{k+1}\|^2}{\left< r_{k+1} - r_k, p_k \right>}
         *   p_k, 1/(<r_{k+1} - r_k, p_k>) r_k \right>,
         * \f]
         * \f[
         *   \beta_{k+1} := max\{ \beta, -1/ \left( \| p_k \| min\{ \mathrm{thresh}, \| r_{k+1} \| \} \right) \},
         * \f]
         *
         * where thresh is a predefined threshold parameter (i.e. 0.01).
         *
         * Note the sign changes due to this being formulated in terms of the defect rather than the function's
         * gradient.
         *
         * \warning This update strategy needs more precision in the linesearch than all others and is not very
         * efficient in the current implementation. This serves more as a starting point for further improvement.
         */
        Status hager_zhang(DataType& beta, DataType& gamma)
        {
          DataType thresh = DataType(0.01);
          DataType eta_old(this->_eta);
          // _vec_p still contains the old search direction
          DataType eta_mid = this->_vec_r.dot(this->_vec_p);
          // _vec_z still contains the old (preconditioned) defect
          DataType gamma_mid = this->_vec_z.dot(this->_vec_r);

          // apply preconditioner
          if(!this->_apply_precond(_vec_z, _vec_r, this->_filter))
            return Status::aborted;

          gamma = this->_vec_r.dot(this->_vec_z);

          beta = -( (gamma - gamma_mid)/(eta_mid - eta_old)
            + DataType(2)*(Math::sqr(gamma) - DataType(2)*gamma*gamma_mid + Math::sqr(gamma_mid))*eta_mid
            /Math::sqr(eta_mid-eta_old));

          beta = Math::max(beta,-DataType(1)/(this->_vec_p.norm2()*Math::min(thresh, Math::sqrt(gamma))));

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The Hestenes-Stiefel update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad z_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, z_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} := \left< r_k, z_{k+1} \right>, \quad
         *   \eta_{k+\frac{1}{2}} := \left<r_{k+1}, p_k \right>, \quad \eta_k := \left<r_k, d_k\right>
         * \f]
         * \f[
         *   \beta_{k+1} := \max \left\{0, \frac{\gamma_{k+1} - \gamma_{k+\frac{1}{2}}}{ - \eta_{k+\frac{1}{2}}
         *   + \eta_k} \right\}
         * \f]
         *
         * Note the sign changes due to this being formulated in terms of the defect rather than the function's
         * gradient.
         */
        Status hestenes_stiefel(DataType& beta, DataType& gamma)
        {
          DataType eta_old = this->_eta;
          // _vec_p still contains the old search direction
          DataType eta_mid = this->_vec_r.dot(this->_vec_p);
          // _vec_z still contains the old (preconditioned) defect
          DataType gamma_mid = this->_vec_z.dot(this->_vec_r);

          // apply preconditioner
          if(!this->_apply_precond(_vec_z, _vec_r, this->_filter))
            return Status::aborted;

          gamma = this->_vec_r.dot(this->_vec_z);

          beta = Math::max(DataType(0), (gamma - gamma_mid)/( -eta_mid + eta_old));

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The Polak-Ribière update sets
         * \f{align*}{
         *   r_k & := -\nabla f(x_k), \quad z_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, z_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} := \left<r_{k+1}, z_k \right> \\
         *   \beta_{k+1} & := \max\left\{ 0,  \frac{\gamma_{k+1} - \gamma_{k+\frac{1}{2}}}{\gamma_k} \right\}
         * \f}
         */
        Status polak_ribiere(DataType& beta, DataType& gamma)
        {
          DataType gamma_old(gamma);
          // vec_z still contains the old (preconditioned) defect
          DataType gamma_mid = this->_vec_z.dot(this->_vec_r);

          // apply preconditioner
          if(!this->_apply_precond(_vec_z, _vec_r, this->_filter))
            return Status::aborted;

          gamma = this->_vec_r.dot(this->_vec_z);

          beta = Math::max((gamma - gamma_mid)/gamma_old, DataType(0));

          return Status::progress;
        }

    }; // class NLCG

    /**
     * \brief Creates a new NLCG solver object
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \param[in] linesearch
     * The linesearch to use.
     *
     * \param[in] direction_update
     * The direction update to use, defaults to Polak-Ribiere.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new NLCG object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Operator_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLCG<Operator_, Filter_>> new_nlcg(
      Operator_& op, Filter_& filter, Linesearch_& linesearch,
      NLCGDirectionUpdate direction_update = NLCG<Operator_, Filter_>::direction_update_default,
      bool keep_iterates = false)
      {
        return std::make_shared<NLCG<Operator_, Filter_>>(op, filter, linesearch, direction_update,
        keep_iterates, nullptr);
      }
    template<typename Operator_, typename Filter_, typename Linesearch_, typename Precond_>
    inline std::shared_ptr<NLCG<Operator_, Filter_>> new_nlcg(
      Operator_& op, Filter_& filter, Linesearch_& linesearch,
      NLCGDirectionUpdate direction_update,
      bool keep_iterates,
      std::shared_ptr<Precond_> precond)
      {
        return std::make_shared<NLCG<Operator_, Filter_>>(op, filter, linesearch, direction_update,
        keep_iterates, precond);
      }
#else
    template<typename Operator_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLCG<Operator_, Filter_>> new_nlcg(
      Operator_& op, Filter_& filter, Linesearch_& linesearch,
      NLCGDirectionUpdate direction_update = NLCG<Operator_, Filter_>::direction_update_default,
      bool keep_iterates = false,
      std::shared_ptr<NLOptPrecond<typename Operator_::VectorTypeL, Filter_>> precond = nullptr)
      {
        return std::make_shared<NLCG<Operator_, Filter_>>(op, filter, linesearch, direction_update,
        keep_iterates, precond);
      }
#endif
  } //namespace Solver
} // namespace FEAT

#endif
