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
        /// descend direction vector, normalised for better numerical stability
        VectorType _vec_pn;
        /// temporary vector: Preconditioned defect
        VectorType _vec_z;
        /// temporary vector: y[k+1] = r[k+1] - r[k]
        VectorType _vec_y;

        /// Tolerance for function improvement
        DataType _tol_fval;
        /// Tolerance for the length of the update step
        DataType _tol_step;

        /// Current functional value
        DataType _fval;
        /// Functional value from the previous iteration
        DataType _fval_prev;

        /// Maximum number of restarts triggered by the strong Wolfe conditions not holding
        Index _max_num_restarts;
        /// Number of these restarts performed consecutively
        Index _num_restarts;
        /// Restart frequency, defaults to problemsize+3
        Index _restart_freq;

      public:
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
          _tol_step(Math::eps<DataType>()),
          _max_num_restarts(10),
          _num_restarts(0),
          _restart_freq(_op.columns() + Index(3)),
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
          _vec_pn = this->_op.create_vector_r();
          _vec_y = this->_op.create_vector_r();
          _vec_z = this->_op.create_vector_r();

          _vec_z.format();
          _linesearch->init_symbolic();
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          if(iterates != nullptr)
            iterates->clear();

          this->_vec_p.clear();
          this->_vec_pn.clear();
          this->_vec_r.clear();
          this->_vec_y.clear();
          this->_vec_z.clear();
          _restart_freq = Index(0);
          _linesearch->done_symbolic();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "NLCG";
        }

        /// \copydoc BaseClass::get_formated_solver_tree()
        virtual String get_formated_solver_tree() const override
        {
          String result(name());
          result += " ( "+stringify(_direction_update)+", "+_linesearch->get_formated_solver_tree();
          if(_precond != nullptr)
            result += " "+_precond->get_formated_solver_tree();
          result += " )";
          return result;
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
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // save defect
          this->_vec_r.copy(vec_def);
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
          //std::cout << std::scientific;
          //std::cout << std::setprecision(16);
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
          /// p[k+1] <- r[k+1] + _beta * p[k+1]
          DataType beta;
          /// eta[k] = <r[k], p[k]>
          DataType eta;
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

          // First scale so that all entries are |.| < 1
          this->_vec_pn.scale(this->_vec_p, DataType(1)/this->_vec_p.max_element());
          // Now scale so that normed_dir.norm2() == 1. Note that this will only be approximately true, so we
          // compute and use its norm later on anyway
          this->_vec_pn.scale(this->_vec_pn, DataType(1)/this->_vec_pn.norm2());

          ////std::cout << "dir " << *_vec_p << std::endl;
          //DataType max_elem(this->_vec_p.max_element());
          //DataType scale1 = DataType(1)/max_elem;
          //std::cout << "max_elem   " << max_elem << " " << scale1 << std::endl;
          //this->_vec_pn.scale(this->_vec_p, scale1);
          //DataType new_norm(_vec_pn.norm2());
          //DataType scale2(DataType(1)/new_norm);
          //std::cout << " new norm  " << new_norm << " " << scale2 << std::endl;
          //this->_vec_pn.scale(_vec_pn, scale2);
          ////std::cout << "normalised dir " << *_vec_pn << std::endl;

          // Compute initial eta = < p, r>
          eta = this->_vec_p.dot(this->_vec_r);

          // If the preconditioner was not pd in the first step, reset the search direction to steepest descent
          if(eta <= DataType(0))
          {
            this->_vec_p.clone(this->_vec_r);
            this->_vec_pn.scale(this->_vec_p, DataType(1)/this->_vec_p.max_element());
            this->_vec_pn.scale(this->_vec_pn, DataType(1)/this->_vec_pn.norm2());
          }

          // compute initial gamma = < z[0], r[0] >
          DataType gamma(this->_vec_z.dot(this->_vec_r));
          DataType gamma_prev(0);

          if(this->_def_init <= this->_tol_rel)
            return Status::success;

          Index its_since_restart(0);
          _num_restarts = Index(0);
          Index first_restart(_restart_freq+Index(1));
          // start iterating
          while(status == Status::progress)
          {
            IterationStats stat(*this);

            _fval_prev = _fval;

            this->_vec_y.clone(this->_vec_r);

            // Copy information to the linesearch
            _linesearch->set_initial_fval(this->_fval);
            _linesearch->set_grad_from_defect(this->_vec_r);

            // Call the linesearch to update vec_sol
            Status linesearch_status = _linesearch->correct(vec_sol, this->_vec_pn);

            // Copy back information from the linesearch
            this->_fval = _linesearch->get_final_fval();
            _linesearch->get_defect_from_grad(this->_vec_r);
            this->_vec_y.axpy(this->_vec_y, this->_vec_r, -DataType(1));

            // Log iterates if necessary
            if(iterates != nullptr)
              iterates->push_back(vec_sol.clone());

            // Compute defect norm. This also performs the convergence/divergence checks.
            status = this->_set_new_defect(this->_vec_r, vec_sol);

            if(status != Status::progress)
              return status;

            // Update preconditioner if necessary
            if(this->_precond != nullptr)
              this->_precond->prepare(vec_sol, this->_filter);

            // apply preconditioner
            if(!this->_apply_precond(_vec_z, _vec_r, this->_filter))
              return Status::aborted;

            // Save old gamma and compute new
            gamma_prev = gamma;
            gamma = this->_vec_r.dot(this->_vec_z);

            //String info("");
            bool restart(false);
            // We need to perform a steepest descent step if the Wolfe conditions do not hold (linesearch status
            // != success).
            if(linesearch_status != Status::success)
            {
              //info += " Wolfe Conditions do not hold";
              restart = true;
              _num_restarts++;
            }
            else
              _num_restarts = Index(0);

            if(_restart_freq > Index(0) && this->_num_iter >= first_restart && its_since_restart%_restart_freq == 0)
                {
                  //info += " Scheduled restart";
                  restart = true;
                  its_since_restart = Index(0);
                }
            // This needs to be done after all the checks
            ++its_since_restart;

            // Set beta = 0 or compute depending on the restart flag
            restart ? beta = DataType(0) : beta = compute_beta(gamma, gamma_prev);

            //std::cout << "Beta " << beta << info << std::endl;

            // We need to check beta again here as some direction updates might set it to zero
            // Discard the old search direction and perform (preconditioned) steepest descent
            if(beta == DataType(0))
              this->_vec_p.clone(this->_vec_z);
            else
              this->_vec_p.axpy(this->_vec_p, this->_vec_z, beta);

            // First scale so that all entries are |.| < 1
            this->_vec_pn.scale(this->_vec_p, DataType(1)/this->_vec_p.max_element());
            // Now scale so that normed_dir.norm2() == 1. Note that this will only be approximately true, so we
            // compute and use its norm later on anyway
            this->_vec_pn.scale(this->_vec_pn, DataType(1)/this->_vec_pn.norm2());

            //max_elem = this->_vec_p.max_element();
            //scale1 = DataType(1)/max_elem;
            //std::cout << "max_elem   " << max_elem << " " << scale1 << std::endl;
            //this->_vec_pn.scale(this->_vec_p, scale1);
            //new_norm = this->_vec_pn.norm2();
            //scale2 = DataType(1)/new_norm;
            //std::cout << " new norm  " << new_norm << " " << scale2 << std::endl;
            //this->_vec_pn.scale(_vec_pn, scale2);

            eta = this->_vec_p.dot(this->_vec_r);

            // Safeguard as this should not happen
            // TODO: Correct the output of the preconditioner if it turns out to not have been positive definite
            if(eta <= DataType(0))
            {
              this->_vec_p.clone(this->_vec_r);
              this->_vec_pn.scale(this->_vec_p, DataType(1)/this->_vec_p.max_element());
              this->_vec_pn.scale(this->_vec_pn, DataType(1)/this->_vec_pn.norm2());
              //return Status::aborted;
            }

          }

          // We should never come to this point
          return Status::undefined;
        }

        /**
         * \brief Computes the parameter beta for the search direction update
         *
         * \param[in] gamma
         * <r[k+1], z[k+1]>
         *
         * \param[in] gamma_prev
         * <r[k], z[k]>
         *
         * If the returned \f$ \beta = 0,\f$ the next search direction is the steepest descent direction.
         *
         * \returns The new parameter beta.
         */
        DataType compute_beta(const DataType& gamma, const DataType& gamma_prev) const
        {
          switch(_direction_update)
          {
            case NLCGDirectionUpdate::DaiYuan:
              return dai_yuan(gamma);
            case NLCGDirectionUpdate::DYHSHybrid:
              return dy_hs_hybrid(gamma);
            case NLCGDirectionUpdate::FletcherReeves:
              return fletcher_reeves(gamma, gamma_prev);
            case NLCGDirectionUpdate::HagerZhang:
              // Warning: This update is not very efficient in its current implementation.
              return hager_zhang(gamma);
            case NLCGDirectionUpdate::HestenesStiefel:
              return hestenes_stiefel();
            case NLCGDirectionUpdate::PolakRibiere:
              return polak_ribiere(gamma_prev);
            default:
              throw InternalError(__func__,__FILE__,__LINE__,
              "Unhandled direction update: "+stringify(_direction_update));
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

            Index ls_iter_digits(Math::ilog10(_linesearch->get_max_iter()));

            std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            <<  " (" << stringify(this->_linesearch->get_num_iter()).pad_front(ls_iter_digits) << ")"
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << " : " << stringify_fp_sci(this->_fval)
            << " : " << stringify_fp_sci(this->_linesearch->get_rel_update())
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

          // maximum number of iterations performed?
          if(this->_num_iter >= this->_max_iter)
            return Status::max_iter;

          // Check for convergence of the gradient norm
          if(this->is_converged())
            return Status::success;

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

          if( (_linesearch->get_rel_update() <= this->_tol_step))
          {
            //std::cout << "update step stagnation: " << _linesearch->get_rel_update() << " <= " << this->_tol_step << std::endl;
            //vec_sol.clone(_linesearch->get_initial_sol());
            return Status::success;
          }

          // If there were too many subsequent restarts, the solver is stagnated
          if(_max_num_restarts > Index(0) && _num_restarts > _max_num_restarts)
            return Status::stagnated;

          // If there were too many stagnated iterations, the solver is stagnated
          if(this->_min_stag_iter > 0 && this->_num_stag_iter > this->_min_stag_iter)
            return Status::stagnated;

          // continue iterating
          return Status::progress;
        }

        /**
         * \brief Dai-Yuan update
         *
         * \param[in] gamma
         * <r[k+1], z[k+1]>
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
         *
         * \returns The new \$f \beta_{DY}\f$
         */
        DataType dai_yuan(const DataType gamma) const
        {
          // <r[k+1] - r[k], p[k]>
          DataType eta_dif = this->_vec_p.dot(this->_vec_y);

          return gamma/(-eta_dif);
        }

        /**
         * \brief Dai-Yuan-Hestenes-Stiefel update
         *
         * \param[in] gamma
         * <r[k+1], z[k+1]>
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
         *
         * \returns The new \$f \beta_{DYHS}\f$
         */
        DataType dy_hs_hybrid(const DataType gamma) const
        {

          // <r[k+1] - r[k], z[k+1]>
          DataType gamma_dif = this->_vec_z.dot(this->_vec_y);
          // <r[k+1] - r[k], p[k]>
          DataType eta_dif = this->_vec_p.dot(this->_vec_y);

          DataType beta(Math::min(gamma/(-eta_dif), gamma_dif/(-eta_dif)));
          beta = Math::max(DataType(0), beta);

          return beta;

        }

        /**
         * \brief Fletcher-Reeves update
         *
         * \param[in] gamma
         * <r[k+1], z[k+1]>
         *
         * \param[in] gamma_prev
         * <r[k], z[k]>
         *
         * The Fletcher-Reeves update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad z_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, z_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} = \left<r_{k+1}, p_k \right>
         * \f]
         * \f[
         *   \beta_{k+1} := \frac{\gamma_{k+1}}{\gamma_k}
         * \f]
         *
         * \returns The new \$f \beta_{FR}\f$
         */
        DataType fletcher_reeves(const DataType gamma, const DataType gamma_prev) const
        {
          return gamma/gamma_prev;
        }

        /**
         * \brief Hager-Zhang update
         *
         * \param[in] gamma
         * <r[k+1], z[k+1]>
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
         *
         * \returns The new \$f \beta_{HZ}\f$
         */
        DataType hager_zhang(const DataType gamma) const
        {

          DataType thresh = DataType(0.01);

          // <r[k+1] - r[k], p[k]>
          DataType eta_dif = this->_vec_p.dot(this->_vec_y);
          // <r[k+1] - r[k], z[k+1]>
          DataType gamma_dif = this->_vec_y.dot(this->_vec_z);
          // <p[k], z[k+1]>
          DataType omega_mid = this->_vec_p.dot(this->_vec_z);
          // || r[k+1] - r[k]
          DataType norm_y = this->_vec_y.norm2();

          DataType beta = -( gamma_dif/(eta_dif) + DataType(2)*(norm_y)*omega_mid/Math::sqr(eta_dif));
          beta = Math::max(beta,-DataType(1)/(this->_vec_p.norm2()*Math::min(thresh, Math::sqrt(gamma))));

          return beta;
        }

        /**
         * \brief Hestenes-Stiefel update
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
         *
         * \returns The new \$f \beta_{HS}\f$
         */
        DataType hestenes_stiefel() const
        {
          // <r[k+1] - r[k], p[k]>
          DataType eta_dif = this->_vec_p.dot(this->_vec_y);
          // <r[k+1] - r[k], z[k+1]>
          DataType gamma_dif = this->_vec_z.dot(this->_vec_y);

          return gamma_dif/( -eta_dif);
        }

        /**
         * \brief Modified Polak-Ribiere-Polyak
         *
         * \param[in] gamma_prev
         * <r[k], z[k]>
         *
         * The Polak-Ribière-Polyak update sets
         * \f{align*}{
         *   r_k & := -\nabla f(x_k), \quad z_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, z_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} := \left<r_{k+1}, z_k \right> \\
         *   \beta_{k+1} & := \max\left\{ 0,  \frac{\gamma_{k+1} - \gamma_{k+\frac{1}{2}}}{\gamma_k} \right\}
         * \f}
         *
         * \returns The new \$f \beta_{PRP+}\f$
         */
        DataType polak_ribiere(const DataType gamma_prev) const
        {
          // <r[k+1] - r[k], z[k+1]>
          DataType gamma_dif = this->_vec_z.dot(this->_vec_y);

          return  Math::max(gamma_dif/gamma_prev, DataType(0));
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
