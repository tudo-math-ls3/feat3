#pragma once
#ifndef FEAST_SOLVER_NLCG
#define FEAST_SOLVER_NLCG 1
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/solver/linesearch.hpp>

#include <deque>
namespace FEAST
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
      HestenesStiefel,
      PolakRibiere
    };

    /// \cond internal
    /**
     * \brief Streamin operator for NLCGDirectionUpdates
     */
    inline std::ostream& operator<<(std::ostream& os, NLCGDirectionUpdate update)
    {
      switch(update)
      {
        case NLCGDirectionUpdate::undefined:
          return os << "undefined";
        case NLCGDirectionUpdate::DaiYuan:
          return os << "Dai-Yuan";
        case NLCGDirectionUpdate::DYHSHybrid:
          return os << "DY-HS-Hybrid";
        case NLCGDirectionUpdate::FletcherReeves:
          return os << "Fletcher-Reeves";
        case NLCGDirectionUpdate::HestenesStiefel:
          return os << "Hestenes-Stiefel";
        case NLCGDirectionUpdate::PolakRibiere:
          return os << "Polak-Ribiere";
        default:
          return os << "-unknown-";
      }
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
     * Hestenes-Stiefel \cite HS52 and Polak-Ribière \cite PR64.
     */
    template<typename Operator_, typename Filter_, typename Linesearch_>
    class NLCG : public PreconditionedIterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        /// The nonlinear operator type
        typedef Operator_ OperatorType;
        /// The filter type
        typedef Filter_ FilterType;
        /// Our type of linesearch
        typedef Linesearch_ LinesearchType;

        /// Type of the operator's gradient has
        typedef typename Operator_::GradientType GradientType;
        /// Input type for the gradient
        typedef typename Operator_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Operator_::DataType DataType;

        /// Our baseclass
        typedef PreconditionedIterativeSolver<VectorType> BaseClass;
        /// Generic preconditioner
        typedef SolverBase<VectorType> PrecondType;
        /// Maximum number of subsequent restarts (meaning steepest descent steps) before aborting
        static constexpr Index max_num_subs_restarts = Index(10);

      protected:
        /// Our nonlinear operator
        Operator_& _op;
        /// The filter we apply to the gradient
        const Filter_& _filter;
        /// The linesearch used along the descent direction
        LinesearchType& _linesearch;

        /// Method to update the search direction, defaults to PolakRibiere
        NLCGDirectionUpdate _direction_update;

        /// defect vector
        VectorType _vec_def;
        /// descend direction vector
        VectorType _vec_dir;
        /// temporary vector
        VectorType _vec_tmp;

        /// Tolerance for function improvement
        DataType _tol_update_func;
        /// Tolerance for the lenght of the update step
        DataType _tol_update_step;

        /// <vec_def, vec_dir>
        DataType _eta;

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
        explicit NLCG(Operator_& op_, const Filter_& filter_, LinesearchType& linesearch_, NLCGDirectionUpdate du_,
        bool keep_iterates = false, std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass("NLCG", precond),
          _op(op_),
          _filter(filter_),
          _linesearch(linesearch_),
          _direction_update(du_),
          _tol_update_func(DataType(0)),
          _tol_update_step(Math::sqrt(Math::eps<DataType>())),
          restart_freq(0),
          iterates(nullptr)
          {
            //this->_min_stag_iter = Index(3);

            // If we use Fletcher-Reeves, frequent restarts are needed
            if(_direction_update == NLCGDirectionUpdate::FletcherReeves)
              restart_freq = _op.columns() + Index(1);

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
          _vec_def = this->_op.create_vector_r();
          _vec_dir = this->_op.create_vector_r();
          _vec_tmp = this->_op.create_vector_r();
          _linesearch.init_symbolic();
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          this->_vec_tmp.clear();
          this->_vec_dir.clear();
          this->_vec_def.clear();
          restart_freq = Index(0);
          _linesearch.done_symbolic();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "NLCG-"+_linesearch.name()+"-"+stringify(_direction_update);
        }

        /// \copydoc BaseClass::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // save defect
          this->_vec_def.copy(vec_def);
          //this->_system_filter.filter_def(this->_vec_def);

          // clear solution vector
          vec_cor.format();

          this->_op.prepare(vec_cor);

          // apply
          return _apply_intern(vec_cor);
        }

        /// \copydoc BaseClass::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& DOXY(vec_rhs)) override
        {
          this->_op.prepare(vec_sol);
          // compute defect
          this->_op.compute_grad(this->_vec_def);
          this->_vec_def.scale(this->_vec_def,DataType(-1));
          this->_filter.filter_def(this->_vec_def);

          // apply
          Status st =_apply_intern(vec_sol);

          return st;
        }

        /**
         * \brief Sets the direction update method
         */
        void set_direction_update(NLCGDirectionUpdate update_)
        {
          _direction_update = update_;
          // If we use Fletcher-Reeves, frequent restarts are needed
          if(_direction_update == NLCGDirectionUpdate::FletcherReeves)
            restart_freq = _vec_def.size() + Index(1);
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
          if(iterates != nullptr)
          {
            auto tmp = vec_sol.clone();
            iterates->push_back(std::move(tmp));
          }

          // compute initial defect
          Status status = this->_set_initial_defect(this->_vec_def, vec_sol);
          if(status != Status::progress)
            return status;

          // apply preconditioner to defect vector
          if(!this->_apply_precond(this->_vec_tmp, this->_vec_def, this->_filter))
            return Status::aborted;

          this->_vec_dir.clone(this->_vec_tmp);

          // compute initial gamma
          DataType gamma = this->_vec_def.dot(this->_vec_dir);
          // Compute initial eta = <d, r>
          this->_eta = gamma;

          // start iterating
          Index its_since_restart(0);
          Index num_subsequent_restarts(0);
          while(status == Status::progress)
          {
            TimeStamp at;
            ++its_since_restart;

            _linesearch.correct(vec_sol, this->_vec_dir);
            if(_linesearch.get_rel_update() < this->_tol_update_step)
              this->_num_stag_iter++;

            // Log iterates if necessary
            if(iterates != nullptr)
            {
              auto tmp = vec_sol.clone();
              iterates->push_back(std::move(tmp));
            }

            this->_op.prepare(vec_sol);
            // update defect vector
            this->_op.compute_grad(this->_vec_def);
            this->_vec_def.scale(this->_vec_def,DataType(-1));
            this->_filter.filter_def(this->_vec_def);

            // compute defect norm
            status = this->_set_new_defect(this->_vec_def, vec_sol);

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return status;
            }

            DataType beta(DataType(0));

            // Compute the new beta for the search direction update. This also checks if the computed beta is valid
            // (e.g. leads to a decrease) and applies the preconditioner to the new defect vector
            status = compute_beta(beta, gamma, at);

            // If a restart is scheduled, reset beta to 0
            if(restart_freq > 0 && its_since_restart%restart_freq == 0)
              beta = DataType(0);

            /// Restarting means discarding the new search direction and setting the new search direction to the
            // (preconditioned) steepest descent direction
            if(beta == DataType(0))
            {
              its_since_restart = Index(0);
              num_subsequent_restarts++;

              this->_vec_dir.clone(this->_vec_tmp);
            }
            else
            {
              num_subsequent_restarts = Index(0);
              this->_vec_dir.axpy(this->_vec_dir, this->_vec_tmp, beta);
            }
            std::cout << "NLCG: beta = " << stringify_fp_sci(beta) << std::endl;

            // Now that we know the new search direction, we can update _eta
            _eta = this->_vec_dir.dot(this->_vec_def);
            if(_eta <= DataType(0))
              std::cout << "New direction is not a descent direction!" << std::endl;

            // If there were too many subsequent restarts, the solver is stagnated
            if(num_subsequent_restarts > max_num_subs_restarts)
              return Status::stagnated;

            // If there were too many stagnated iterations, the solver is stagnated
            if(this->_min_stag_iter > 0 && this->_num_stag_iter > this->_min_stag_iter)
              return Status::stagnated;

          }

          // We should never come to this point
          return Status::undefined;
        }

        /**
         * \brief Computes the parameter beta for the search direction update
         *
         * \param[out] beta
         * Parameter for updating the search direction: \f$ d_{k+1} = s_{k+1} + \beta_{k+1} d_k \f$
         *
         * \param[in, out] gamma
         * This is the scalar product of the defect and the preconditioned defect
         *
         * \param[in, out] at
         * Timestamp for solver timings
         *
         * If the returned \f$ \beta = 0,\f$ the next search direction is the steepest descent direction. Many of the
         * update formulas need the old (preconditioned) defect we keep in _vec_temp, so the preconditioner is always
         * called in these routines.
         *
         * \returns A solver status code.
         */
        Status compute_beta(DataType& beta, DataType& gamma, TimeStamp& at)
        {
          switch(_direction_update)
          {
            case NLCGDirectionUpdate::DaiYuan:
              return dai_yuan(beta, gamma, at);
            case NLCGDirectionUpdate::DYHSHybrid:
              return dy_hs_hybrid(beta, gamma, at);
            case NLCGDirectionUpdate::FletcherReeves:
              return fletcher_reeves(beta, gamma, at);
            case NLCGDirectionUpdate::HestenesStiefel:
              return hestenes_stiefel(beta, gamma, at);
            case NLCGDirectionUpdate::PolakRibiere:
              return polak_ribiere(beta, gamma, at);
            default:
              return Status::aborted;
          }
        }

        /**
         * \brief Internal function: sets the new defect norm
         *
         * This function computes the defect vector's norm, increments the iteration count,
         * plots an output line to std::cout and checks whether any of the stopping criterions is fulfilled.
         *
         * \param[in] vec_def
         * The new defect vector.
         *
         * \param[in] vec_sol
         * The current solution vector approximation.
         *
         * \returns
         * A solver status code.
         */
        virtual Status _set_new_defect(const VectorType& vec_def, const VectorType& vec_sol) override
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
            this->_def_cur = this->_calc_def_norm(vec_def, vec_sol);

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

          // continue iterating
          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The Dai-Yuan update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad s_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, s_k \right>, \quad
         *   \eta_{k+\frac{1}{2}} := \left<r_{k+1}, d_k \right>, \quad \eta_k := \left<s_k, d_k\right>
         * \f]
         * \f[
         *   \beta_{k+1} := \frac{\gamma_{k+1}}{ - \eta_{k+\frac{1}{2}} - \eta_k}
         * \f]
         *
         * Note the sign changes due to this being formulated in terms of the defect rather than the function's
         * gradient.
         */
        Status dai_yuan(DataType& beta, DataType& gamma, TimeStamp& at)
        {

          // _vec_tmp still contais the old (preconditioned) defect
          DataType eta_old(this->_eta);
          DataType eta_mid(this->_vec_def.dot(this->_vec_dir));

          // apply preconditioner
          if(!this->_apply_precond(_vec_tmp, _vec_def, this->_filter))
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            return Status::aborted;
          }

          gamma = this->_vec_def.dot(this->_vec_tmp);

          beta = gamma/( - eta_mid + eta_old);

          std::cout << "Dai-Yuan:  def = " << this->_vec_def << " dir = " << this->_vec_dir << std::endl;
          std::cout << "Dai-Yuan:  eta = " << stringify_fp_sci(eta_mid) << " eta_old = " << stringify_fp_sci(eta_old) << std::endl;
          std::cout << "Dai-Yuan: beta = " << stringify_fp_sci(beta) << " gamma   = " << stringify_fp_sci(gamma) << std::endl;

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The hybrid Dai-Yuan/Hestenes-Stiefel update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad s_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, s_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} := \left< r_k, s_{k+1} \right>, \quad
         *   \eta_{k+\frac{1}{2}} := \left<r_{k+1}, d_k \right>, \quad \eta_k := \left<s_k, d_k\right>
         * \f]
         * \f[
         *   \beta_{k+1} := \max \left\{0, \min \left\{ \frac{\gamma_{k+1} - \gamma_{k+\frac{1}{2}}}{ - \eta_{k+\frac{1}{2}}+ \eta_k}, \frac{\gamma_{k+1}}{ - \eta_{k+\frac{1}{2}}+ \eta_k} \right\} \right\}
         * \f]
         *
         * Note the sign changes due to this being formulated in terms of the defect rather than the function's
         * gradient.
         */
        Status dy_hs_hybrid(DataType& beta, DataType& gamma, TimeStamp& at)
        {
          DataType eta_old = this->_eta;
          // _vec_dir still contains the old search direction
          DataType eta_mid = this->_vec_def.dot(this->_vec_dir);
          // _vec_tmp still contains the old (preconditioned) defect
          DataType gamma_mid = this->_vec_tmp.dot(this->_vec_def);

          // apply preconditioner
          if(!this->_apply_precond(_vec_tmp, _vec_def, this->_filter))
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            return Status::aborted;
          }

          gamma = this->_vec_def.dot(this->_vec_tmp);

          beta = Math::max(DataType(0), Math::min(gamma, gamma - gamma_mid)/( -eta_mid + eta_old));

          //std::cout << "DY-HS-Hybrid:  def = " << this->_vec_def << " dir = " << this->_vec_dir << std::endl;
          //std::cout << "DY-HS-Hybrid:  eta = " << stringify_fp_sci(_eta) << " eta_old = " << stringify_fp_sci(eta_old) << std::endl;
          //std::cout << "DY-HS-Hybrid: beta = " << stringify_fp_sci(beta) << " gamma   = " << stringify_fp_sci(gamma) << std::endl;

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         * The Fletcher-Reeves update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad s_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, s_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} = \left<r_{k+1}, d_k \right>
         * \f]
         * \f[
         *   \beta_{k+1} :=
         *   \begin{cases}
         *     0, & \gamma_{k+\frac{1}{2}} \leq 0 \\
         *     \frac{\gamma_{k+1}}{\gamma_k}, & \mathrm{else}
         *   \end{cases}
         * \f]
         */
        Status fletcher_reeves(DataType& beta, DataType& gamma, TimeStamp& at)
        {

          // apply preconditioner
          if(!this->_apply_precond(_vec_tmp, _vec_def, this->_filter))
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            return Status::aborted;
          }

          DataType gamma_old = gamma;
          DataType gamma_mid = this->_vec_tmp.dot(this->_vec_dir);

          if(gamma_mid < -gamma_old)
            beta = DataType(0);
          else
          {
            gamma = this->_vec_def.dot(this->_vec_tmp);
            beta = gamma/gamma_old;
          }

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The Hestenes-Stiefel update sets
         * \f[
         *   r_k := -\nabla f(x_k), \quad s_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, s_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} := \left< r_k, s_{k+1} \right>, \quad
         *   \eta_{k+\frac{1}{2}} := \left<r_{k+1}, d_k \right>, \quad \eta_k := \left<s_k, d_k\right>
         * \f]
         * \f[
         *   \beta_{k+1} := \max \left\{0, \frac{\gamma_{k+1} - \gamma_{k+\frac{1}{2}}}{ - \eta_{k+\frac{1}{2}}
         *   + \eta_k} \right\}
         * \f]
         *
         * Note the sign changes due to this being formulated in terms of the defect rather than the function's
         * gradient.
         */
        Status hestenes_stiefel(DataType& beta, DataType& gamma, TimeStamp& at)
        {
          DataType eta_old = this->_eta;
          // _vec_dir still contains the old search direction
          DataType eta_mid = this->_vec_def.dot(this->_vec_dir);
          // _vec_tmp still contains the old (preconditioned) defect
          DataType gamma_mid = this->_vec_tmp.dot(this->_vec_def);

          // apply preconditioner
          if(!this->_apply_precond(_vec_tmp, _vec_def, this->_filter))
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            return Status::aborted;
          }

          gamma = this->_vec_def.dot(this->_vec_tmp);

          beta = Math::max(DataType(0), (gamma - gamma_mid)/( -eta_mid + eta_old));

          //std::cout << "Hestenes-Stiefel:  def = " << this->_vec_def << " dir = " << this->_vec_dir << std::endl;
          //std::cout << "Hestenes-Stiefel:  eta = " << stringify_fp_sci(eta_mid) << " eta_old = " << stringify_fp_sci(eta_old) << std::endl;
          //std::cout << "Hestenes-Stiefel: beta = " << stringify_fp_sci(beta) << " gamma   = " << stringify_fp_sci(gamma) << std::endl;

          return Status::progress;
        }

        /**
         * \copydoc compute_beta()
         *
         * The Polak-Ribière update sets
         * \f{align*}{
         *   r_k & := -\nabla f(x_k), \quad s_k := M^{-1} r_k, \quad  \gamma_k := \left< r_k, s_k \right>, \quad
         *   \gamma_{k+\frac{1}{2}} := \left<r_{k+1}, s_k \right> \\
         *   \beta_{k+1} & := \max\left\{ 0,  \frac{\gamma_{k+1} - \gamma_{k+\frac{1}{2}}}{\gamma_k} \right\}
         * \f}
         */
        Status polak_ribiere(DataType& beta, DataType& gamma, TimeStamp& at)
        {
          DataType gamma_old(gamma);
          // vec_tmp still contains the old (preconditioned) defect
          DataType gamma_mid = this->_vec_tmp.dot(this->_vec_def);

          // apply preconditioner
          if(!this->_apply_precond(_vec_tmp, _vec_def, this->_filter))
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            return Status::aborted;
          }

          gamma = this->_vec_def.dot(this->_vec_tmp);

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
#if defined(FEAST_COMPILER_GNU) && (FEAST_COMPILER_GNU < 40900)
    template<typename Operator_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLCG<Operator_, Filter_, Linesearch_>> new_nlcg(
      Operator_& op, const Filter_& filter, Linesearch_& linesearch,
      NLCGDirectionUpdate direction_update = NLCGDirectionUpdate::PolakRibiere,
      bool keep_iterates = false)
      {
        return std::make_shared<NLCG<Operator_, Filter_, Linesearch_>>(op, filter, linesearch, direction_update,
        keep_iterates, nullptr);
      }
    template<typename Operator_, typename Filter_, typename Linesearch_, typename Precond_>
    inline std::shared_ptr<NLCG<Operator_, Filter_, Linesearch_>> new_nlcg(
      Operator_& op, const Filter_& filter, Linesearch_& linesearch,
      NLCGDirectionUpdate direction_update,
      bool keep_iterates,
      std::shared_ptr<Precond_> precond)
      {
        return std::make_shared<NLCG<Operator_, Filter_, Linesearch_>>(op, filter, linesearch, direction_update,
        keep_iterates, precond);
      }
#else
    template<typename Operator_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLCG<Operator_, Filter_, Linesearch_>> new_nlcg(
      Operator_& op, const Filter_& filter, Linesearch_& linesearch,
      NLCGDirectionUpdate direction_update = NLCGDirectionUpdate::PolakRibiere,
      bool keep_iterates = false,
      std::shared_ptr<SolverBase<typename Operator_::VectorTypeL>> precond = nullptr)
      {
        return std::make_shared<NLCG<Operator_, Filter_, Linesearch_>>(op, filter, linesearch, direction_update,
        keep_iterates, precond);
      }
#endif
  } //namespace Solver
} // namespace FEAST

#endif
