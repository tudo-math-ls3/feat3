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
        static constexpr Index max_num_subs_restarts = Index(3);

      protected:
        /// Our nonlinear operator
        Operator_& _op;
        /// The filter we apply to the gradient
        const Filter_& _filter;
        /// The linesearch used along the descent direction
        LinesearchType& _linesearch;

        /// defect vector
        VectorType _vec_def;
        /// descend direction vector
        VectorType _vec_dir;
        /// temporary vector
        VectorType _vec_tmp;

        /// The mininum update we require the linesearch to make
        DataType _min_update;

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
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         * \param[in, out] precond
         * Preconditioner, defaults to nullptr. Cannot be const as internal data changes
         *
         */
        explicit NLCG(Operator_& op_, const Filter_& filter_, LinesearchType& linesearch_,
        bool keep_iterates = false, std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass("NLCG", precond),
          _op(op_),
          _filter(filter_),
          _linesearch(linesearch_),
          restart_freq(0),
          iterates(nullptr)
          {
            this->_min_stag_iter = Index(3);

            if(keep_iterates)
              iterates = new std::deque<VectorType>;

            _min_update = Math::sqrt(Math::eps<DataType>());
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
          restart_freq = _vec_def.size() + Index(1);
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
          return "NLCG-"+_linesearch.name();
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
          this->_op.compute_gradient(this->_vec_def);
          this->_vec_def.scale(this->_vec_def,DataType(-1));
          this->_filter.filter_def(this->_vec_def);

          // apply
          Status st =_apply_intern(vec_sol);

          return st;
        }

      protected:
        /**
         * \brief Internal function, applies the solver
         *
         * \param[in, out] vec_sol
         * The initial guess, gets overwritten by the solution
         *
         * This does not have a right hand side because that is contained in the gradient of the operator and we
         * always seek grad operator(vec_sol) = 0
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

          // start iterating
          Index its_since_restart(0);
          Index num_subsequent_restarts(0);
          while(status == Status::progress)
          {
            TimeStamp at;
            ++its_since_restart;

            _linesearch.correct(vec_sol, this->_vec_dir);
            if(_linesearch.get_rel_update() < _min_update)
            {
              this->_num_stag_iter++;
            }

            if(iterates != nullptr)
            {
              auto tmp = vec_sol.clone();
              iterates->push_back(std::move(tmp));
            }

            this->_op.prepare(vec_sol);
            // update defect vector
            this->_op.compute_gradient(this->_vec_def);
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
            bool restart(false);

            status = polak_ribiere(beta, gamma, restart, at);

            if(restart_freq > 0 && its_since_restart%restart_freq == 0)
              restart = true;

            if(restart)
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


            if(num_subsequent_restarts > max_num_subs_restarts)
              return Status::stagnated;

            if(this->_min_stag_iter > 0 && this->_num_stag_iter > this->_min_stag_iter)
              return Status::stagnated;

          }

          return Status::undefined;
        }

        Status polak_ribiere(DataType& beta, DataType& gamma, bool& restart, TimeStamp& at)
        {
            DataType gamma_old(gamma);
            DataType gamma_mid = this->_vec_tmp.dot(this->_vec_def);

            // apply preconditioner
            if(!this->_apply_precond(_vec_tmp, _vec_def, this->_filter))
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              return Status::aborted;
            }

            gamma = this->_vec_def.dot(this->_vec_tmp);

            beta = (gamma - gamma_mid)/gamma_old;

            restart = (beta < DataType(0));
          return Status::progress;
        }

        Status fletcher_reeves(DataType& beta, DataType& gamma, bool& restart, TimeStamp& at)
        {
          DataType gamma_old = gamma;
          DataType gamma_mid = this->_vec_dir.dot(this->vec_def);

          // apply preconditioner
          if(!this->_apply_precond(_vec_tmp, _vec_def, this->_filter))
          {
            TimeStamp bt;
            Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
            return Status::aborted;
          }

          gamma = this->_vec_def.dot(this->_vec_tmp);

          beta = gamma/gamma_old;

          restart = (gamma_mid < DataType(0));

          return Status::progress;
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
       * A Status code.
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
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new PCG object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAST_COMPILER_GNU) && (FEAST_COMPILER_GNU < 40900)
    template<typename Operator_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLCG<Operator_, Filter_, Linesearch_>> new_nlcg(
      Operator_& op, const Filter_& filter, Linesearch_& linesearch, bool keep_iterates = false)
      {
        return std::make_shared<NLCG<Operator_, Filter_, Linesearch_>>(op, filter, linesearch, keep_iterates, nullptr);
      }
    template<typename Operator_, typename Filter_, typename Linesearch_, typename Precond_>
    inline std::shared_ptr<NLCG<Operator_, Filter_, Linesearch_>> new_nlcg(
      Operator_& op, const Filter_& filter, Linesearch_& linesearch, bool keep_iterates,
      std::shared_ptr<Precond_> precond)
      {
        return std::make_shared<NLCG<Operator_, Filter_, Linesearch_>>(op, filter, linesearch, precond);
      }
#else
    template<typename Operator_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLCG<Operator_, Filter_, Linesearch_>> new_nlcg(
      Operator_& op, const Filter_& filter, Linesearch_& linesearch, bool keep_iterates = false,
      std::shared_ptr<SolverBase<typename Operator_::VectorTypeL>> precond = nullptr)
      {
        return std::make_shared<NLCG<Operator_, Filter_, Linesearch_>>(op, filter, linesearch, keep_iterates, precond);
      }
#endif
  } //namespace Solver
} // namespace FEAST

#endif
