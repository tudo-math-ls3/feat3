#pragma once
#ifndef FEAT_KERNEL_SOLVER_QPENALTY
#define FEAT_KERNEL_SOLVER_QPENALTY 1
namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Quadratic penalty iteration
     *
     * \tparam Functional_
     * The functional (or functional) for the constrained optimisation problem.
     *
     * This class implements an outer solver for a quadratic penalty iteration. Assume we have an optimisation
     * problem of the form
     * \f[
     *   x^* = \mathrm{argmin}_{x \in D} f(x) ~ \mathrm{subject~to~} \forall i \in E: c_i(x) = 0
     * \f]
     * \f$ c_i \f$ are <em>equality constraints</em> turning the minimisation of \f$ f \f$ into a
     * <em>constrained optimisation problem</em>. To solve this constrained problem, one can solve a series of
     * unconstrained problems with the quadratic penalty iteration by defining the quadratic penalty function
     * \f[
     *   Q(x,\mu) := f(x) + \frac{\mu}{2} \sum_{i \in E} c_i^2(x)
     * \f]
     * and solving it for a series
     * \f$ (\mu_k)_k \in \mathbb{N}, \mu_k \stackrel{k \to \infty}{\longrightarrow} \infty \f$.
     * Because the penalty terms are smooth, tools from unconstrained optimisation can be applied (see
     * \cite NW06 Chapter 17.1). The greatest advantage of this is that this needs no constraint qualifications and
     * is still applicable if \f$ \mathrm{dim}( \mathrm{ker} (\nabla c_{i_1}, \dots, \nabla c_{i_M}) ) > 0 \f$,
     * which is required by all more sophisticated methods like augmented Lagragian etc.
     *
     * All of this means that the Functional_ \f$ Q \f$ must be assembled anew in every iteration and the inner solver
     * needs to take care of this.
     *
     * One critical part is the choice of \f$\mu_k\f$. It needs to start small (i.e. \f$\mu_1 = 1\f$) because the initial state
     * might violate the constraints quite strongly. The rate of increase has to be fast enough that the number of
     * penalty iterations remains small (i.e. \f$ \leq 10 \f$), but the systematic ill-conditioning of the penalty
     * function means that we cannot be too quick about this. If the increase is too quick, the previous iterate will
     * be too far away of a minimiser of the new problem and the inner solver will fail to converge. So take care.
     *
     * \author Jordi Paul
     */
    template<typename Functional_>
    class QPenalty : public IterativeSolver<typename Functional_::VectorTypeR>
    {
      public:
        /// The functional type
        typedef Functional_ FunctionalType;
        /// The input vector type
        typedef typename Functional_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename VectorType::DataType DataType;
        /// Our base class
        typedef IterativeSolver<VectorType> BaseClass;

      private:
        /// The functional that takes the penalty parameters and is assembled in the inner solver
        FunctionalType& _functional;
        /// The inner solver for the penalised, unconstrained optimisation problem
        std::shared_ptr<IterativeSolver<VectorType>> _inner_solver;
        /// We start with this initial penalty parameter
        DataType _initial_penalty_param;
        /// Maximum value the penalty parameter is allowed to take
        DataType _tol_penalty;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] functional
         * The functional for the inner constrained optimisation problem.
         *
         * \param[in] inner_solver
         * The solver for the inner penalised unconstrained optimisation problem.
         *
         * \param[in] initial_penalty_param
         * The starting penalty parameter, defaults to 1.
         */
        explicit QPenalty(FunctionalType& functional, std::shared_ptr<IterativeSolver<VectorType>> inner_solver,
        DataType initial_penalty_param = DataType(1)) :
          BaseClass("QPenalty"),
          _functional(functional),
          _inner_solver(inner_solver),
          _initial_penalty_param(initial_penalty_param),
          _tol_penalty(Math::pow(Math::huge<DataType>(), DataType(0.25)))
        {
        }

        /**
         * \brief Constructor using a PropertyMap
         *
         * \tparam Functional_
         * The type of the quadratic penalty function \f$ Q \f$.
         *
         * \param[in] section_name
         * The name of the config section, which it does not know by itself.
         *
         * \param[in] section
         * A pointer to the PropertyMap section configuring this solver.
         *
         * \param[in] functional
         * The quadratic penalty function \f$ Q \f$.
         *
         * \param[in] inner_solver
         * The inner solver for solving the penalised unconstrained optimisation problem.
         */
        explicit QPenalty(const String& section_name, PropertyMap* section,
        FunctionalType& functional, std::shared_ptr<IterativeSolver<VectorType>> inner_solver) :
          BaseClass("QPenalty", section_name, section),
          _functional(functional),
          _inner_solver(inner_solver),
          _initial_penalty_param(1),
          _tol_penalty(Math::pow(Math::huge<DataType>(), DataType(0.25)))
        {
          auto initial_penalty_param_p = section->query("initial_penalty_param");
          if(initial_penalty_param_p.second)
          {
            set_initial_penalty_param(DataType(std::stod(initial_penalty_param_p.first)));
          }

          auto tol_penalty_p = section->query("tol_penalty");
          if(tol_penalty_p.second)
          {
            set_tol_penalty(DataType(std::stod(tol_penalty_p.first)));
          }
        }

        /// Explicitly delete the copy constructor
        QPenalty(const QPenalty&) = delete;

        /**
         * \brief Virtual destructor
         */
        virtual ~QPenalty()
        {
          _inner_solver = nullptr;
        }

        /// \copydoc SolverBase::name()
        virtual String name() const override
        {
          return "QPenalty";
        }

        /// \copydoc SolverBase::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();
          _inner_solver->init_symbolic();
        }

        /// \copydoc SolverBase::init_symbolic()
        virtual void done_symbolic() override
        {
          _inner_solver->done_symbolic();
          BaseClass::done_symbolic();
        }

        /**
         * \brief Sets the initial penalty parameter
         */
        void set_initial_penalty_param(DataType initial_penalty_param)
        {
          XASSERT(initial_penalty_param > DataType(0));

          _initial_penalty_param = initial_penalty_param;
        }

        /**
         * \brief Sets the initial penalty parameter
         */
        void set_tol_penalty(DataType tol_penalty)
        {
          XASSERT(tol_penalty > DataType(0));

          _tol_penalty = tol_penalty;
        }

        /// \copydoc SolverBase::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          _functional.set_penalty_param(_initial_penalty_param);

          // clear solution vector
          vec_cor.format();

          // apply
          Status st(_apply_intern(vec_cor, vec_def));
          this->plot_summary(st);
          return st;
        }

        /// \copydoc IterativeSolver::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
        {
          _functional.set_penalty_param(_initial_penalty_param);

          // apply
          Status st(_apply_intern(vec_sol, vec_rhs));
          this->plot_summary(st);
          return st;
        }

      protected:
        /**
         * \brief Internal function, applies the solver
         *
         * \param[in, out] vec_sol
         * The initial guess, gets overwritten by the solution
         *
         * \param[in] vec_rhs
         * The right-hand side.
         *
         * \returns
         * A solver status code.
         *
         */
        virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

          const Index inner_iter_digits(Math::ilog10(_inner_solver->get_max_iter()));

          Status st(_set_initial_defect(_functional.get_constraint()));

          IterationStats stat(*this);

          // The Penalty Iteration(TM)
          while(st == Status::progress)
          {
            ++(this->_num_iter);

            DataType def_old(this->_def_cur);
            //DataType penalty_param_old(penalty_param);

            Status inner_st(_inner_solver->correct(vec_sol, vec_rhs));
            if(inner_st == Status::aborted)
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            // After the inner solver was called, the constraint is already there so we can just get() it
            this->_def_cur = _functional.get_constraint();
            DataType penalty_param(_functional.get_penalty_param());

            if(this->_plot_iter())
            {
              std::cout << this->_plot_name
              << ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
              << " (" << stringify(this->_inner_solver->get_num_iter()).pad_front(inner_iter_digits)
              << ":" << inner_st << ")"
              << " : " << stringify_fp_sci(this->_def_cur)
              << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
              << " / " << stringify_fp_sci(DataType(0.5)*Math::sqr(this->_def_cur))
              << " : " << stringify_fp_sci(penalty_param)
              << std::endl;
            }

            // ensure that the defect is neither NaN nor infinity
            if(!Math::isfinite(this->_def_cur))
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            // is diverged?
            if(this->is_diverged())
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::diverged, this->get_num_iter()));
              return Status::diverged;
            }

            // minimum number of iterations performed?
            if(this->_num_iter < this->_min_iter)
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::progress, this->get_num_iter()));
              return Status::progress;
            }

            // maximum number of iterations performed?
            if(this->_num_iter >= this->_max_iter)
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::max_iter, this->get_num_iter()));
              return Status::max_iter;
            }

            // Check for convergence
            if(this->is_converged())
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            }

            // Increment for the penalty factor by at least factor 5, very arbitrary
            penalty_param = DataType(5)*Math::sqr(penalty_param)
              *Math::max(Math::sqr( def_old/this->_def_cur),DataType(1));

            if(penalty_param >= _tol_penalty)
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::stagnated, this->get_num_iter()));
              return Status::stagnated;
            }

            _functional.set_penalty_param(penalty_param);
          }

          // We should never come to this point
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return Status::undefined;
        }

      /**
       * \brief Internal function: sets the initial defect
       *
       * \param[in] def_init_
       * The initial constraint violation.
       *
       * \returns
       * A Status code.
       */
      Status _set_initial_defect(const DataType def_init_)
      {
        // store new defect
        this->_def_init = this->_def_cur = def_init_;
        // insert special toe to signal new start of solver
        //Statistics::add_solver_toe(this->_branch, double(-1));
        //Statistics::add_solver_mpi_execute(this->_branch, double(-1));
        //Statistics::add_solver_mpi_wait(this->_branch, double(-1));
        //insert -1 as first defect, to signalize a new starting solver iteration run
        //Statistics::add_solver_defect(this->_branch, double(-1));
        //Statistics::add_solver_defect(this->_branch, double(this->_def_init));
        this->_num_iter = Index(0);
        this->_num_stag_iter = Index(0);
        Statistics::add_solver_expression(std::make_shared<ExpressionDefect>(this->name(), this->_def_init, this->get_num_iter()));

        // Plot?
        if(this->_plot_iter())
        {
          std::cout << this->_plot_name
            <<  ": " << stringify(0).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_init) << std::endl;
        }

        // Ensure that the defect is neither NaN nor infinity
        if(!Math::isfinite(this->_def_init))
        {
          return Status::aborted;
        }

        // Check if the initial defect lower than the absolute tolerance
        if(this->_def_init <= this->_tol_abs)
        {
          return Status::success;
        }

        // continue iterating
        return Status::progress;
      }

      /// Makes the BaseClass' routine available, which cannot be used but shuts up the overloaded virtual warning.
      using BaseClass::_set_initial_defect;

    };

    /**
     * \brief Creates a new QPenalty object.
     *
     * \tparam Functional_
     * The type of the quadratic penalty function \f$ Q \f$.
     *
     * \param[in] functional
     * The quadratic penalty function \f$ Q \f$.
     *
     * \param[in] inner_solver
     * The inner solver for solving the penalised unconstrained optimisation problem.
     *
     * \returns An std::shared_ptr to the new QPenalty solver object
     */
    template<typename Functional_>
    inline std::shared_ptr<QPenalty<Functional_>> new_qpenalty( Functional_& functional,
    std::shared_ptr<IterativeSolver<typename Functional_::VectorTypeR>> inner_solver,
    typename Functional_::VectorTypeR::DataType initial_penalty_param =
    typename Functional_::VectorTypeR::DataType(1))
    {
      return std::make_shared<QPenalty<Functional_>>(functional, inner_solver, initial_penalty_param);
    }

    /**
     * \brief Creates a new QPenalty object using a PropertyMap
     *
     * \tparam Functional_
     * The type of the quadratic penalty function \f$ Q \f$.
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] functional
     * The quadratic penalty function \f$ Q \f$.
     *
     * \param[in] inner_solver
     * The inner solver for solving the penalised unconstrained optimisation problem.
     *
     * \returns An std::shared_ptr to the new QPenalty solver object
     */
    template<typename Functional_>
    inline std::shared_ptr<QPenalty<Functional_>> new_qpenalty(
      const String& section_name, PropertyMap* section,
      Functional_& functional, std::shared_ptr<IterativeSolver<typename Functional_::VectorTypeR>> inner_solver)
    {
      return std::make_shared<QPenalty<Functional_>>(section_name, section, functional, inner_solver);
    }

  }
} // namespace FEAT

#endif // FEAT_KERNEL_SOLVER_QPENALTY
