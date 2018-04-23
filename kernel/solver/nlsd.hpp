#pragma once
#ifndef FEAT_SOLVER_NLSD
#define FEAT_SOLVER_NLSD 1
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nloptls.hpp>
#include <kernel/solver/nlopt_precond.hpp>

#include <deque>
namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Nonlinear Steepest Descent method for finding a minimum of an functional's gradient
     *
     * \tparam Functional_
     * Nonlinear Functional to minimise the gradient of
     *
     * \tparam Filter_
     * Filter to apply to the functional's gradient
     *
     * See \cite NW06 for an overview of optimisation techniques.
     *
     */
    template<typename Functional_, typename Filter_>
    class NLSD : public NLOptLS<Functional_, Filter_>
    {
      public:
        /// The nonlinear functional type
        typedef Functional_ FunctionalType;
        /// The filter type
        typedef Filter_ FilterType;
        /// The baseclass for all applicable linesearches
        typedef Linesearch<FunctionalType, FilterType> LinesearchType;

        /// Type of the functional's gradient has
        typedef typename Functional_::GradientType GradientType;
        /// Input type for the gradient
        typedef typename Functional_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Functional_::DataType DataType;

        /// Our baseclass
        typedef NLOptLS<Functional_, Filter_> BaseClass;
        /// Generic preconditioner
        typedef NLOptPrecond<typename Functional_::VectorTypeL, Filter_> PrecondType;

      protected:
        /// The linesearch used along the descent direction
        std::shared_ptr<LinesearchType> _linesearch;
        /// This will be the preconditioner, or a nullptr. We need to save it ourselves because we cannot access the
        /// prepare() routine through the SolverBase pointer in our BaseClass
        std::shared_ptr<PrecondType> _precond;

        /// defect vector
        VectorType _vec_r;
        /// descend direction vector
        VectorType _vec_p;

      public:
        /// For debugging purposes, all iterates can be logged to here
        std::deque<VectorType>* iterates;

      public:
        /**
         * \brief Standard constructor
         *
         * \param[in, out] functional_
         * The (nonlinear) functional. Cannot be const because it saves its own state
         *
         * \param[in] filter_
         * Filter to apply to the functional's gradient
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
        explicit NLSD(Functional_& functional_, Filter_& filter_,
        std::shared_ptr<LinesearchType> linesearch_,
        bool keep_iterates = false, std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass("NLSD", functional_, filter_, precond),
          _linesearch(linesearch_),
          _precond(precond),
          iterates(nullptr)
          {
            XASSERT(_linesearch != nullptr);

            this->set_ls_iter_digits(Math::ilog10(_linesearch->get_max_iter()));

            if(keep_iterates)
            {
              iterates = new std::deque<VectorType>;
            }
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
         * \param[in] functional_
         * The functional.
         *
         * \param[in] filter_
         * The system filter.
         *
         * \param[in] precond
         * The preconditioner. May be \c nullptr.
         *
         * \param[in] linesearch_
         * The linesearch to use.
         *
         */
        explicit NLSD(const String& section_name, PropertyMap* section,
        Functional_& functional_, Filter_& filter_, std::shared_ptr<LinesearchType> linesearch_,
        std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass("NLSD", section_name, section, functional_, filter_, precond),
          _linesearch(linesearch_),
          _precond(precond),
          iterates(nullptr)
          {
            XASSERT(_linesearch != nullptr);

            this->set_ls_iter_digits(Math::ilog10(_linesearch->get_max_iter()));

            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second && std::stoul(keep_iterates_p.first) == 1)
            {
              iterates = new std::deque<VectorType>;
            }
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~NLSD()
        {
          if(iterates != nullptr)
          {
            delete iterates;
          }
        }

        /// \copydoc SolverBase::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();
          // create three temporary vectors
          _vec_r = this->_functional.create_vector_r();
          _vec_p = this->_functional.create_vector_r();
          _linesearch->init_symbolic();
        }

        /// \copydoc SolverBase::done_symbolic()
        virtual void done_symbolic() override
        {
          if(iterates != nullptr)
            iterates->clear();

          //this->_vec_tmp.clear();
          this->_vec_p.clear();
          this->_vec_r.clear();
          _linesearch->done_symbolic();
          BaseClass::done_symbolic();
        }

        /// \copydoc SolverBase::name()
        virtual String name() const override
        {
          return "NLSD";
        }

        /**
         * \brief Sets the iterates deque according to a bool
         */
        void set_keep_iterates(bool keep_iterates)
        {
          if(iterates != nullptr)
          {
            delete iterates;
          }

          if(keep_iterates)
          {
            iterates = new std::deque<VectorType>;
          }
        }

        /// \copydoc SolverBase::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // Clear solution vector
          vec_cor.format();

          this->_functional.prepare(vec_cor, this->_filter);
          this->_functional.eval_fval_grad(this->_fval, this->_vec_r);

          // Copy back given defect
          this->_vec_r.copy(vec_def);

          if(this->_precond != nullptr)
          {
            this->_precond->prepare(vec_cor, this->_filter);
          }

          // apply
          this->_status = _apply_intern(vec_cor);
          this->plot_summary();
          return this->_status;
        }

        /// \copydoc IterativeSolver::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& DOXY(vec_rhs)) override
        {
          this->_functional.prepare(vec_sol, this->_filter);
          // Compute functional value and gradient
          this->_functional.eval_fval_grad(this->_fval, this->_vec_r);
          this->_vec_r.scale(this->_vec_r,DataType(-1));
          this->_filter.filter_def(this->_vec_r);

          if(this->_precond != nullptr)
          {
            this->_precond->prepare(vec_sol, this->_filter);
          }

          // apply
          this->_status = _apply_intern(vec_sol);
          this->plot_summary();
          return this->_status;
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
         * This does not have a right hand side because that is contained in the gradient of the functional and we
         * always seek grad functional(vec_sol) = 0
         *
         */
        virtual Status _apply_intern(VectorType& vec_sol)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

          // Reset member variables in the LineSearch
          _linesearch->reset();

          if(iterates != nullptr)
          {
            iterates->push_back(std::move(vec_sol.clone()));
          }

          // compute initial defect
          Status status = this->_set_initial_defect(this->_vec_r, vec_sol);
          if(status != Status::progress)
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }

          // The first direction has to be the steepest descent direction
          this->_vec_p.copy(this->_vec_r);

          // apply preconditioner to defect vector
          if(!this->_apply_precond(this->_vec_p, this->_vec_r, this->_filter))
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }

          // Compute initial eta = <d, r>
          DataType eta(this->_vec_p.dot(this->_vec_r));

          // If the preconditioned search direction is not a descent direction, reset it to steepest descent
          // TODO: Correct the output of the preconditioner if it turns out to not have been positive definite
          if(eta <= DataType(0))
          {
            this->_vec_p.copy(this->_vec_r);
          }

          // start iterating
          while(status == Status::progress)
          {
            IterationStats stat(*this);

            this->_fval_prev = this->_fval;

            // Copy information to the linesearch
            _linesearch->set_initial_fval(this->_fval);
            _linesearch->set_grad_from_defect(this->_vec_r);
            // If we are using a preconditioner, use the additional information about the preconditioned step length in the line search
            if(this->_precond != nullptr)
            {
              _linesearch->set_dir_scaling(true);
            }

            // Call the linesearch to update vec_sol
            status = _linesearch->correct(vec_sol, this->_vec_p);

            // Copy back information from the linesearch
            this->_fval = _linesearch->get_final_fval();
            this->_ls_its = _linesearch->get_num_iter();
            this->_steplength = _linesearch->get_rel_update();
            _linesearch->get_defect_from_grad(this->_vec_r);

            // Log iterates if necessary
            if(iterates != nullptr)
            {
              iterates->push_back(vec_sol.clone());
            }

            // Compute defect norm. This also performs the convergence/divergence checks.
            status = this->_set_new_defect(this->_vec_r, vec_sol);

            if(status != Status::progress)
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
              return status;
            }

            // Re-assemble preconditioner if necessary
            if(this->_precond != nullptr)
            {
              this->_precond->prepare(vec_sol, this->_filter);
            }

            // apply preconditioner
            if(!this->_apply_precond(_vec_p, _vec_r, this->_filter))
            {
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            }

            // Compute new eta
            eta = this->_vec_p.dot(this->_vec_r);

            // If the preconditioned search direction is not a descent direction, reset it to steepest descent
            // TODO: Correct the output of the preconditioner if it turns out to not have been positive definite
            if(eta <= DataType(0))
            {
              this->_vec_p.copy(this->_vec_r);
            }
          }

          // We should never come to this point
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
          return Status::undefined;
        }

    }; // class NLSD

    /**
     * \brief Creates a new NLSD solver object
     *
     * \param[in] functional
     * The nonlinear functional.
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
     * A shared pointer to a new NLSD object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Functional_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLSD<Functional_, Filter_>> new_nlsd(
      Functional_& functional, Filter_& filter, Linesearch_& linesearch, bool keep_iterates = false)
      {
        return std::make_shared<NLSD<Functional_, Filter_>>(functional, filter, linesearch,
        keep_iterates, nullptr);
      }
    template<typename Functional_, typename Filter_, typename Linesearch_, typename Precond_>
    inline std::shared_ptr<NLSD<Functional_, Filter_>> new_nlsd(
      Functional_& functional, Filter_& filter, Linesearch_& linesearch, bool keep_iterates,
      std::shared_ptr<Precond_> precond)
      {
        return std::make_shared<NLSD<Functional_, Filter_>>(functional, filter, linesearch,
        keep_iterates, precond);
      }
#else
    template<typename Functional_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLSD<Functional_, Filter_>> new_nlsd(
      Functional_& functional, Filter_& filter, Linesearch_& linesearch, bool keep_iterates = false,
      std::shared_ptr<NLOptPrecond<typename Functional_::VectorTypeL, Filter_>> precond = nullptr)
      {
        return std::make_shared<NLSD<Functional_, Filter_>>(functional, filter, linesearch,
        keep_iterates, precond);
      }
#endif

    /**
     * \brief Creates a new NLSD solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself.
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver.
     *
     * \param[in] functional
     * The nonlinear functional.
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
     * \returns
     * A shared pointer to a new NLSD object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Functional_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLSD<Functional_, Filter_>> new_nlsd(
      const String& section_name, PropertyMap* section,
      Functional_& functional, Filter_& filter, Linesearch_& linesearch)
      {
        return std::make_shared<NLSD<Functional_, Filter_>>(section_name, section, functional, filter, linesearch,
        nullptr);
      }

    template<typename Functional_, typename Filter_, typename Linesearch_, typename Precond_>
    inline std::shared_ptr<NLSD<Functional_, Filter_>> new_nlsd(
      const String& section_name, PropertyMap* section,
      Functional_& functional, Filter_& filter, Linesearch_& linesearch,
      std::shared_ptr<Precond_> precond)
      {
        return std::make_shared<NLSD<Functional_, Filter_>>(section_name, section, functional, filter, linesearch,
        precond);
      }
#else
    template<typename Functional_, typename Filter_, typename Linesearch_>
    inline std::shared_ptr<NLSD<Functional_, Filter_>> new_nlsd(
      const String& section_name, PropertyMap* section,
      Functional_& functional, Filter_& filter, Linesearch_& linesearch,
      std::shared_ptr<NLOptPrecond<typename Functional_::VectorTypeL, Filter_>> precond = nullptr)
      {
        return std::make_shared<NLSD<Functional_, Filter_>>(section_name, section, functional, filter, linesearch,
        precond);
      }
#endif
  } //namespace Solver
} // namespace FEAT

#endif
