// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef FEAT_SOLVER_ALGLIB_WRAPPER
#define FEAT_SOLVER_ALGLIB_WRAPPER 1
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/solver/nlcg.hpp>

#include <deque>

#ifdef FEAT_HAVE_ALGLIB
FEAT_DISABLE_WARNINGS
#include <optimization.h>
FEAT_RESTORE_WARNINGS
#endif // FEAT_HAVE_ALGLIB

#if defined(FEAT_HAVE_ALGLIB) || defined(DOXYGEN)

namespace FEAT
{
  namespace Solver
  {
    /// \cond internal
    namespace Intern
    {
      /// returns the object, if T_ has a GateType, i.e. is a GlobalVector - SFINAE at its best
      template <typename Evaluator_, typename T_>
      auto derefer(T_ & object, typename Evaluator_::GateType *) -> decltype(object.local())
      {
        return object.local();
      }

      /// returns the dereferenced object, if T_ has no GateType, i.e. is a LocalVector - SFINAE at its best
      template <typename Evaluator_, typename T_>
      T_& derefer(T_ & object, ...)
      {
        return object;
      }
    } // namespace Intern
    /// \endcond

    /**
     * \brief Wrapper around ALGLIB's lBFGS implementation for minimising an functional's gradient
     *
     * \tparam Functional_
     * Nonlinear Functional to minimize the gradient of
     *
     * \tparam Filter_
     * Filter to apply to the functional's gradient
     *
     * \note ALGLIB's algorithms always run in Mem::Main in double precision. Although the Functional can specify other
     * types, this will just to excess type conversions with no added benefit (like the speedup from computing in
     * float) and a general slow-down. It is not prohibited at this point so that these classes can be instantiated
     * so check if the implementation is clean.
     *
     */
    template<typename Functional_, typename Filter_>
    class ALGLIBMinLBFGS: public NLOptLS<Functional_, Filter_>
    {
      public:
        /// The nonlinear functional type
        typedef Functional_ FunctionalType;
        /// The filter type
        typedef Filter_ FilterType;

        /// Type of the functional's gradient has
        typedef typename Functional_::GradientType GradientType;
        /// Input type for the gradient
        typedef typename Functional_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Functional_::DataType DataType;

        /// Our baseclass
        typedef NLOptLS<Functional_, Filter_> BaseClass;
        /// Generic preconditioner
        typedef SolverBase<VectorType> PrecondType;

      protected:

        /// defect vector
        VectorType _vec_def;
        /// temporary vector
        VectorType _vec_tmp;

        /// Optimization variable for ALGLIB
        alglib::real_1d_array _functionalt_var;
        /// This will hold the state of the optimization problem in ALGLIB
        alglib::minlbfgsstate _state;
        /// Convergence report etc.
        alglib::minlbfgsreport _report;
        /// Dimension for lBFGS Hessian update
        alglib::ae_int_t _lbfgs_dim;

      public:
        /// Can hold all iterates for debugging purposes
        std::deque<VectorType>* iterates;

      public:
        /**
         * \brief Standard constructor
         *
         * \param[in, out] functional_
         * The (nonlinear) functional. Cannot be const because it saves its own state.
         *
         * \param[in] filter_
         * Filter to apply to the functional's gradient.
         *
         * \param[in] lbfgs_dim_
         * How many vectors to keep for the lBFGS hessian update. Defaults to min(7, functional_.columns()).
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit ALGLIBMinLBFGS(
          Functional_& functional_, Filter_& filter_, const alglib::ae_int_t lbfgs_dim_ = alglib::ae_int_t(0),
          const bool keep_iterates = false) :
          BaseClass("ALGLIBMinLBFGS", functional_, filter_, nullptr),
          _lbfgs_dim(lbfgs_dim_),
          iterates(nullptr)
          {

            this->set_ls_iter_digits(2);

            if(keep_iterates)
            {
              iterates = new std::deque<VectorType>;
            }

            if(_lbfgs_dim == alglib::ae_int_t(0))
            {
              _lbfgs_dim = alglib::ae_int_t(Math::min(Index(7), this->_functional.columns()));
            }
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
         * \param[in] functional_
         * The functional.
         *
         * \param[in] filter_
         * The system filter.
         *
         */
        explicit ALGLIBMinLBFGS(const String& section_name, PropertyMap* section,
        Functional_& functional_, Filter_& filter_) :
          BaseClass("ALGLIBMinLBFGS", section_name, section, functional_, filter_, nullptr),
          _lbfgs_dim(0),
          iterates(nullptr)
          {

            this->set_ls_iter_digits(2);

            // Get direction update
            auto lbfgs_dim_p = section->query("lbfgs_dim");
            if(lbfgs_dim_p.second)
            {
              _lbfgs_dim = alglib::ae_int_t(std::stoul(lbfgs_dim_p.first));
              XASSERT(_lbfgs_dim > alglib::ae_int_t(0));
            }

            // Check if we have to keep the iterates
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second && std::stoul(keep_iterates_p.first) == 1)
            {
              iterates = new std::deque<VectorType>;
            }
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~ALGLIBMinLBFGS()
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
          // create two temporary vectors
          _vec_def = this->_functional.create_vector_r();
          _vec_tmp = this->_functional.create_vector_r();

          // The length of the optimization variable is the raw length of a temporary vector
          _functionalt_var.setlength(alglib::ae_int_t(
            Intern::derefer<VectorType>(_vec_def, nullptr).template size<LAFEM::Perspective::pod>()));

          for(alglib::ae_int_t i(0); i < _functionalt_var.length(); ++i)
          {
            _functionalt_var[i] = double(0);
          }

          alglib::minlbfgscreate(_lbfgs_dim, _functionalt_var, _state);
          alglib::minlbfgssetxrep(_state, true);
          // Set stopping criteria: absolute tolerance, function improvement, length of update step, max iterations
          // Since we do not want the solver to stop based on the absolute criterion alone, we always pass 0 as the
          // second argument
          alglib::minlbfgssetcond(_state, double(0), double(this->_tol_fval), double(this->_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }

        /// \copydoc SolverBase::done_symbolic()
        virtual void done_symbolic() override
        {
          this->_vec_def.clear();
          this->_vec_tmp.clear();
          BaseClass::done_symbolic();
        }

        /// \copydoc SolverBase::name()
        virtual String name() const override
        {
          return "ALGLIBMinLBFGS";
        }

        /**
         * \brief Sets the relative tolerance for the norm of the defect vector
         *
         * \param[in] tol_abs
         * New relative tolerance for the norm of the defect vector.
         *
         */
        void set_tol_abs(DataType tol_abs)
        {
          BaseClass::set_tol_abs(tol_abs);
          // Since we do not want the solver to stop based on the absolute criterion alone, we always pass 0 as the
          // second argument
          alglib::minlbfgssetcond(_state, double(0), double(this->_tol_fval), double(this->_tol_step),
          alglib::ae_int_t(this->_max_iter));
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
        virtual void set_tol_fval(DataType tol_fval) override
        {
          BaseClass::set_tol_fval(tol_fval);
          // Since we do not want the solver to stop based on the absolute criterion alone, we always pass 0 as the
          // second argument
          alglib::minlbfgssetcond(_state, double(0), double(this->_tol_fval), double(this->_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }

        /**
         * \copydoc IterativeSolver::set_max_iter()
         */
        void set_max_iter(Index max_iter)
        {
          BaseClass::set_max_iter(max_iter);
          // Since we do not want the solver to stop based on the absolute criterion alone, we always pass 0 as the
          // second argument
          alglib::minlbfgssetcond(_state, double(0), double(this->_tol_fval), double(this->_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }

        /**
         * \brief Sets the dimension for the LBFGS update
         */
        void set_lbfgs_dim(alglib::ae_int_t lbfgs_dim)
        {
          XASSERT(lbfgs_dim > alglib::ae_int_t(0));
          _lbfgs_dim = lbfgs_dim;
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

          // Unused here
          DataType fval(0);

          // clear solution vector
          vec_cor.format();

          this->_functional.prepare(vec_cor, this->_filter);
          this->_functional.eval_fval_grad(fval, this->_vec_def);

          // Copy back defect
          this->_vec_def.copy(vec_def);
          this->_filter.filter_def(this->_vec_def);

          // apply
          this->_status = _apply_intern(vec_cor);
          this->plot_summary();
          return this->_status;
        }

        /// \copydoc IterativeSolver::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& DOXY(vec_rhs)) override
        {

          this->_functional.prepare(vec_sol, this->_filter);
          // compute defect
          this->_functional.eval_fval_grad(this->_fval, this->_vec_def);
          this->_vec_def.scale(this->_vec_def,DataType(-1));
          this->_filter.filter_def(this->_vec_def);

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

          // Write initial guess to iterates if desired
          if(iterates != nullptr)
          {
            auto tmp = vec_sol.clone();
            iterates->push_back(std::move(tmp));
          }

          // compute initial defect
          Status status = this->_set_initial_defect(this->_vec_def, vec_sol);

          if(status != Status::progress)
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }

          // Copy initial guess to optimization state variable
          auto vec_sol_elements = Intern::derefer<VectorType>(vec_sol, nullptr).template elements<LAFEM::Perspective::pod>();
          for(alglib::ae_int_t i(0); i < _functionalt_var.length(); ++i)
          {
            _functionalt_var[i] = vec_sol_elements[i];
          }

          alglib::minlbfgsrestartfrom(_state, _functionalt_var);

          IterationStats stat(*this);
          alglib::minlbfgsoptimize(_state, _func_grad, _log, this);
          alglib::minlbfgsresults(_state, _functionalt_var, _report);

          for(alglib::ae_int_t i(0); i < _functionalt_var.length(); ++i)
          {
            vec_sol_elements[i] = DataType(_functionalt_var[i]);
          }

          switch(_report.terminationtype)
          {
            case(-8):
              //std::cout << "ALGLIB: Got inf or NaN in function/gradient evaluation." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            case(-7):
              //std::cout << "ALGLIB: Gradient verification failed." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            case(1):
              //std::cout << "ALGLIB: Function value improvement criterion fulfilled." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            case(2):
              //std::cout << "ALGLIB: Update step size stagnated." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            case(4):
              //std::cout << "ALGLIB: Gradient norm criterion fulfilled." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            case(5):
              //std::cout << "ALGLIB: Maximum number of iterations" << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::max_iter, this->get_num_iter()));
              return Status::max_iter;
            case(7):
              //std::cout << "ALGLIB: Stopping criteria too stringent, further improvement impossible." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::stagnated, this->get_num_iter()));
              return Status::stagnated;
            case(8):
              //std::cout << "ALGLIB: Stopped by user" << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            default:
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
              return Status::undefined;
          }
        }

        /**
         * \brief Internal function for logging/plotting
         *
         * \param[in] x
         * The state variable of the optimization problem
         *
         * \param[in] func
         * The functional value
         *
         * \param[in] ptr
         * this, as the function needs to be static because it gets passed around as a C-style function pointer
         *
         */
        static void _log(const alglib::real_1d_array& DOXY(x), double DOXY(func), void* ptr)
        {
          ALGLIBMinLBFGS<FunctionalType, FilterType>* me =
            reinterpret_cast<ALGLIBMinLBFGS<FunctionalType, FilterType>*>(ptr);

          // Because of how ALGLIB counts its iterations, we have to make sure we do not call this before the first
          // functional evaluation (repnfev is the reported number of functional evaluations)
          if(me->_state.c_ptr()->repnfev> 0)
          {
            me->_steplength = DataType(me->_state.c_ptr()->stp);
            me->_ls_its = Index(me->_state.c_ptr()->nfev);
            Status st = me->_set_new_defect(me->_vec_def, me->_vec_tmp);
            // Because ALGLIB knows no relative tolerance, we have to do it here
            if( st != Status::progress  )
            {
              alglib::minlbfgsrequesttermination(me->_state);
            }
          }

        }

        /**
         * \brief Internal function for computing functional value and gradient
         *
         * \param[in] x
         * The state variable of the optimization problem
         *
         * \param[out] fval
         * The functional value
         *
         * \param[out] grad
         * The functional's gradient
         *
         * \param[in] ptr
         * this, as the function needs to be static because it gets passed around as a C-style function pointer
         *
         */
        static void _func_grad(const alglib::real_1d_array& x, double& fval, alglib::real_1d_array& grad, void* ptr)
        {
          // Downcast because we know what we are doing, right?
          ALGLIBMinLBFGS<FunctionalType, FilterType>* me = reinterpret_cast<ALGLIBMinLBFGS<FunctionalType, FilterType>*>
            (ptr);

          auto vec_tmp_elements = Intern::derefer<VectorType>(me->_vec_tmp, nullptr).template elements<LAFEM::Perspective::pod>();

          // Copy back ALGLIB's state variable to our solver
          for(alglib::ae_int_t i(0); i < x.length(); ++i)
            vec_tmp_elements[i] = DataType(x[i]);

          me->_fval_prev = me->_fval;
          // Prepare the functional
          me->_functional.prepare(me->_vec_tmp, me->_filter);
          // Compute functional value and gradient
          me->_functional.eval_fval_grad(me->_fval, me->_vec_def);
          me->_filter.filter_def(me->_vec_def);

          fval = double(me->_fval);

          // Copy the functional's gradient to ALGLIB's grad variable
          auto vec_def_elements = Intern::derefer<VectorType>(me->_vec_def, nullptr).template elements<LAFEM::Perspective::pod>();
          for(alglib::ae_int_t i(0); i < grad.length(); ++i)
          {
            grad[i] = double(vec_def_elements[i]);
          }

        }

    }; // class ALGLIBMinLBFGS

    /**
     * \brief Creates a new ALGLIBMinLBFGS solver object
     *
     * \param[in] functional_
     * The functional
     *
     * \param[in] filter_
     * The system filter.
     *
     * \param[in] lbfgs_dim_
     * How many vectors to keep for the lBFGS hessian update. Defaults to min(7, functional_.columns()).
     *
     * \param[in] keep_iterates
     * Keep all iterates in a std::deque. Defaults to false.
     *
     * \returns
     * A shared pointer to a new ALGLIBMinLBFGS object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<ALGLIBMinLBFGS<Functional_, Filter_>> new_alglib_minlbfgs(
      Functional_& functional_, Filter_& filter_,
      alglib::ae_int_t lbfgs_dim_ = alglib::ae_int_t(0), bool keep_iterates = false)
      {
        return std::make_shared<ALGLIBMinLBFGS<Functional_, Filter_>>(functional_, filter_, lbfgs_dim_, keep_iterates);
      }

    /**
     * \brief Creates a new ALGLIBMinLBFGS solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] functional_
     * The functional
     *
     * \param[in] filter_
     * The system filter.
     *
     * \returns
     * A shared pointer to a new ALGLIBMinLBFGS object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<ALGLIBMinLBFGS<Functional_, Filter_>> new_alglib_minlbfgs(
      const String& section_name, PropertyMap* section,
      Functional_& functional_, Filter_& filter_)
      {
        return std::make_shared<ALGLIBMinLBFGS<Functional_, Filter_>>(section_name, section, functional_, filter_);
      }

    /**
     * \brief Wrapper around ALGLIB's mincg implementation for minimising an functional's gradient
     *
     * \tparam Functional_
     * Nonlinear Functional to minimize the gradient of
     *
     * \tparam Filter_
     * Filter to apply to the functional's gradient
     *
     * \note ALGLIB's mincg appearantly supports preconditioning with the diagonal of the Hessian which can be set by
     * calling alglib::mincgsetprecdiag(state, [diagonal_of_hessian]). Every call like this triggers an internal
     * reset, seemingly turning the algorithm into steepest descent if called in every iteration. For strongly
     * nonlinear problems (like the ones we are interested in), this makes preconditioning awkward to useless.
     *
     * \note ALGLIB's algorithms always run in Mem::Main in double precision. Although the Functional can specify other
     * types, this will just to excess type conversions with no added benefit (like the speedup from computing in
     * float) and a general slow-down. It is not prohibited at this point so that these classes can be instantiated
     * so check if the implementation is clean.
     *
     */
    template<typename Functional_, typename Filter_>
    class ALGLIBMinCG :  public NLOptLS<Functional_, Filter_>
    {
      public:
        /// The nonlinear functional type
        typedef Functional_ FunctionalType;
        /// The filter type
        typedef Filter_ FilterType;

        /// Type of the functional's gradient has
        typedef typename Functional_::GradientType GradientType;
        /// Input type for the gradient
        typedef typename Functional_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Functional_::DataType DataType;

        /// Our baseclass
        typedef NLOptLS<Functional_, Filter_> BaseClass;
        /// Generic preconditioner
        typedef SolverBase<VectorType> PrecondType;

        /// Default NLCG search direction update
        static constexpr NLCGDirectionUpdate direction_update_default = NLCGDirectionUpdate::DYHSHybrid;

      protected:
        /// Method to update the search direction, defaults to DYHSHybrid
        NLCGDirectionUpdate _direction_update;

        /// defect vector
        VectorType _vec_def;
        /// temporary vector
        VectorType _vec_tmp;

        /// Optimization variable for ALGLIB
        alglib::real_1d_array _functionalt_var;
        /// This will hold the state of the optimization problem in ALGLIB
        alglib::mincgstate _state;
        /// Convergence report etc.
        alglib::mincgreport _report;

      public:
        /// Can hold all iterates for debugging purposes
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
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit ALGLIBMinCG(Functional_& functional_, Filter_& filter_,
        NLCGDirectionUpdate du_ = direction_update_default, bool keep_iterates = false) :
          BaseClass("ALGLIBMinCG", functional_, filter_, nullptr),
          _direction_update(du_),
          iterates(nullptr)
          {

            this->set_ls_iter_digits(2);

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
         * The functional
         *
         * \param[in] filter_
         * The system filter.
         *
         */
        explicit ALGLIBMinCG(const String& section_name, PropertyMap* section,
        Functional_& functional_, Filter_& filter_) :
          BaseClass("ALGLIBMinCG", section_name, section, functional_, filter_, nullptr),
          _direction_update(direction_update_default),
          iterates(nullptr)
          {

            this->set_ls_iter_digits(2);

            // Get direction update
            auto direction_update_p = section->query("direction_update");
            if(direction_update_p.second)
            {
              NLCGDirectionUpdate my_update;
              my_update << direction_update_p.first;
              set_direction_update(my_update);
            }

            // Check if we have to keep the iterates
            auto keep_iterates_p = section->query("keep_iterates");
            if(keep_iterates_p.second && std::stoul(keep_iterates_p.first) == 1)
            {
              iterates = new std::deque<VectorType>;
            }
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~ALGLIBMinCG()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc SolverBase::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();
          // Create two temporary vectors
          _vec_def = this->_functional.create_vector_r();
          _vec_tmp = this->_functional.create_vector_r();

          // The length of the optimization variable is the raw length of a temporary vector
          _functionalt_var.setlength(alglib::ae_int_t(
            Intern::derefer<VectorType>(_vec_def, nullptr).template size<LAFEM::Perspective::pod>()));

          for(alglib::ae_int_t i(0); i < _functionalt_var.length(); ++i)
          {
            _functionalt_var[i] = double(0);
          }

          alglib::mincgcreate(_functionalt_var, _state);
          alglib::mincgsetxrep(_state, true);
          // Set stopping criteria: Absolute tolerance, function improvement, length of update step, max iterations
          // Since we do not want the solver to stop based on the absolute criterion alone, we always pass 0 as the
          // second argument
          alglib::mincgsetcond(_state, double(0), double(this->_tol_fval), double(this->_tol_step),
          alglib::ae_int_t(this->_max_iter));
          // Set the direction update
          set_direction_update(this->_direction_update);
        }

        /**
         * \brief Sets the direction update method
         */
        void set_direction_update(NLCGDirectionUpdate update_)
        {
          switch(update_)
          {
            case NLCGDirectionUpdate::DaiYuan:
              alglib::mincgsetcgtype(_state, 0);
              break;
            case NLCGDirectionUpdate::DYHSHybrid:
              alglib::mincgsetcgtype(_state, 1);
              break;
            default:
              XABORTM("got invalid direction update: " +stringify(update_));
          }
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

        /**
         * \brief Sets the tolerance for function value improvement
         *
         * \param[in] tol_fval
         * New tolerance for function value improvement.
         *
         * The convergence check is against the maximum of the absolute and relative function value.
         *
         */
        virtual void set_tol_fval(DataType tol_fval) override
        {
          BaseClass::set_tol_fval(tol_fval);
          // Since we do not want the solver to stop based on the absolute criterion alone, we always pass 0 as the
          // second argument
          alglib::mincgsetcond(_state, double(0), double(this->_tol_fval), double(this->_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }

        /**
         * \copydoc IterativeSolver::set_max_iter()
         */
        void set_max_iter(Index max_iter)
        {
          BaseClass::set_max_iter(max_iter);
          // Since we do not want the solver to stop based on the absolute criterion alone, we always pass 0 as the
          // second argument
          alglib::mincgsetcond(_state, double(0), double(this->_tol_fval), double(this->_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }

        /**
         * \brief Sets the relative tolerance for the norm of the defect vector
         *
         * \param[in] tol_abs
         * New relative tolerance for the norm of the defect vector.
         *
         */
        void set_tol_abs(DataType tol_abs)
        {
          BaseClass::set_tol_abs(tol_abs);
          // Since we do not want the solver to stop based on the absolute criterion alone, we always pass 0 as the
          // second argument
          alglib::mincgsetcond(_state, double(0), double(this->_tol_fval), double(this->_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }


        /// \copydoc SolverBase::done_symbolic()
        virtual void done_symbolic() override
        {
          this->_vec_def.clear();
          this->_vec_tmp.clear();
          BaseClass::done_symbolic();
        }

        /// \copydoc SolverBase::name()
        virtual String name() const override
        {
          return "ALGLIBMinCG";
        }

        /// \copydoc SolverBase::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // clear solution vector
          vec_cor.format();

          this->_functional.prepare(vec_cor, this->_filter);
          this->_functional.eval_fval_grad(this->_fval, this->_vec_def);

          // Copy back defect
          this->_vec_def.copy(vec_def);
          this->_filter.filter_def(this->_vec_def);

          // apply
          this->_status = _apply_intern(vec_cor);
          this->plot_summary();
          return this->_status;
        }

        /// \copydoc IterativeSolver::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& DOXY(vec_rhs)) override
        {

          this->_functional.prepare(vec_sol, this->_filter);
          // compute defect
          this->_functional.eval_fval_grad(this->_fval, this->_vec_def);
          this->_vec_def.scale(this->_vec_def,DataType(-1));
          this->_filter.filter_def(this->_vec_def);

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

          // Write initial guess to iterates if desired
          if(iterates != nullptr)
          {
            auto tmp = vec_sol.clone();
            iterates->push_back(std::move(tmp));
          }
          // compute initial defect
          Status status = this->_set_initial_defect(this->_vec_def, vec_sol);

          if(status != Status::progress)
            return status;

          // Copy initial guess to optimization state variable
          auto vec_sol_elements = Intern::derefer<VectorType>(vec_sol, nullptr).template elements<LAFEM::Perspective::pod>();
          for(alglib::ae_int_t i(0); i < _functionalt_var.length(); ++i)
          {
            _functionalt_var[i] = vec_sol_elements[i];
          }

          alglib::mincgrestartfrom(_state, _functionalt_var);
          IterationStats stat(*this);
          alglib::mincgoptimize(_state, _func_grad, _log, this);
          alglib::mincgresults(_state, _functionalt_var, _report);

          for(alglib::ae_int_t i(0); i < _functionalt_var.length(); ++i)
          {
            vec_sol_elements[i] = DataType(_functionalt_var[i]);
          }

          switch(_report.terminationtype)
          {
            case(-8):
              //std::cout << "ALGLIB: Got inf or NaN in function/gradient evaluation." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            case(-7):
              //std::cout << "ALGLIB: Gradient verification failed." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
              return Status::aborted;
            case(1):
              //std::cout << "ALGLIB: Function value improvement criterion fulfilled." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            case(2):
              //std::cout << "ALGLIB: Update step size stagnated." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            case(4):
              //std::cout << "ALGLIB: Gradient norm criterion fulfilled." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            case(5):
              //std::cout << "ALGLIB: Maximum number of iterations" << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::max_iter, this->get_num_iter()));
              return Status::max_iter;
            case(7):
              //std::cout << "ALGLIB: Stopping criteria too stringent, further improvement impossible." << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::stagnated, this->get_num_iter()));
              return Status::stagnated;
            case(8):
              //std::cout << "ALGLIB: Stopped by user" << std::endl;
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::success, this->get_num_iter()));
              return Status::success;
            default:
              Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::undefined, this->get_num_iter()));
              return Status::undefined;
          }
        }

        /**
         * \brief Internal function for logging/plotting
         *
         * \param[in] x
         * The state variable of the optimization problem
         *
         * \param[in] func
         * The functional value
         *
         * \param[in] ptr
         * this, as the function needs to be static because it gets passed around as a C-style function pointer
         *
         */
        static void _log(const alglib::real_1d_array& DOXY(x), double DOXY(func), void* ptr)
        {
          ALGLIBMinCG<FunctionalType, FilterType>* me = reinterpret_cast<ALGLIBMinCG<FunctionalType, FilterType>*>(ptr);

          // Because of how ALGLIB counts its iterations, we have to make sure we do not call this before the first
          // functional evaluation (repnfev is the reported number of functional evaluations)
          if(me->_state.c_ptr()->repnfev > 0)
          {
            me->_steplength = DataType(me->_state.c_ptr()->stp);
            me->_ls_its = Index(me->_state.c_ptr()->nfev);
            Status st = me->_set_new_defect(me->_vec_def, me->_vec_tmp);
            // Because ALGLIB knows no relative tolerance, we have to do it here
            if( st != Status::progress  )
            {
              alglib::mincgrequesttermination(me->_state);
            }
          }

        }

        /**
         * \brief Internal function for computing functional value and gradient
         *
         * \param[in] x
         * The state variable of the optimization problem
         *
         * \param[out] fval
         * The functional value
         *
         * \param[out] grad
         * The functional's gradient
         *
         * \param[in] ptr
         * this, as the function needs to be static because it gets passed around as a C-style function pointer
         *
         */
        static void _func_grad(const alglib::real_1d_array& x, double& fval, alglib::real_1d_array& grad, void* ptr)
        {
          // Downcast because we know what we are doing, right?
          ALGLIBMinCG<FunctionalType, FilterType>* me = reinterpret_cast<ALGLIBMinCG<FunctionalType, FilterType>*>(ptr);

          auto vec_tmp_elements = Intern::derefer<VectorType>(me->_vec_tmp, nullptr).template
            elements<LAFEM::Perspective::pod>();

          // Copy back ALGLIB's state variable to our solver
          for(alglib::ae_int_t i(0); i < x.length(); ++i)
          {
            vec_tmp_elements[i] = DataType(x[i]);
          }

          me->_fval_prev = me->_fval;
          // Prepare the functional
          me->_functional.prepare(me->_vec_tmp, me->_filter);
          // Compute functional value and gradient
          me->_functional.eval_fval_grad(me->_fval, me->_vec_def);
          me->_filter.filter_def(me->_vec_def);

          fval = double(me->_fval);

          // Copy the functional's gradient to ALGLIB's grad variable
          auto vec_def_elements = Intern::derefer<VectorType>(me->_vec_def, nullptr).template
            elements<LAFEM::Perspective::pod>();

          for(alglib::ae_int_t i(0); i < grad.length(); ++i)
          {
            grad[i] = double(vec_def_elements[i]);
          }

        }

    }; // class ALGLIBMinCG

    /**
     * \brief Creates a new ALGLIBMinCG solver object
     *
     * \param[in] functional_
     * The functional
     *
     * \param[in] filter_
     * The system filter.
     *
     * \param[in] keep_iterates_
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new ALGLIBMinCG object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<ALGLIBMinCG<Functional_, Filter_>> new_alglib_mincg(
      Functional_& functional_, Filter_& filter_,
      NLCGDirectionUpdate du_ = ALGLIBMinCG<Functional_, Filter_>::direction_update_default,
      bool keep_iterates_ = false)
    {
      return std::make_shared<ALGLIBMinCG<Functional_, Filter_>>(functional_, filter_, du_, keep_iterates_);
    }

    /**
     * \brief Creates a new ALGLIBMinCG solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] functional_
     * The functional
     *
     * \param[in] filter_
     * The system filter.
     *
     * \returns
     * A shared pointer to a new ALGLIBMinCG object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<ALGLIBMinCG<Functional_, Filter_>> new_alglib_mincg(
      const String& section_name, PropertyMap* section,
      Functional_& functional_, Filter_& filter_)
    {
      return std::make_shared<ALGLIBMinCG<Functional_, Filter_>>(section_name, section, functional_, filter_);
    }
  } // namespace Solver
} // namespace FEAT
#endif // defined(FEAT_HAVE_ALGLIB) || defined(DOXYGEN)
#endif // FEAT_KERNEL_SOLVER_ALGLIB_WRAPPER
