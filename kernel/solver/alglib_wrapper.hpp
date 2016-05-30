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
#include <thirdparty/alglib/cpp/src/optimization.h>
FEAT_RESTORE_WARNINGS
namespace FEAT
{
  namespace Solver
  {
    namespace Intern
    {
      /// returns the object, if T_ has a GateType, i.e. is a GlobalVector - SFINAE at its best
      template <typename Evaluator_, typename T_>
      static auto derefer(T_ & object, typename Evaluator_::GateType *) -> decltype(*object)
      {
        return *object;
      }

      /// returns the dereferenced object, if T_ has no GateType, i.e. is a LocalVector - SFINAE at its best
      template <typename Evaluator_, typename T_>
      static T_& derefer(T_ & object, ...)
      {
        return object;
      }

    } // namespace Intern

    /**
     * \brief Wrapper around ALGLIB's lBFGS implementation for minimising an operator's gradient
     *
     * \tparam Operator_
     * Nonlinear Operator to minimise the gradient of
     *
     * \tparam Filter_
     * Filter to apply to the operator's gradient
     *
     * \note ALGLIB's algorithms always run in Mem::Main in double precision. Although the Operator can specify other
     * types, this will just to excess type conversions with no added benefit (like the speedup from computing in
     * float) and a general slow-down. It is not prohibited at this point so that these classes can be instantiated
     * so check if the implementation is clean.
     *
     */
    template<typename Operator_, typename Filter_>
    class ALGLIBMinLBFGS: public PreconditionedIterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        /// The nonlinear operator type
        typedef Operator_ OperatorType;
        /// The filter type
        typedef Filter_ FilterType;

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

      protected:
        /// Our nonlinear operator
        Operator_& _op;
        /// The filter we apply to the gradient
        Filter_& _filter;

        /// defect vector
        VectorType _vec_def;
        /// temporary vector
        VectorType _vec_tmp;

        /// Tolerance for function improvement
        DataType _tol_fval;
        /// Tolerance for gradient improvement
        DataType _tol_step;

        /// Optimisation variable for ALGLIB
        alglib::real_1d_array _opt_var;
        /// This will hold the state of the optimisation problem in ALGLIB
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
         * \param[in, out] op_
         * The (nonlinear) operator. Cannot be const because it saves its own state.
         *
         * \param[in] filter_
         * Filter to apply to the operator's gradient.
         *
         * \param[in] lbfgs_dim_
         * How many vectors to keep for the lBFGS hessian update. Defaults to min(7, op_.columns()).
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit ALGLIBMinLBFGS(
          Operator_& op_, Filter_& filter_, const alglib::ae_int_t lbfgs_dim_ = alglib::ae_int_t(0),
          const bool keep_iterates = false) :
          BaseClass("ALGLIBMinLBFGS"),
          _op(op_),
          _filter(filter_),
          _tol_fval(DataType(0)),
          _tol_step(Math::sqrt(Math::eps<DataType>())),
          _lbfgs_dim(lbfgs_dim_),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;

            if(_lbfgs_dim == alglib::ae_int_t(0))
              _lbfgs_dim = alglib::ae_int_t(Math::min(Index(7), _op.columns()));

            // ALGLIB's optimiser only know an absolute tolerance. We usually want to work with a relative one, which
            // is done through minlbfgsrequesttermination() from the _log callback function. So we set tol_abs to
            // zero, which can be overwritten with the set_tol_abs function of this class.
            set_tol_abs(DataType(0));
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~ALGLIBMinLBFGS()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();
          // create two temporary vectors
          _vec_def = this->_op.create_vector_r();
          _vec_tmp = this->_op.create_vector_r();

          // The length of the optimisation variable is the raw length of a temporary vector
          _opt_var.setlength(alglib::ae_int_t(
            Intern::derefer<VectorType>(_vec_def, nullptr).template size<LAFEM::Perspective::pod>()));

          for(alglib::ae_int_t i(0); i < _opt_var.length(); ++i)
            _opt_var[i] = double(0);

          alglib::minlbfgscreate(_lbfgs_dim, _opt_var, _state);
          alglib::minlbfgssetxrep(_state, true);
          // Set stopping criteria: absolute tolerance, function improvement, length of update step, max iterations
          alglib::minlbfgssetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
          alglib::ae_int_t(this->_max_iter));
          // Set the direction update
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          this->_vec_def.clear();
          this->_vec_tmp.clear();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "ALGLIBMinLBFGS";
        }

        /**
         * \copydoc BaseClass::set_max_iter()
         */
        void set_max_iter(Index max_iter)
        {
          this->_max_iter = max_iter;
          alglib::minlbfgssetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
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
        void set_tol_fval(DataType tol_fval)
        {
          _tol_fval = tol_fval;
          alglib::minlbfgssetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
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
          this->_tol_abs = tol_abs;
          alglib::minlbfgssetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
          alglib::ae_int_t(this->_max_iter));
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
          alglib::minlbfgssetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }

        /// \copydoc BaseClass::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // save defect
          this->_vec_def.copy(vec_def);
          //this->_system_filter.filter_def(this->_vec_def);

          // clear solution vector
          vec_cor.format();

          this->_op.prepare(vec_cor, this->_filter);

          // apply
          return _apply_intern(vec_cor);
        }

        /// \copydoc BaseClass::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& DOXY(vec_rhs)) override
        {
          this->_op.prepare(vec_sol, this->_filter);
          // compute defect
          this->_op.compute_grad(this->_vec_def);
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
         * \returns
         * A solver status code.
         *
         * This does not have a right hand side because that is contained in the gradient of the operator and we
         * always seek grad operator(vec_sol) = 0
         *
         */
        virtual Status _apply_intern(VectorType& vec_sol)
        {

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

          // Copy initial guess to optimisation state variable
          auto vec_sol_elements = Intern::derefer<VectorType>(vec_sol, nullptr).template elements<LAFEM::Perspective::pod>();
          for(alglib::ae_int_t i(0); i < _opt_var.length(); ++i)
            _opt_var[i] = vec_sol_elements[i];

          alglib::minlbfgsrestartfrom(_state, _opt_var);
          alglib::minlbfgsoptimize(_state, _func_grad, _log, this);
          alglib::minlbfgsresults(_state, _opt_var, _report);

          for(alglib::ae_int_t i(0); i < _opt_var.length(); ++i)
            vec_sol_elements[i] = DataType(_opt_var[i]);

          switch(_report.terminationtype)
          {
            case(-8):
              //std::cout << "ALGLIB: Got inf or NaN in function/gradient evaluation." << std::endl;
              return Status::aborted;
            case(-7):
              //std::cout << "ALGLIB: Gradient verification failed." << std::endl;
              return Status::aborted;
            case(1):
              //std::cout << "ALGLIB: Function value improvement criterion fulfilled." << std::endl;
              return Status::stagnated;
            case(2):
              //std::cout << "ALGLIB: Update step size stagnated." << std::endl;
              return Status::stagnated;
            case(4):
              //std::cout << "ALGLIB: Gradient norm criterion fulfilled." << std::endl;
              return Status::success;
            case(5):
              return Status::max_iter;
            case(7):
              //std::cout << "ALGLIB: Stopping criteria too stringent, further improvement impossible." << std::endl;
              return Status::stagnated;
            case(8):
              //std::cout << "ALGLIB: Stopped by user" << std::endl;
              return Status::success;
            default:
              return Status::undefined;
          }
        }

        /**
         * \brief Internal function for logging/plotting
         *
         * \param[in] x
         * The state variable of the optimisation problem
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
          ALGLIBMinLBFGS<OperatorType, FilterType>* me =
            reinterpret_cast<ALGLIBMinLBFGS<OperatorType, FilterType>*>(ptr);
          if(me->_num_iter>0)
          {
            // first, let's see if we have to compute the defect at all
            bool calc_def = false;
            calc_def = calc_def || (me->_min_iter < me->_max_iter);
            calc_def = calc_def || me->_plot;
            calc_def = calc_def || (me->_min_stag_iter > Index(0));

            // compute new defect
            if(calc_def)
              me->_def_cur = me->_calc_def_norm(me->_vec_def, me->_vec_tmp);

            Statistics::add_solver_defect(me->_branch, double(me->_def_cur));

            // plot?
            if(me->_plot)
            {
              std::cout << me->_plot_name
              <<  ": " << stringify(me->_num_iter).pad_front(me->_iter_digits)
              << " : " << stringify_fp_sci(me->_def_cur)
              << " / " << stringify_fp_sci(me->_def_cur / me->_def_init)
              << std::endl;
            }
            // Log iterates if necessary
            if(me->iterates != nullptr)
            {
              auto tmp = me->_vec_tmp.clone();
              me->iterates->push_back(std::move(tmp));
            }

            // Because ALGLIB knows no relative tolerance, we have to do it here
            if(me->_def_cur <= (me->_tol_rel * me->_def_init))
              alglib::minlbfgsrequesttermination(me->_state);

          }

          // increase iteration count
          ++me->_num_iter;
        }

        /**
         * \brief Internal function for computing functional value and gradient
         *
         * \param[in] x
         * The state variable of the optimisation problem
         *
         * \param[out] func
         * The functional value
         *
         * \param[out] grad
         * The operator's gradient
         *
         * \param[in] ptr
         * this, as the function needs to be static because it gets passed around as a C-style function pointer
         *
         */
        static void _func_grad(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)
        {
          // Downcast because we know what we are doing, right?
          ALGLIBMinLBFGS<OperatorType, FilterType>* me = reinterpret_cast<ALGLIBMinLBFGS<OperatorType, FilterType>*>
            (ptr);

          auto vec_tmp_elements = Intern::derefer<VectorType>(me->_vec_tmp, nullptr).template elements<LAFEM::Perspective::pod>();

          // Copy back ALGLIB's state variable to our solver
          for(alglib::ae_int_t i(0); i < x.length(); ++i)
            vec_tmp_elements[i] = DataType(x[i]);

          // Prepare the operator
          me->_op.prepare(me->_vec_tmp, me->_filter);
          // Compute functional value and gradient
          func = me->_op.compute_func();
          me->_op.compute_grad(me->_vec_def);
          me->_filter.filter_def(me->_vec_def);

          // Copy the operator's gradient to ALGLIB's grad variable
          auto vec_def_elements = Intern::derefer<VectorType>(me->_vec_def, nullptr).template elements<LAFEM::Perspective::pod>();
          for(alglib::ae_int_t i(0); i < grad.length(); ++i)
            grad[i] = double(vec_def_elements[i]);

        }

    }; // class ALGLIBMinLBFGS

    /**
     * \brief Creates a new ALGLIBMinLBFGS solver object
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new ALGLIBMinLBFGS object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<ALGLIBMinLBFGS<Operator_, Filter_>> new_alglib_minlbfgs(
      Operator_& op_, Filter_& filter_, alglib::ae_int_t lbfgs_dim_ = alglib::ae_int_t(0),
      bool keep_iterates_ = false)
      {
        return std::make_shared<ALGLIBMinLBFGS<Operator_, Filter_>>(op_, filter_, lbfgs_dim_, keep_iterates_);
      }
    /**
     * \brief Wrapper around ALGLIB's mincg implementation for minimising an operator's gradient
     *
     * \tparam Operator_
     * Nonlinear Operator to minimise the gradient of
     *
     * \tparam Filter_
     * Filter to apply to the operator's gradient
     *
     * \note ALGLIB's mincg appearantly supports preconditioning with the diagonal of the Hessian which can be set by
     * calling alglib::mincgsetprecdiag(state, [diagonal_of_hessian]). Every call like this triggers an internal
     * reset, seemingly turning the algorithm into steepest descent if called in every iteration. For strongly
     * nonlinear problems (like the ones we are interested in), this makes preconditioning awkward to useless.
     *
     * \note ALGLIB's algorithms always run in Mem::Main in double precision. Although the Operator can specify other
     * types, this will just to excess type conversions with no added benefit (like the speedup from computing in
     * float) and a general slow-down. It is not prohibited at this point so that these classes can be instantiated
     * so check if the implementation is clean.
     *
     */
    template<typename Operator_, typename Filter_>
    class ALGLIBMinCG : public PreconditionedIterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        /// The nonlinear operator type
        typedef Operator_ OperatorType;
        /// The filter type
        typedef Filter_ FilterType;

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

        /// Default NLCG search direction update
        static constexpr NLCGDirectionUpdate direction_update_default = NLCGDirectionUpdate::DYHSHybrid;

      protected:
        /// Our nonlinear operator
        Operator_& _op;
        /// The filter we apply to the gradient
        Filter_& _filter;

        /// Method to update the search direction, defaults to DYHSHybrid
        NLCGDirectionUpdate _direction_update;

        /// defect vector
        VectorType _vec_def;
        /// temporary vector
        VectorType _vec_tmp;

        /// Tolerance for function improvement
        DataType _tol_fval;
        /// Tolerance for gradient improvement
        DataType _tol_step;

        /// Optimisation variable for ALGLIB
        alglib::real_1d_array _opt_var;
        /// This will hold the state of the optimisation problem in ALGLIB
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
        explicit ALGLIBMinCG(Operator_& op_, Filter_& filter_,
        NLCGDirectionUpdate du_ = direction_update_default, bool keep_iterates = false) :
          BaseClass("ALGLIBMinCG"),
          _op(op_),
          _filter(filter_),
          _direction_update(du_),
          _tol_fval(DataType(0)),
          _tol_step(Math::sqrt(Math::eps<DataType>())),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;

            // ALGLIB's optimiser only know an absolute tolerance. We usually want to work with a relative one, which
            // is done through mincgrequesttermination() from the _log callback function. So we set tol_abs to zero,
            // which can be overwritten with the set_tol_abs function of this class.
            set_tol_abs(DataType(0));
          }

        /**
         * \brief Virtual destructor
         */
        virtual ~ALGLIBMinCG()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        /// \copydoc BaseClass::init_symbolic()
        virtual void init_symbolic() override
        {
          BaseClass::init_symbolic();
          // Create two temporary vectors
          _vec_def = this->_op.create_vector_r();
          _vec_tmp = this->_op.create_vector_r();

          // The length of the optimisation variable is the raw length of a temporary vector
          _opt_var.setlength(alglib::ae_int_t(
            Intern::derefer<VectorType>(_vec_def, nullptr).template size<LAFEM::Perspective::pod>()));

          for(alglib::ae_int_t i(0); i < _opt_var.length(); ++i)
            _opt_var[i] = double(0);

          alglib::mincgcreate(_opt_var, _state);
          alglib::mincgsetxrep(_state, true);
          // Set stopping criteria: Absolute tolerance, function improvement, length of update step, max iterations
          alglib::mincgsetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
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
                throw InternalError(__func__,__FILE__, __LINE__, name()+" got invalid direction update: "
                    +stringify(update_));

          }
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          this->_vec_def.clear();
          this->_vec_tmp.clear();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "ALGLIBMinCG-"+stringify(_direction_update);
        }

        /**
         * \copydoc BaseClass::set_max_iter()
         */
        void set_max_iter(Index max_iter)
        {
          this->_max_iter = max_iter;
          alglib::mincgsetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
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
        void set_tol_fval(DataType tol_fval)
        {
          _tol_fval = tol_fval;
          alglib::mincgsetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
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
          this->_tol_abs = tol_abs;
          alglib::mincgsetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }

        /**
         * \brief Sets the tolerance for the linesearch step size.
         *
         * \param[in] tol_step
         * New tolerance for the linesearch step size.
         *
         * If the linesearch fails to find a new iterate because its relative update is too small, the direction
         * update will fail to produce a new search direction so ALGLIBMinCG has to be terminated.
         *
         */
        void set_tol_step(DataType tol_step)
        {
          _tol_step = tol_step;
          alglib::mincgsetcond(_state, double(this->_tol_abs), double(_tol_fval), double(_tol_step),
          alglib::ae_int_t(this->_max_iter));
        }

        /// \copydoc BaseClass::apply()
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
        {
          // save defect
          this->_vec_def.copy(vec_def);
          //this->_system_filter.filter_def(this->_vec_def);

          // clear solution vector
          vec_cor.format();

          this->_op.prepare(vec_cor, this->_filter);

          // apply
          return _apply_intern(vec_cor);
        }

        /// \copydoc BaseClass::correct()
        virtual Status correct(VectorType& vec_sol, const VectorType& DOXY(vec_rhs)) override
        {
          this->_op.prepare(vec_sol, this->_filter);
          // compute defect
          this->_op.compute_grad(this->_vec_def);
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
         * \returns
         * A solver status code.
         *
         * This does not have a right hand side because that is contained in the gradient of the operator and we
         * always seek grad operator(vec_sol) = 0
         *
         */
        virtual Status _apply_intern(VectorType& vec_sol)
        {

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

          // Copy initial guess to optimisation state variable
          auto vec_sol_elements = Intern::derefer<VectorType>(vec_sol, nullptr).template elements<LAFEM::Perspective::pod>();
          for(alglib::ae_int_t i(0); i < _opt_var.length(); ++i)
            _opt_var[i] = vec_sol_elements[i];

          alglib::mincgrestartfrom(_state, _opt_var);
          alglib::mincgoptimize(_state, _func_grad, _log, this);
          alglib::mincgresults(_state, _opt_var, _report);

          for(alglib::ae_int_t i(0); i < _opt_var.length(); ++i)
            vec_sol_elements[i] = DataType(_opt_var[i]);

          switch(_report.terminationtype)
          {
            case(-8):
              //std::cout << "ALGLIB: Got inf or NaN in function/gradient evaluation." << std::endl;
              return Status::aborted;
            case(-7):
              //std::cout << "ALGLIB: Gradient verification failed." << std::endl;
              return Status::aborted;
            case(1):
              //std::cout << "ALGLIB: Function value improvement criterion fulfilled." << std::endl;
              return Status::stagnated;
            case(2):
              //std::cout << "ALGLIB: Update step size stagnated." << std::endl;
              return Status::stagnated;
            case(4):
              //std::cout << "ALGLIB: Gradient norm criterion fulfilled." << std::endl;
              return Status::success;
            case(5):
              return Status::max_iter;
            case(7):
              //std::cout << "ALGLIB: Stopping criteria too stringent, further improvement impossible." << std::endl;
              return Status::stagnated;
            case(8):
              //std::cout << "ALGLIB: Stopped by user" << std::endl;
              return Status::success;
            default:
              return Status::undefined;
          }
        }

        /**
         * \brief Internal function for logging/plotting
         *
         * \param[in] x
         * The state variable of the optimisation problem
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
          ALGLIBMinCG<OperatorType, FilterType>* me = reinterpret_cast<ALGLIBMinCG<OperatorType, FilterType>*>(ptr);

          if(me->_num_iter>0)
          {
            // first, let's see if we have to compute the defect at all
            bool calc_def = false;
            calc_def = calc_def || (me->_min_iter < me->_max_iter);
            calc_def = calc_def || me->_plot;
            calc_def = calc_def || (me->_min_stag_iter > Index(0));

            // compute new defect
            if(calc_def)
              me->_def_cur = me->_calc_def_norm(me->_vec_def, me->_vec_tmp);

            Statistics::add_solver_defect(me->_branch, double(me->_def_cur));

            // plot?
            if(me->_plot)
            {
              std::cout << me->_plot_name
              <<  ": " << stringify(me->_num_iter).pad_front(me->_iter_digits)
              << " : " << stringify_fp_sci(me->_def_cur)
              << " / " << stringify_fp_sci(me->_def_cur / me->_def_init)
              << std::endl;
            }
            // Log iterates if necessary
            if(me->iterates != nullptr)
            {
              auto tmp = me->_vec_tmp.clone();
              me->iterates->push_back(std::move(tmp));
            }

            // Because ALGLIB knows no relative tolerance, we have to do it here
            if(me->_def_cur <= (me->_tol_rel * me->_def_init))
              alglib::mincgrequesttermination(me->_state);

          }

          // increase iteration count
          ++me->_num_iter;

        }

        /**
         * \brief Internal function for computing functional value and gradient
         *
         * \param[in] x
         * The state variable of the optimisation problem
         *
         * \param[out] func
         * The functional value
         *
         * \param[out] grad
         * The operator's gradient
         *
         * \param[in] ptr
         * this, as the function needs to be static because it gets passed around as a C-style function pointer
         *
         */
        static void _func_grad(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)
        {
          // Downcast because we know what we are doing, right?
          ALGLIBMinCG<OperatorType, FilterType>* me = reinterpret_cast<ALGLIBMinCG<OperatorType, FilterType>*>(ptr);

          auto vec_tmp_elements = Intern::derefer<VectorType>(me->_vec_tmp, nullptr).template elements<LAFEM::Perspective::pod>();

          // Copy back ALGLIB's state variable to our solver
          for(alglib::ae_int_t i(0); i < x.length(); ++i)
            vec_tmp_elements[i] = DataType(x[i]);

          // Prepare the operator
          me->_op.prepare(me->_vec_tmp, me->_filter);
          // Compute functional value and gradient
          func = me->_op.compute_func();
          me->_op.compute_grad(me->_vec_def);
          me->_filter.filter_def(me->_vec_def);

          // Copy the operator's gradient to ALGLIB's grad variable
          auto vec_def_elements = Intern::derefer<VectorType>(me->_vec_def, nullptr).template elements<LAFEM::Perspective::pod>();

          for(alglib::ae_int_t i(0); i < grad.length(); ++i)
            grad[i] = double(vec_def_elements[i]);

        }

    }; // class ALGLIBMinCG

    /**
     * \brief Creates a new ALGLIBMinCG solver object
     *
     * \param[in] op
     * The operator
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false
     *
     * \returns
     * A shared pointer to a new ALGLIBMinCG object.
     */
    /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<ALGLIBMinCG<Operator_, Filter_>> new_alglib_mincg(
      Operator_& op_, Filter_& filter_,
      NLCGDirectionUpdate du_ = ALGLIBMinCG<Operator_, Filter_>::direction_update_default,
      bool keep_iterates_ = false)
      {
        return std::make_shared<ALGLIBMinCG<Operator_, Filter_>>(op_, filter_, du_, keep_iterates_);
      }


  } // namespace Solver
} // namespace FEAT
#endif // FEAT_HAVE_ALGLIB
#endif // FEAT_KERNEL_SOLVER_ALGLIB_WRAPPER
