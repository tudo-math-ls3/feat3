#pragma once
#ifndef FEAST_SOLVER_ALGLIB_WRAPPER
#define FEAST_SOLVER_ALGLIB_WRAPPER 1
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>

#include <deque>

#ifdef FEAST_HAVE_ALGLIB
FEAST_DISABLE_WARNINGS
#include <thirdparty/alglib/cpp/src/optimization.h>
FEAST_RESTORE_WARNINGS
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

      protected:
        /// Our nonlinear operator
        Operator_& _op;
        /// The filter we apply to the gradient
        const Filter_& _filter;

        /// defect vector
        VectorType _vec_def;
        /// temporary vector
        VectorType _vec_tmp;

        alglib::real_1d_array _opt_var;
        alglib::mincgstate _state;
        alglib::mincgreport _report;

      public:
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
        explicit ALGLIBMinCG(Operator_& op_, const Filter_& filter_, bool keep_iterates = false, std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass("ALGLIBMinCG", precond),
          _op(op_),
          _filter(filter_),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;
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
          double epsg(1e-10);
          double epsf(0);
          double epsx(0);
          BaseClass::init_symbolic();
          //// create three temporary vectors
          _vec_def = this->_op.create_vector_r();
          _vec_tmp = this->_op.create_vector_r();

          _opt_var.setlength(alglib::ae_int_t(_op.columns()));
          alglib::mincgcreate(_opt_var, _state);
          // Set the algorithm type: -1 (automatic selection of best algorithm), 0 (Dai-Yuan), 1 (hybrid Day-Yuan
          // and Hestenes-Stiefel)
          alglib::mincgsetcgtype(_state, -1);
          alglib::mincgsetcond(_state, epsg, epsf, epsx, alglib::ae_int_t(this->_max_iter));
        }

        /// \copydoc BaseClass::done_symbolic()
        virtual void done_symbolic() override
        {
          this->_vec_tmp.clear();
          this->_vec_def.clear();
          BaseClass::done_symbolic();
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "ALGLIBMinCG";
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

          auto vec_sol_elements = vec_sol.template elements<LAFEM::Perspective::pod>();
          for(alglib::ae_int_t i(0); i < _opt_var.length(); ++i)
            _opt_var[i] = vec_sol_elements[i];

          alglib::mincgoptimize(_state, _func_grad, nullptr, this);
          alglib::mincgresults(_state, _opt_var, _report);

        std::cout << "mincg: terminationtype " << _report.terminationtype << ", " << _report.iterationscount << " its, " << _report.nfev << " grad evals" << std::endl;

          for(alglib::ae_int_t i(0); i < _opt_var.length(); ++i)
            vec_sol_elements[i] = _opt_var[i];

          switch(_report.terminationtype)
          {
            case(-8):
              std::cout << "ALGLIB: Got inf or NaN in function/gradient evaluation." << std::endl;
              return Status::aborted;
            case(-7):
              std::cout << "ALGLIB: Gradient verification failed." << std::endl;
              return Status::aborted;
            case(1):
              std::cout << "ALGLIB: Function value improvement criterion fulfilled." << std::endl;
              return Status::stagnated;
            case(2):
              std::cout << "ALGLIB: Update step size stagnated." << std::endl;
              return Status::stagnated;
            case(4):
              std::cout << "ALGLIB: Gradient norm criterion fulfilled." << std::endl;
              return Status::success;
            case(5):
              return Status::max_iter;
            case(7):
              std::cout << "ALGLIB: Stopping criteria too stringent, further improvement impossible." << std::endl;
              return Status::stagnated;
            default:
              return Status::undefined;
          }
        }

        static void _log(const alglib::real_1d_array& DOXY(x), double DOXY(func), void* ptr)
        {
          std::cout << "log called" << std::endl;
          ALGLIBMinCG<OperatorType, FilterType>* me = reinterpret_cast<ALGLIBMinCG<OperatorType, FilterType>*>(ptr);

          // increase iteration count
          ++me->_num_iter;

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
        }

        static void _func_grad(const alglib::real_1d_array& x, double& func, alglib::real_1d_array& grad, void* ptr)
        {
          ALGLIBMinCG<OperatorType, FilterType>* me = reinterpret_cast<ALGLIBMinCG<OperatorType, FilterType>*>(ptr);

          auto vec_tmp_elements = me->_vec_tmp.template elements<LAFEM::Perspective::pod>();

          for(alglib::ae_int_t i(0); i < x.length(); ++i)
            vec_tmp_elements[i] = x[i];

          me->_op.prepare(me->_vec_tmp);
          func = me->_op.compute_fval();
          me->_op.compute_grad(me->_vec_def);
          me->_filter.filter_def(me->_vec_def);

          auto vec_def_elements = me->_vec_def.template elements<LAFEM::Perspective::pod>();
          for(alglib::ae_int_t i(0); i < grad.length(); ++i)
            grad[i] = vec_def_elements[i];

        }

    }; // class ALGLIBMinCG

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
    template<typename Operator_, typename Filter_>
    inline std::shared_ptr<ALGLIBMinCG<Operator_, Filter_>> new_alglib_mincg(
      Operator_& op, const Filter_& filter, bool keep_iterates = false)
      {
        return std::make_shared<ALGLIBMinCG<Operator_, Filter_>>(op, filter, keep_iterates, nullptr);
      }
#endif // FEAST_HAVE_ALGLIB
  } // namespace Solver
} // namespace FEAST
#endif // FEAST_KERNEL_SOLVER_ALGLIB_WRAPPER
