#pragma once
#ifndef FEAT_SOLVER_NLOPTLS
#define FEAT_SOLVER_NLOPTLS 1
#include <kernel/base_header.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/solver/linesearch.hpp>
#include <kernel/solver/nlopt_precond.hpp>

namespace FEAT
{
  namespace Solver
  {

    /**
     * \brief Base class for line search based nonlinear optimisers
     *
     * \tparam Operator_
     * The nonlinear functional to be minimised
     *
     * \tparam Filter_
     * The filters representing essential (boundary) conditions for the functional's state
     *
     * This is the base class for all nonlinear optimisers which use line searches. These solvers try to minimise a
     * nonlinear functional (given by Op_) by finding critical points \f$x \f$, meaning
     * \f$ \mathrm{grad} \mathcal{F}(x) = 0.\f$
     * The line search methods can be based on Solver::Linesearch, or they can be implicit like the ones used by
     * Solver::ALGLIBMinLBFGS.
     *
     * Due to using a line search, solvers derived from this have additional stopping criteria based on
     * - the functional value improvement
     * - the step size of the line search
     * - the line search not terminating successfully
     *
     * \author Jordi Paul
     *
     */
    template<typename Operator_, typename Filter_>
    class NLOptLS : public PreconditionedIterativeSolver<typename Operator_::VectorTypeR>
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

      protected:
        /// Our nonlinear operator
        Operator_& _op;
        /// The filter we apply to the gradient
        Filter_& _filter;

        /// Tolerance for function improvement
        DataType _tol_fval;
        /// Tolerance for the length of the update step
        DataType _tol_step;

        /// Initial function value
        DataType _fval_init;
        /// Current functional value
        DataType _fval;
        /// Functional value from the previous iteration
        DataType _fval_prev;
        /// The last step length of the line search
        DataType _steplength;
        /// The number of digits used to plot line search iteration numbers
        Index _ls_iter_digits;
        /// The number of iterations the linesearch performed in the last iteration of this solver
        Index _ls_its;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] plot_name_
         * The String identifier used for plots.
         *
         * \param[in] op_
         * The nonlinear functional. Cannot be const as internal variables change upon functional evaluation.
         *
         * \param[in] filter_
         * The filter for essential boundary conditions. Cannot be const, see op_.
         *
         * \param[in] precond
         * The preconditioner, defaults to nullptr. Cannot be const, see op_.
         *
         */
        explicit NLOptLS(const String& plot_name_, Operator_& op_, Filter_& filter_,
        std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass(plot_name_, precond),
          _op(op_),
          _filter(filter_),
          _tol_fval(DataType(0)),
          _tol_step(Math::eps<DataType>()),
          _fval_init(-Math::huge<DataType>()),
          _fval(Math::huge<DataType>()),
          _fval_prev(-Math::huge<DataType>()),
          _steplength(1),
          _ls_iter_digits(2),
          _ls_its(~Index(0))
          {
          }

        /**
         * \brief Constructor
         *
         * \param[in] plot_name_
         * The String identifier used for plots.
         *
         * \param[in] op_
         * The nonlinear functional. Cannot be const as internal variables change upon functional evaluation.
         *
         * \param[in] filter_
         * The filter for essential boundary conditions. Cannot be const, see op_.
         *
         * \param[in] precond
         * The preconditioner, defaults to nullptr. Cannot be const, see op_.
         *
         */
        explicit NLOptLS(const String& plot_name, const String& section_name, PropertyMap* section,
        Operator_& op_, Filter_& filter_, std::shared_ptr<PrecondType> precond = nullptr) :
          BaseClass(plot_name, section_name, section, precond),
          _op(op_),
          _filter(filter_),
          _tol_fval(DataType(0)),
          _tol_step(Math::eps<DataType>()),
          _fval_init(-Math::huge<DataType>()),
          _fval(Math::huge<DataType>()),
          _fval_prev(-Math::huge<DataType>()),
          _steplength(1),
          _ls_iter_digits(2),
          _ls_its(~Index(0))
          {
            auto tol_fval_p = section->get_entry("tol_fval");
            if (tol_fval_p.second)
            {
              set_tol_fval((DataType)std::stod(tol_fval_p.first));
            }

            auto tol_step_p = section->get_entry("tol_step");
            if (tol_step_p.second)
            {
              set_tol_step((DataType)std::stod(tol_step_p.first));
            }
          }

        /**
         * \brief Empty virtual destructor
         */
        virtual ~NLOptLS()
        {
        }

        /**
         * \brief Plot a summary of the last solver run
         */
        virtual void plot_summary(Status st) const override
        {
          // Print solver summary
          Dist::Comm comm_world(Dist::Comm::world());

          String msg(this->get_plot_name()+ ": its: "+stringify(this->get_num_iter())+" ("+ stringify(st)+")"
              +", evals: "+stringify(_op.get_num_func_evals())+" (func) "
              + stringify(_op.get_num_grad_evals()) + " (grad) "
              + stringify(_op.get_num_hess_evals()) + " (hess)"
              +" last step: "+stringify_fp_sci(_steplength)+"\n");
          msg +=this->get_plot_name()+": fval: "+stringify_fp_sci(_fval_init)
            + " -> "+stringify_fp_sci(_fval)
            + ", reduction factor "+stringify_fp_sci(_fval/_fval_init)+"\n";
          msg += this->get_plot_name()  +": grad: "+stringify_fp_sci(this->_def_init)
            + " -> "+stringify_fp_sci(this->_def_cur)
            + ", reduction factor " +stringify_fp_sci(this->_def_cur/this->_def_init);
          comm_world.print(msg);

        }

        /// \copydoc BaseClass::write_config()
        virtual PropertyMap* write_config(PropertyMap* parent, const String& section_name) const override
        {
          XASSERT(parent != nullptr);

          Dist::Comm comm(Dist::Comm::world());

          PropertyMap* my_section = BaseClass::write_config(parent, section_name);

          my_section->add_entry("tol_fval", stringify_fp_sci(_tol_fval));
          my_section->add_entry("tol_step", stringify_fp_sci(_tol_step));

          return my_section;

        }

        /**
         * \brief Gets the tolerance for function value improvement
         *
         * The convergence check is against the maximum of the absolute and relative function value.
         *
         * \returns The function value improvement tolerance
         */
        virtual DataType get_tol_fval()
        {
          return _tol_fval;
        }

        /**
         * \brief Gets the tolerance for the linesearch step size
         *
         * If the linesearch fails to find a new iterate because its relative update is too small, the direction
         * update will fail to produce a new search direction so the NLSD has to be terminated.
         *
         * \returns The function value improvement tolerance
         */
        virtual DataType get_tol_step()
        {
          return _tol_step;
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
        virtual void set_tol_fval(DataType tol_fval)
        {
          _tol_fval = tol_fval;
        }

        /**
         * \brief Sets the tolerance for the linesearch step size
         *
         * \param[in] tol_step
         * New tolerance for the linesearch step size.
         *
         * If the linesearch fails to find a new iterate because its relative update is too small, the direction
         * update will fail to produce a new search direction so the NLCG has to be terminated.
         *
         */
        virtual void set_tol_step(DataType tol_step)
        {
          _tol_step = tol_step;
        }

        /**
         * \brief Sets the number of digits used to plot line search iteration numbers
         */
        virtual void set_ls_iter_digits(Index digits)
        {
          _ls_iter_digits = digits;
        }

      protected:
        /**
         * \copydoc BaseClass::_set_initial_defect()
         */
        virtual Status _set_initial_defect(const VectorType& vec_def, const VectorType& vec_sol) override
        {
          // store new defect
          this->_def_init = this->_def_cur = this->_calc_def_norm(vec_def, vec_sol);
          this->_num_iter = Index(0);
          this->_num_stag_iter = Index(0);

          _fval_init = _fval;
          _fval_prev = Math::huge<DataType>();
          _steplength = DataType(1);
          _ls_its = Index(0);

          Statistics::add_solver_expression(
            std::make_shared<ExpressionDefect>(this->name(), this->_def_init, this->get_num_iter()));

          if(this->_plot)
          {
            std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            <<  " (" << stringify(this->_ls_its).pad_front(_ls_iter_digits) << ")"
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << " : " << stringify_fp_sci(this->_fval)
            << std::endl;
          }

          // Ensure that the initial fval and defect are neither NaN nor infinity
          if( !(Math::isfinite(this->_def_init) && Math::isfinite(this->_fval_init)) )
          {
            return Status::aborted;
          }

          // Check if the initial defect is already low enough
          if(this->_def_init <= Math::sqr(Math::eps<DataType>()))
          {
            return Status::success;
          }

          // continue iterating
          return Status::progress;
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
          {
            this->_def_cur = this->_calc_def_norm(vec_r, vec_sol);
            Statistics::add_solver_expression(
              std::make_shared<ExpressionDefect>(this->name(), this->_def_cur, this->get_num_iter()));
          }

          // plot?
          if(this->_plot)
          {
            std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            <<  " (" << stringify(this->_ls_its).pad_front(_ls_iter_digits) << ")"
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << " : " << stringify_fp_sci(this->_fval)
            << " : " << stringify_fp_sci(this->_steplength)
            << std::endl;
          }

          // ensure that the defect is neither NaN nor infinity
          if(!Math::isfinite(this->_def_cur))
          {
            return Status::aborted;
          }

          // is diverged?
          if(this->is_diverged())
          {
            return Status::diverged;
          }

          // minimum number of iterations performed?
          if(this->_num_iter < this->_min_iter)
          {
            return Status::progress;
          }

          // maximum number of iterations performed?
          if(this->_num_iter >= this->_max_iter)
          {
            return Status::max_iter;
          }

          // Check for convergence of the gradient norm (relative AND absolute criterion)
          if(this->is_converged())
          {
            return Status::success;
          }

          // Check for convergence wrt. the function value improvement if _tol_fval says so
          if(_tol_fval > DataType(0))
          {
            // This is the factor for the relative function value
            DataType scale(Math::max(_fval, _fval_prev));
            // Make sure it is at least 1
            scale = Math::max(scale, DataType(1));
            // Check for success
            if(Math::abs(_fval_prev - _fval) <= _tol_fval*scale)
            {
              return Status::success;
            }
          }

          if(_steplength <= _tol_step)
          {
            return Status::success;
          }

          // If there were too many stagnated iterations, the solver is stagnated
          if(this->_min_stag_iter > 0 && this->_num_stag_iter > this->_min_stag_iter)
          {
            return Status::stagnated;
          }

          // continue iterating
          return Status::progress;
        }

    };
  } // namespace Solver
} //namespace FEAT
#endif // FEAT_SOLVER_NLOPTLS
