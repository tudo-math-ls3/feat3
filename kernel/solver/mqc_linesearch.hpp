// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef FEAT_KERNEL_SOLVER_MQC_LINESEARCH
#define FEAT_KERNEL_SOLVER_MQC_LINESEARCH 1
#include <kernel/base_header.hpp>
#include <kernel/solver/linesearch.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Mixed quadratic-cubic line search
     *
     * \tparam Functional_
     * The (nonlinear) functional to be evaluated
     *
     * \tparam Filter_
     * The filter to be applied to the functional's gradient
     *
     * This class implements a linesearch which approximately finds
     * \f[
     *   \alpha^* = \mathrm{argmin} \nabla f(x + \alpha d) \cdot d
     * \f]
     *
     * This is quite complicated code in some places, and very sensitive to annihilation as lots of information is
     * condensed into very few floating point values (i.e. the dot products of the gradient and the search direction
     * for a highly nonlinear and high dimensional problem).
     *
     * The line search works with two potentially overlapping intervals: The interval of uncertainty
     * \f$ [\alpha_{\mathrm{soft~min}}, \alpha_{\mathrm{soft~max}}] \f$ and the current search interval
     * \f$ [\alpha_{\mathrm{lo}}, \alpha_{\mathrm{hi}}] \f$.
     *
     * Intitially, we have an initial step length \f$ \alpha_0 \f$ which defaults to 1 if it was not set explicitly
     * before and set \f$ \alpha_{\mathrm{soft~min}} = 0, \alpha_{\mathrm{soft~max}} = \alpha_0 + c(\alpha_0 -
     * \alpha_{\mathrm{soft~min}}) \f$, where the constant is \f$ c = 4\f$ which is pretty arbitrary.
     * The search interval is the degenerate interval \f$ [\alpha_0, \alpha_0]\f$. Not that we know that the
     * derivative \f$ \left< d, grad(\alpha_0) \right> < 0 \f$ because \f$ d \f$ is a descent direction.
     *
     * At this point we do not know if the minimum we are looking for is in either interval (it is not, in general),
     * so we first need to expand the interval of uncertainty until we are sure that it contains the minimum.
     * This is the case if i.e. the derivative changes sign for some step size \f$ \alpha > \alpha_0 \f$, or if the
     * function value decreases so the sufficient decrease condition is satisfied along with the curvature condition.
     *
     *
     */
    template<typename Functional_, typename Filter_>
    class MQCLinesearch : public Linesearch<Functional_, Filter_>
    {
      public:
        /// Filter type to be applied to the gradient of the functional
        typedef Filter_ FilterType;
        /// Input vector type for the functional's gradient
        typedef typename Functional_::VectorTypeR VectorType;
        /// Underlying floating point type
        typedef typename Functional_::DataType DataType;
        /// Our base class
        typedef Linesearch<Functional_, Filter_> BaseClass;

      protected:
        /// Hard maximum for the step length
        DataType _alpha_hard_max;
        /// Hard minimum for the step length
        DataType _alpha_hard_min;
        /// Lower bound of the interval of uncertainty
        DataType _alpha_soft_max;
        /// Upper bound of the interval of uncertainty
        DataType _alpha_soft_min;

      public:
        /**
         * \brief Standard constructor
         *
         * \param[in, out] functional
         * The (nonlinear) functional. Cannot be const because it saves its own state
         *
         * \param[in] filter
         * Filter to apply to the functional's gradient
         *
         * \param[in] keep_iterates
         * Keep all iterates in a std::deque. Defaults to false.
         *
         */
        explicit MQCLinesearch(Functional_& functional, Filter_& filter, bool keep_iterates = false) :
          BaseClass("MQC-LS", functional, filter, keep_iterates),
          _alpha_hard_max(DataType(0)),
          _alpha_hard_min(DataType(0)),
          _alpha_soft_max(DataType(0)),
          _alpha_soft_min(DataType(0))
          {
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
         * \param[in] functional
         * The functional.
         *
         * \param[in] filter
         * The system filter.
         *
         */
        explicit MQCLinesearch(const String& section_name, PropertyMap* section,
        Functional_& functional, Filter_& filter) :
          BaseClass("MQC-LS", section_name, section, functional, filter),
          _alpha_hard_max(DataType(0)),
          _alpha_hard_min(DataType(0)),
          _alpha_soft_max(DataType(0)),
          _alpha_soft_min(DataType(0))
          {
          }

        /// \copydoc ~BaseClass()
        virtual ~MQCLinesearch()
        {
        }

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "MQCLinesearch";
        }

        /// \copydoc BaseClass::reset()
        virtual void reset() override
        {
          BaseClass::reset();
          _alpha_hard_max = DataType(0);
          _alpha_hard_min = DataType(0);
          _alpha_soft_max = DataType(0);
          _alpha_soft_min = DataType(0);
        }

        /**
         * \brief Applies the solver, setting the initial guess to zero.
         *
         * \param[out] vec_cor
         * Solution, gets zeroed
         *
         * \param[in] vec_dir
         * Search direction
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        virtual Status apply(VectorType& vec_cor, const VectorType& vec_dir) override
        {
          // clear solution vector
          vec_cor.format();
          this->_functional.prepare(vec_cor, this->_filter);

          // apply
          this->_status = _apply_intern(vec_cor, vec_dir);
          this->plot_summary();
          return this->_status;
        }

        /**
         * \brief Applies the solver, making use of an initial guess
         *
         * \param[out] vec_sol
         * Initial guess, gets overwritten by the solution
         *
         * \param[in] vec_dir
         * Search direction
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        virtual Status correct(VectorType& vec_sol, const VectorType& vec_dir) override
        {
          this->_functional.prepare(vec_sol, this->_filter);

          // apply
          this->_status = _apply_intern(vec_sol, vec_dir);
          this->plot_summary();
          return this->_status;
        }

      protected:
        /**
         * \brief Internal function: Applies the solver
         *
         * \param[in, out] vec_sol
         * Initial guess, gets overwritten by solution
         *
         * \param[in] vec_dir
         * Search direction
         *
         * \note This assumes that the initial functional value _fval_0 and the gradient were already set from the
         * calling solver!
         *
         * \returns The solver status (success, max_iter, stagnated)
         */
        virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_dir)
        {
          Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

          static constexpr DataType extrapolation_width = DataType(4);
          Status status(Status::progress);

          // The step length wrt. to the NORMALISED search direction
          DataType alpha(0);
          // The functional value
          DataType fval(0);
          /// <vec_pn, vec_grad>. We want to find the minimum of the functional value along vec_pn
          DataType df(0);

          // Perform initialisations and checks
          Status st = this->_startup(alpha, fval, df, vec_sol, vec_dir);
          //alpha = this->_alpha_0;
          if(st != Status::progress)
          {
            return st;
          }

          // Scaling if we are to use step sizes wrt. to the non-normalised search direction
          // This appears to be the right thing theoretically, but in practice not using the search direction norm as
          // initial guess for the step size works better. Stupid reality...
          //if(this->_dir_scaling)
          //{
          //  alpha = this->_norm_dir;
          //}

          // Set hard limits to default values if they have not been set
          _alpha_hard_min = DataType(0);
          // This bogus value should be around 1e50 for double precision. It is chosen to make comparing results with
          // ALGLIB easier
          if(_alpha_hard_max < Math::eps<DataType>())
          {
            _alpha_hard_max = Math::pow(Math::huge<DataType>(), DataType(0.1622));
          }

          // Set the intervall of uncertainty
          _alpha_soft_min = DataType(0);
          _alpha_soft_max = this->_alpha_0 + extrapolation_width*this->_alpha_0;
          _alpha_soft_max = Math::min(_alpha_soft_max, _alpha_hard_max);

          //std::cout << "Linesearch initial alpha " << this->_alpha_0 << std::endl;

          DataType alpha_lo(0);
          // It is critical that _f_0 was set from the outside!
          DataType fval_lo(this->_fval_0);
          DataType df_lo(this->_delta_0);

          DataType alpha_hi(0);
          DataType fval_hi(this->_fval_0);
          DataType df_hi(this->_delta_0);

          // This is the width of the search interval
          DataType width(Math::abs(_alpha_hard_max - _alpha_hard_min));
          DataType width_old(DataType(2)*width);

          // Do we know the interval of uncertainty?
          bool interval_known(false);
          // Is the minimum in the interval of uncertainty?
          bool min_in_interval(false);
          // Does the minimum lie outside the current search interval, requiring us to drive the step to its
          // boundary?
          bool drive_to_bndry(false);

          // start iterating
          while(status == Status::progress)
          {

            IterationStats stat(*this);
            ++(this->_num_iter);

            // If we know the minimum is in the search interval, the interval of uncertainty is the search interval
            if(min_in_interval)
            {
              _alpha_soft_min = Math::min(alpha_lo, alpha_hi);
              _alpha_soft_max = Math::max(alpha_lo, alpha_hi);
            }
            // Enlarge the interval of uncertainty
            else
            {
              _alpha_soft_min = alpha_lo;
              _alpha_soft_max = alpha + extrapolation_width*Math::abs(alpha - alpha_lo);
            }
            //std::cout << "Set alpha_smin " << _alpha_soft_min << " alpha_smax " << _alpha_soft_max << std::endl;

            // Update solution: sol <- initial_sol + _alpha*dir
            vec_sol.axpy(this->_vec_pn, this->_vec_initial_sol, alpha);
            //std::cout << "Linesearch alpha " << alpha << std::endl;
            //std::cout << "initial_sol " << *(this->_vec_initial_sol) << std::endl;
            //std::cout << "dir " << *this->_vec_pn << std::endl;
            //std::cout << "sol " << *vec_sol << std::endl;

            // Prepare and evaluate
            this->_functional.prepare(vec_sol, this->_filter);

            // Compute and filter the gradient
            this->_functional.eval_fval_grad(fval, this->_vec_grad);
            this->trim_func_grad(fval);
            this->_filter.filter_def(this->_vec_grad);

            // New directional derivative and new defect
            df = this->_vec_pn.dot(this->_vec_grad);
            status = this->_check_convergence(fval, df, alpha);

            // If success is reported, check if it is a real success or if something fishy is going on
            if(status == Status::success)
            {
              this->_vec_tmp.axpy(vec_sol, this->_vec_initial_sol, -DataType(1));
              if(fval >= this->_fval_0 || this->_vec_tmp.norm2() == DataType(0))
              {
                status = Status::stagnated;
              }
            }
            // If we are not successful, check if the interval of uncertainty has become too small
            else if(min_in_interval && (_alpha_soft_max - _alpha_soft_min) <= this->_tol_step*_alpha_soft_max)
            {
              //std::cout << "interval width " << _alpha_soft_max - _alpha_soft_min << " : " << this->_tol_step*_alpha_soft_max << std::endl;
              status = Status::stagnated;
            }
            // Stagnation due to rounding errors
            //if(min_in_interval && (alpha <= _alpha_soft_min || alpha >= _alpha_soft_max))
            //{
            //  std::cout << "Rounding errors" << std::endl;
            //  status = Status::stagnated;
            //}
            // This is not used at the moment because it is only relevant if there are constraints limiting the
            // step length
            //if( alpha == _alpha_hard_max)
            //if( alpha == _alpha_hard_min)

            if(status != Status::progress)
            {
              break;
            }

            if(!interval_known
                && (fval < this->_fval_0 + this->_tol_decrease*alpha*this->_delta_0)
                && (df >= Math::min(this->_tol_decrease, this->_tol_curvature)))
                {
                  interval_known = true;
                }

            //std::cout << "interval known " << interval_known << std::endl;
            // If we do not know that the minimum was in the previous interval of uncertainty, we need to compute a
            // new step size to expand the interval of uncertainty at the start of the next iteration.
            if(!interval_known && (fval <= fval_lo)
                && fval > this->_fval_0 + alpha*this->_tol_decrease*this->_delta_0)
            {
              DataType fval_m(fval - alpha*this->_tol_decrease*this->_delta_0);
              DataType df_m(df - this->_tol_decrease*this->_delta_0);
              DataType fval_lo_m(fval_lo - alpha_lo*this->_tol_decrease*this->_delta_0);
              DataType df_lo_m(df_lo - this->_tol_decrease*this->_delta_0);
              DataType fval_hi_m(fval_hi - alpha_hi*this->_tol_decrease*this->_delta_0);
              DataType df_hi_m(df_hi - this->_tol_decrease*this->_delta_0);

              // Note that the expansion step might already give us the information that the minimum is the the
              // new search interval
              _polynomial_fit(
                alpha, fval_m, df_m, alpha_lo, fval_lo_m, df_lo_m,
                alpha_hi, fval_hi_m, df_hi_m, min_in_interval, drive_to_bndry);

              fval_lo = fval_lo_m + alpha_lo*this->_tol_decrease*this->_delta_0;
              df_lo = df_lo_m + this->_tol_decrease*this->_delta_0;
              fval_hi = fval_hi_m * alpha_hi*this->_tol_decrease*this->_delta_0;
              df_hi = df_hi_m + this->_tol_decrease*this->_delta_0;

            }
            else
            {
              _polynomial_fit(alpha, fval, df, alpha_lo, fval_lo, df_lo, alpha_hi, fval_hi, df_hi, min_in_interval,
              drive_to_bndry);

              //std::cout << "width " << width << " width_old " << width_old << " min_in_interval " << min_in_interval <<std::endl;

              if(min_in_interval)
              {
                if(Math::abs(alpha_hi - alpha_lo) >= DataType(0.66)*width_old)
                {
                  //std::cout << "Forcing " << alpha << " to ";
                  alpha = alpha_lo + DataType(0.5)*(alpha_hi - alpha_lo);
                  //std::cout << alpha << std::endl;
                }
                width_old = width;
                width = Math::abs(alpha_lo - alpha_hi);
              }
            }

          } // while(status == Status:progress)

          this->_alpha_min = alpha;
          this->_fval_min = fval;

          // If we are successful, we save the last step length as the new initial step length
          if(status == Status::success)
          {
            this->_alpha_0 = this->_alpha_min;
          }
          // If we are not successful, we update the best step length and need to re-evaluate everything for that
          // step
          else
          {
            this->_alpha_min = alpha_lo;
            this->_fval_min = fval_lo;
            //std::cout << "Unusual termination alpha_lo " << alpha_lo << "fval_lo " << fval_lo << std::endl;
            vec_sol.axpy(this->_vec_pn, this->_vec_initial_sol, alpha_lo);

            // Prepare and evaluate
            this->_functional.prepare(vec_sol, this->_filter);
            this->_functional.eval_fval_grad(fval, this->_vec_grad);
            this->trim_func_grad(fval);
            this->_filter.filter_def(this->_vec_grad);
          }

          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
          return status;
        }

        /**
         * \brief The great magick trick to find a minimum of a 1d function
         *
         * \param[out] alpha
         * The trial step.
         *
         * \param[in] fval
         * Functional value at alpha
         *
         * \param[in] df
         * Derivative value at alpha
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] fval_lo
         * Functional value at alpha_lo
         *
         * \param[in] df_lo
         * Derivative at alpha_lo
         *
         * \param[in] alpha_hi
         * The step size for the lower functional value
         *
         * \param[in] fval_hi
         * Functional value at alpha_hi
         *
         * \param[in] df_hi
         * Derivative at alpha_hi
         *
         * \param[in,out] min_in_interval
         * Do we know the minimum has to be in the interval (or any superset we worked on before)? Because this gets
         * called several times, this knowledge has to be passed back and forth.
         *
         * \param[in,out] drive_to_bndry
         * In the cases 1 and 3 this is set to true, in cases 2 and 4 to false. This determines if convergence can be
         * improved by driving the step size to the appropriate boundary more quickly.
         *
         * \returns
         * Status::success as there is nothing that can fail.
         *
         * This routine does a quadratic interpolation using fval_lo, df_lo and df, and a cubic interpolation using
         * fval_lo, df_lo, fval_hi and df_hi. The great magick is to determine which is better, and this depends on a
         * lot of factors, which leads to 4 distinct cases.
         *
         */
        Status _polynomial_fit(
          DataType& alpha, DataType& fval, DataType& df,
          DataType& alpha_lo, DataType& fval_lo, DataType& df_lo,
          DataType& alpha_hi, DataType& fval_hi, DataType& df_hi,
          bool& min_in_interval, bool& drive_to_bndry)
        {
          DataType alpha_c(0);
          DataType alpha_q(0);
          DataType alpha_new(0);

          // Case 1: The function value increased.
          if( fval > fval_lo)
          {
            // Computation of cubic step
            alpha_c  = _argmin_cubic(alpha_lo, fval_lo, df_lo, alpha, fval, df);
            // The first derivative was negative, but the function value increases, so the minimum has to be in this
            // interval (or have been in the original interval, which this interval is a subset of)
            min_in_interval = true;
            drive_to_bndry = true;
            // Use 2nd order approximation NOT using the derivative df because it might be garbage (i.e. if we hit a
            // singularity) because the function value increased
            alpha_q = _argmin_quadratic(alpha_lo, fval_lo, df_lo, alpha, fval, df, false);

            // If the cubic minimum is closer to the old best step length, take it. Otherwise, take the mean
            if(Math::abs(alpha_c - alpha_lo) < Math::abs(alpha_q - alpha_lo))
              alpha_new = alpha_c;
            else
              alpha_new = (alpha_q + alpha_c)/DataType(2);
            //std::cout << "Case 1: alpha " << alpha_new << " q " << alpha_q << " c " << alpha_c << std::endl;

          }
          // Case 2: The first derivative changes sign
          // We also know that the function value increased in the last trial step
          else if( df_lo*df < DataType(0) )
          {
            alpha_c  = _argmin_cubic(alpha, fval, df, alpha_lo, fval_lo, df_lo);
            // Because the derivative changed sign, it has to be zero somewhere in between, so the minimum is in the
            // current interval
            min_in_interval = true;
            drive_to_bndry = false;
            // Default: Quadratic interpolant using fval_lo, df_lo and df
            alpha_q = _argmin_quadratic(alpha_lo, fval_lo, df_lo, alpha, fval, df, true);
            // Take the step closer to the new step alpha
            Math::abs(alpha - alpha_c) > Math::abs(alpha - alpha_q) ? alpha_new = alpha_c : alpha_new = alpha_q;
            //std::cout << "Case 2: alpha " << alpha_new << " q " << alpha_q << " c " << alpha_c << std::endl;
          }
          // Case 3: The absolute value of the derivative increases
          // We also know that the function value increased in the last trial step and that the derivative did not
          // change sign.
          // This means we have to further search in the direction of alpha.
          else if(Math::abs(df) < Math::abs(df_lo))
          {
            alpha_c  = _argmin_cubic_case_3(alpha, fval, df, alpha_lo, fval_lo, df_lo);
            drive_to_bndry = true;

            // Quadratic interpolant using fval_lo, df_lo and df
            alpha_q = _argmin_quadratic(alpha, fval, df, alpha_lo, fval_lo, df_lo, true);

            // If the cubic step is closer to the trial step
            if(Math::abs(alpha - alpha_c) < Math::abs(alpha - alpha_q))
            {
              min_in_interval ? alpha_new = alpha_c : alpha_new = alpha_q;
            }
            else
            {
              min_in_interval ? alpha_new = alpha_q : alpha_new = alpha_c;
            }
            //std::cout << "Case 3: alpha " << alpha_new << " q " << alpha_q << " c " << alpha_c << std::endl;
          }
          // Case 4: The absolute value of the derivative did not increase
          // We also know that the function value increased in the last trial step and that the derivative did not
          // change sign.
          else
          {
            drive_to_bndry = false;
            // Note that the arguments for argmin_cubic are different in this case
            // Here we can just take the cubic step because the cubic interpolation sets it to the alpha_soft_max
            // if it recognises that the minimum is outside the current interval.
            if(min_in_interval)
            {
              alpha_c  = _argmin_cubic(alpha, fval , df, alpha_hi, fval_hi, df_hi);
              alpha_new = alpha_c;
            }
            else
            {
              alpha > alpha_lo ? alpha_new = _alpha_soft_max : alpha_new = _alpha_soft_min;
            }
            //std::cout << "Case 4: alpha " << alpha_new << std::endl;
          }

            // Update the inverval of uncertainty. Has to happen befor we clamp the step to the admissible interval.
            // Check which step(s) we need to replace
            if(fval > fval_lo)
            {
              alpha_hi = alpha;
              fval_hi = fval;
              df_hi = df;
            }
            else
            {
              if(df * df_lo < DataType(0))
              {
                alpha_hi = alpha_lo;
                fval_hi = fval_lo;
                df_hi = df_lo;
              }
              alpha_lo = alpha;
              fval_lo = fval;
              df_lo = df;
            }

            _clamp_step(alpha_new);

            // If we know the minimum is in the interval, we can drive alpha to the corresponding end more quickly
            if(min_in_interval && drive_to_bndry)
            {
              //std::cout << "Bracketing adjust " << alpha_lo << " " << alpha_hi << " " << alpha_new << std::endl;
              //std::cout << "Bracketing adjust alpha_new from " << alpha_new;
              if(alpha_lo < alpha_hi)
                alpha_new = Math::min(alpha_lo + DataType(0.66)*(alpha_hi - alpha_lo), alpha_new);
              else
                alpha_new = Math::max(alpha_lo + DataType(0.66)*(alpha_hi - alpha_lo), alpha_new);
              //std::cout << " to " << alpha_new << std::endl;
            }

            alpha = alpha_new;

            //std::cout << "Polynomial fit: new " << alpha << " " << fval << " " << df <<std::endl;
            //std::cout << "Polynomial fit: lo  " << alpha_lo << " " << fval_lo << " " << df_lo <<std::endl;
            //std::cout << "Polynomial fit: hi  " << alpha_hi << " " << fval_hi << " " << df_hi <<std::endl;

          return Status::success;
        }

        /**
         * \brief Enforces hard and soft step limits, adjusting the soft limits if necessary
         *
         * \param[in,out] alpha_new
         * The step size to be clamped.
         *
         */
        void _clamp_step(DataType& alpha_new) const
        {
          alpha_new = Math::max(alpha_new, _alpha_hard_min);
          alpha_new = Math::min(alpha_new, _alpha_hard_max);

          alpha_new = Math::max(alpha_new, _alpha_soft_min);
          alpha_new = Math::min(alpha_new, _alpha_soft_max);
        }

        /**
         * \brief Computes the minimum of a quadratic interpolation polynomial
         *
         * This has two variants:
         *  a) Interpolate fval_lo, df_lo and fval_hi
         *  b) Interpolate fval_lo, df_lo and df_hi
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] alpha_hi
         * The step size for the higher functional value
         *
         * \param[in] fval_lo
         * Functional value at alpha_lo
         *
         * \param[in] fval_hi
         * Functional value at alpha_hi
         *
         * \param[in] df_lo
         * Derivative at alpha_lo
         *
         * \param[in] df_hi
         * Derivative at alpha_hi
         *
         * \param[in] interpolate_derivative
         * If this is true, variant a) is chosen
         *
         * \returns
         * The step corresponding to the minimum
         *
         */
        DataType _argmin_quadratic(const DataType alpha_lo, const DataType fval_lo, const DataType df_lo, const DataType alpha_hi, const DataType fval_hi, const DataType df_hi, const bool interpolate_derivative) const
        {
          DataType alpha(alpha_lo);

          // Quadratic interpolation using fval_lo, df_lo, df_hi
          if(interpolate_derivative)
            alpha += df_lo /(df_lo - df_hi) * (alpha_hi - alpha_lo);
          // Quadratic interpolation using fval_lo, fval_hi, df_lo
          else
            alpha += df_lo/( (fval_hi - fval_lo)/( alpha_hi - alpha_lo) - df_lo)/DataType(2)*(alpha_lo - alpha_hi);

          return alpha;
        }

        /**
         * \brief Computes the minimum of a cubic interpolation polynomial
         *
         * If the minimum of the cubic interpolation polynomial is in the interiour, it is computed here. If it is
         * not in the interiour, it lies across one of the endpoints depending on the functional and derivative
         * values. In this case, we set the return value to _alpha_soft_{min, max} accordingly.
         *
         * \param[in] alpha_lo
         * The step size for the lower functional value
         *
         * \param[in] alpha_hi
         * The step size for the higher functional value
         *
         * \param[in] fval_lo
         * Functional value at alpha_lo
         *
         * \param[in] fval_hi
         * Functional value at alpha_hi
         *
         * \param[in] df_lo
         * Derivative at alpha_lo
         *
         * \param[in] df_hi
         * Derivative at alpha_hi
         *
         * \returns
         * The step corresponding to the minimum or _alpha_soft_{min, max} if there is no minimum in the interiour.
         *
         * \note The code is extremely ugly and sensitive to annihilation.
         *
         */
        DataType _argmin_cubic(DataType alpha_lo, DataType fval_lo, DataType df_lo,
        DataType alpha_hi, DataType fval_hi, DataType df_hi) const
        {
          DataType alpha(alpha_lo);

          DataType d1 = DataType(3)*(fval_lo - fval_hi)/(alpha_hi - alpha_lo) + df_lo + df_hi;
          // Scale the computation of r for better numerical stability
          DataType scale = Math::max(Math::abs(df_lo), Math::abs(df_hi));
          scale = Math::max(scale, Math::abs(d1));

          DataType r = Math::sqr(d1/scale) - (df_lo/scale) * (df_hi/scale);

          //std::cout << "fval_hi " << fval_hi << " fval_lo "<< fval_lo << std::endl;
          //std::cout << "alpha_hi " << alpha_hi << " alpha_lo "<< alpha_lo << std::endl;
          //std::cout << "df_hi " << df_hi << " df_lo "<< df_lo << std::endl;
          //std::cout << "d1 " << d1 << std::endl;
          //std::cout << "scale " << scale << std::endl;

          DataType d2(0);

          if(r > DataType(0))
          {
            d2 = Math::signum(alpha_hi - alpha_lo) * scale * Math::sqrt(r);

            DataType p(d2 - df_lo + d1);

            DataType q(0);
            // This sorting is done to avoid annihilation
            df_hi*df_lo > DataType(0) ?  q = d2 +(df_hi-df_lo) + d2 : q = d2 - df_lo + d2 + df_hi;

            //std::cout << "d2 " << d2 << std::endl;
            //std::cout << "p " << p << " q " << q << std::endl;

            alpha += (alpha_hi - alpha_lo)*(p/q);
          }
          // r <= 0 means that the minimum is not in the interior, so it has to lie across the endpoint with the
          // lower function value
          else
          {
            //d2 = Math::signum(alpha_hi - alpha_lo) * scale *
            //  Math::sqrt(Math::max(r, DataType(0)));
            //DataType q(0);
            //if(df_hi*df_lo > 0)
            //  q = d2 +(df_hi-df_lo) + d2;
            //else
            //  q = d2 - df_lo + d2 +df_hi;
            //std::cout << "d2 " << d2 << std::endl;
            //std::cout << "p " << (d2 - df_lo + d1) << " q " << q << std::endl;

            (alpha_lo < alpha_hi && df_lo > DataType(0) )? alpha = _alpha_soft_min: alpha = _alpha_soft_max;
          }

          return alpha;
        }

        /**
         * \copydoc _argmin_cubic
         * This is the weird check if the function value does not tend to infinity in the direction of the cubic
         * step.
         */
        DataType _argmin_cubic_case_3(
          const DataType alpha_lo, const DataType fval_lo, const DataType df_lo,
          const DataType alpha_hi, const DataType fval_hi, DataType df_hi) const
        {
          DataType alpha_c(alpha_lo);

          if(alpha_lo == alpha_hi)
          {
            return alpha_lo;
          }

          DataType d1 = DataType(3)*(fval_hi - fval_lo)/(alpha_lo - alpha_hi) + df_hi + df_lo;
          DataType scale = Math::max(Math::abs(df_lo), Math::abs(df_hi));
          scale = Math::max(scale, Math::abs(d1));

          DataType r = Math::sqr(d1/scale) - (df_lo/scale) * (df_hi/scale);

          DataType d2 = Math::signum(alpha_hi - alpha_lo) * scale * Math::sqrt(Math::max(r, DataType(0)));

          DataType p(d2 - df_lo + d1);
          DataType q(d2 +(df_hi - df_lo) + d2);

          //std::cout << "fval_hi " << fval_hi << " fval_lo "<< fval_lo << std::endl;
          //std::cout << "alpha_hi " << alpha_hi << " alpha_lo "<< alpha_lo << std::endl;
          //std::cout << "df_hi " << df_hi << " df_lo "<< df_lo << std::endl;
          //std::cout << "d1 " << d1 << std::endl;
          //std::cout << "scale " << scale << std::endl;
          //std::cout << "d2 " << d2 << std::endl;
          //std::cout << "p " << (d2 - df_lo + d1) << " q " << q << std::endl;

          if( p/q < DataType(0) && d2 != DataType(0))
            alpha_c += p/q*(alpha_hi - alpha_lo);
          else
          {
            //std::cout << "Weird ";
            alpha_lo > alpha_hi ? alpha_c = _alpha_soft_max : alpha_c = _alpha_soft_min;
          }

          return alpha_c;
        }

    }; // class MQCLinesearch

    /**
     * \brief Creates a new MQCLinesearch object
     *
     * \param[in] functional
     * The functional.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] keep_iterates
     * Flag for keeping the iterates, defaults to false.
     *
     * \returns
     * A shared pointer to a new MQCLinesearch object.
     */
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<MQCLinesearch<Functional_, Filter_>> new_mqc_linesearch(
      Functional_& functional, Filter_& filter, bool keep_iterates = false)
      {
        return std::make_shared<MQCLinesearch<Functional_, Filter_>>(functional, filter, keep_iterates);
      }

    /**
     * \brief Creates a new MQCLinesearch object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself.
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver.
     *
     * \param[in] functional
     * The functional.
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new MQCLinesearch object.
     */
    template<typename Functional_, typename Filter_>
    inline std::shared_ptr<MQCLinesearch<Functional_, Filter_>> new_mqc_linesearch(
      const String& section_name, PropertyMap* section, Functional_& functional, Filter_& filter)
      {
        return std::make_shared<MQCLinesearch<Functional_, Filter_>> (section_name, section, functional, filter);
      }

  } // namespace Solver
} // namespace FEAT
#endif // FEAT_KERNEL_SOLVER_MQC_LINESEARCH
