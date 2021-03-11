// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ANALYTIC_AUTO_DERIVE_HPP
#define KERNEL_ANALYTIC_AUTO_DERIVE_HPP 1

// includes, FEAT
#include <kernel/analytic/function.hpp>

// includes, system
#include <vector>

namespace FEAT
{
  namespace Analytic
  {
    /// \cond internal
    namespace Intern
    {
      template<bool wrap_>
      struct AutoDeriveGradWrapper
      {
        template<typename Point_, typename FuncEval_, typename AutoDer_>
        static typename AutoDer_::GradientType wrap(const Point_& point, FuncEval_& func_eval, AutoDer_&)
        {
          return func_eval.gradient(point);
        }
      };

      template<>
      struct AutoDeriveGradWrapper<false>
      {
        template<typename Point_, typename FuncEval_, typename AutoDer_>
        static typename AutoDer_::GradientType wrap(const Point_& point, FuncEval_&, AutoDer_& auto_der)
        {
          return auto_der.extrapol_grad(point);
        }
      };

      template<bool wrap_>
      struct AutoDeriveHessWrapper
      {
        template<typename Point_, typename FuncEval_, typename AutoDer_>
        static typename AutoDer_::HessianType wrap(const Point_& point, FuncEval_& func_eval, AutoDer_&)
        {
          return func_eval.hessian(point);
        }
      };

      template<>
      struct AutoDeriveHessWrapper<false>
      {
        template<typename Point_, typename FuncEval_, typename AutoDer_>
        static typename AutoDer_::HessianType wrap(const Point_& point, FuncEval_&, AutoDer_& auto_der)
        {
          return auto_der.extrapol_hess(point);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Auto-Derive function wrapper class template
     *
     * This class extends another given function by adding the capability of computing
     * gradients and hessians via Richardson extrapolation applied onto second-order
     * central difference quotients. The initial 'h' for the difference quotient as well
     * as the maximum number of extrapolation steps can be adjusted by using the
     * #config_grad_extrapol() and #config_hess_extrapol() functions, respectively.
     *
     * \note
     * If the original function already supports the computation of gradients and/or
     * hessians, then the original function's computation is used instead of the
     * extrapolation scheme. Example: If the input function supports computation of
     * function values and gradients but not hessians, then this class template
     * additionally provides the numerical computation of hessians, but uses the
     * original function's evaluator for the computation of gradients.
     *
     * \tparam Function_
     * The function to which the numerical computation of gradients and/or hessians is
     * to be added. This function will be used as a base-class for this template instance.
     *
     * \tparam DataType_
     * The floating point data type to be used to store the values of the initial h.
     *
     * \author Peter Zajac
     */
    template<typename Function_, typename DataType_ = Real>
    class AutoDerive :
      public Function_
    {
    public:
      /// the input function must support value computation
      static_assert(Function_::can_value, "function cannot compute values");

      /// specify our domain dimension
      static constexpr int domain_dim = Function_::domain_dim;

      /// specify our image type
      typedef typename Function_::ImageType ImageType;

      /// our base class provides function values
      static constexpr bool can_value = true;
      /// we provide function gradients
      static constexpr bool can_grad = true;
      /// we provide function hessiants
      static constexpr bool can_hess = true;

      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      public:
        /// coefficient data type
        typedef typename Traits_::DataType DataType;
        /// evaluation point type
        typedef typename Traits_::PointType PointType;
        /// value type
        typedef typename Traits_::ValueType ValueType;
        /// gradient type
        typedef typename Traits_::GradientType GradientType;
        /// hessian type
        typedef typename Traits_::HessianType HessianType;
        /// get our dimension
        static constexpr int domain_dim = Traits_::domain_dim;

      private:
        /// the original function's evaluator
        typename Function_::template Evaluator<Traits_> _func_eval;
        /// our gradient extrapolation table
        std::vector<GradientType> _grad;
        /// our hessian extrapolation table
        std::vector<HessianType> _hess;
        /// initial h for gradient extrapolation
        const DataType _init_grad_h;
        /// initial h for hessian extrapolation
        const DataType _init_hess_h;

      public:
        /// mandatory CTOR
        explicit Evaluator(const AutoDerive& function) :
          _func_eval(function),
          _grad(std::size_t(function._max_grad_steps)),
          _hess(std::size_t(function._max_hess_steps)),
          _init_grad_h(function._init_grad_h),
          _init_hess_h(function._init_hess_h)
        {
        }

        ValueType value(const PointType& point)
        {
          return _func_eval.value(point);
        }

        GradientType gradient(const PointType& point)
        {
          // Depending on whether the original function can compute gradients,
          // either call the original evaluator's "gradient" function or apply
          // our extrapolation scheme implemented in "extrapol_grad".
          return Intern::AutoDeriveGradWrapper<Function_::can_grad>::wrap(point, _func_eval, *this);
        }

        HessianType hessian(const PointType& point)
        {
          // Depending on whether the original function can compute hessians,
          // either call the original evaluator's "hessian" function or apply
          // our extrapolation scheme implemented in "extrapol_hess".
          return Intern::AutoDeriveHessWrapper<Function_::can_hess>::wrap(point, _func_eval, *this);
        }

        /**
         * \brief Computes the gradient by a Richardson extrapolation scheme.
         */
        GradientType extrapol_grad(const PointType& point)
        {
          // first, create a mutable copy of our point
          PointType v(point);

          // next, choose the initial h
          DataType h(_init_grad_h);

          // Note: the '_grad' vector was already allocated to
          //       the correct size by our constructor

          // evaluate gradient
          this->_eval_grad_quot(_grad[0], v, h);

          // Now comes the Richardson extrapolation loop
          const std::size_t n(_grad.size()-1);
          DataType def(Math::huge<DataType>());
          for(std::size_t i(0); i < n; ++i)
          {
            // evaluate next gradient
            this->_eval_grad_quot(_grad[i+1], v, h *= DataType(0.5));

            // initialize scaling fator
            DataType q = DataType(1);

            // perform extrapolation steps except for the last one
            for(std::size_t k(i); k > std::size_t(0); --k)
            {
              q *= DataType(4);
              (_grad[k] -= q*_grad[k+1]) *= DataType(1) / (DataType(1) - q);
            }

            // compute and check our defect
            DataType d = def_norm_sqr(_grad[1] - _grad[0]);
            if(def <= d)
            {
              // The defect has increased, so we return the previous extrapolation result
              return _grad.front();
            }

            // remember current defect
            def = d;

            // perform last extrapolation step
            q *= DataType(4);
            (_grad[0] -= q*_grad[1]) *= DataType(1) / (DataType(1) - q);
          }

          // return our extrapolated gradient
          return _grad.front();
        }

        /**
         * \brief Computes the hessian by a Richardson extrapolation scheme.
         */
        HessianType extrapol_hess(const PointType& point)
        {
          // first, create a mutable copy of our point
          PointType v(point);

          // next, choose the initial h
          DataType h(_init_hess_h);

          // Note: the '_hess' vector was already allocated to
          //       the correct size by our constructor

          // evaluate hessian
          this->_eval_hess_quot(_hess[0], v, h);

          // Now comes the Richardson extrapolation loop
          const std::size_t n(_hess.size()-1);
          DataType def(Math::huge<DataType>());
          for(std::size_t i(0); i < n; ++i)
          {
            // evaluate next hessian
            this->_eval_hess_quot(_hess[i+1], v, h *= DataType(0.5));

            // initialize scaling fator
            DataType q = DataType(1);

            // perform extrapolation steps except for the last one
            for(std::size_t k(i); k > std::size_t(0); --k)
            {
              q *= DataType(4);
              (_hess[k] -= q*_hess[k+1]) *= DataType(1) / (DataType(1) - q);
            }

            // compute and check our defect
            DataType d = def_norm_sqr(_hess[1] - _hess[0]);
            if(def <= d)
            {
              // The defect has increased, so we return the previous extrapolation result
              return _hess.front();
            }

            // remember current defect
            def = d;

            // perform last extrapolation step
            q *= DataType(4);
            (_hess[0] -= q*_hess[1]) *= DataType(1) / (DataType(1) - q);
          }

          // return our extrapolated hessian
          return _hess.front();
        }

      protected:
        template<int n_, int s_>
        static DataType def_norm_sqr(const Tiny::Vector<DataType, n_, s_>& x)
        {
          return x.norm_euclid_sqr();
        }

        template<int m_, int n_, int sm_, int sn_>
        static DataType def_norm_sqr(const Tiny::Matrix<DataType, m_, n_, sm_, sn_>& x)
        {
          return x.norm_hessian_sqr();
        }

        template<int l_, int m_, int n_, int sl_, int sm_, int sn_>
        static DataType def_norm_sqr(const Tiny::Tensor3<DataType, l_, m_, n_, sl_, sm_, sn_>& x)
        {
          DataType r(DataType(0));
          for(int i(0); i < l_; ++i)
          {
            for(int j(0); j < m_; ++j)
            {
              for(int k(0); k < n_; ++k)
              {
                r += Math::sqr(x(i,j,k));
              }
            }
          }
          return r;
        }

        /// evaluates the first-order difference quotients
        void _eval_grad_quot(GradientType& x, PointType& v, const DataType h)
        {
          // difference quotient denominator
          const DataType denom = DataType(1) / (DataType(2) * h);
          ValueType vr, vl;

          // loop over all dimensions
          for(int i(0); i < domain_dim; ++i)
          {
            // backup coordinate i
            const DataType vi(v[i]);

            // evaluate f(v + h*e_i)
            v[i] = vi + h;
            vr = _func_eval.value(v);

            // evaluate f(v - h*e_i)
            v[i] = vi - h;
            vl = _func_eval.value(v);

            // compute difference quotient
            x[i] = denom * (vr - vl);

            // restore coord
            v[i] = vi;
          }
        }

        /// evaluates the second-order difference quotients
        void _eval_hess_quot(HessianType& x, PointType& v, const DataType h)
        {
          // difference quotient denominators
          const DataType denom1 = DataType(1) / (h * h);
          const DataType denom2 = DataType(1) / (DataType(4) * h * h);
          ValueType vc, vr, vl, vne, vnw, vse, vsw;

          // evaluate at point
          vc = _func_eval.value(v);

          // loop over all dimensions
          for(int i(0); i < domain_dim; ++i)
          {
            // backup coord
            const DataType vi(v[i]);

            // eval f(x + h*e_i)
            v[i] = vi + h;
            vr = _func_eval.value(v);

            // eval f(x-h)
            v[i] = vi - h;
            vl = _func_eval.value(v);

            // compute difference quotient
            x[i][i] = denom1 * (vr + vl - DataType(2)*vc);

            // now the mixed derivatives
            for(int j(i+1); j < domain_dim; ++j)
            {
              // backup coord
              const DataType vj(v[j]);

              // we need four points here:
              // north-east: f(v + h*e_i + h*e_j)
              v[i] = vi + h;
              v[j] = vj + h;
              vne = _func_eval.value(v);

              // north-west: f(v - h*e_i + h*e_j)
              v[i] = vi - h;
              v[j] = vj + h;
              vnw = _func_eval.value(v);

              // south-east: f(v + h*e_i - h*e_j)
              v[i] = vi + h;
              v[j] = vj - h;
              vse = _func_eval.value(v);

              // south-west: f(v - h*e_i - h*e_j)
              v[i] = vi - h;
              v[j] = vj - h;
              vsw = _func_eval.value(v);

              // combine into difference quotient
              x[i][j] = x[j][i] = denom2 * ((vne + vsw) - (vnw + vse));

              // restore coord
              v[j] = vj;
            }

            // restore coord
            v[i] = vi;
          }
        }
      }; // class AutoDeriveBase::Evaluator<...>

    protected:
      /// maximum number of gradient extrapolation steps
      int _max_grad_steps;
      /// maximum number of hessian extrapolation steps
      int _max_hess_steps;
      /// initial h for gradient extrapolation
      DataType_ _init_grad_h;
      /// initial h for hessian extrapolation
      DataType_ _init_hess_h;

    public:
      /**
       * \brief Standard constructor
       *
       * This constructor configures the Richardson extrapolation steps
       * for both the gradient and hessian evaluation to:
       * - initial h = 1E-2
       * - maximum steps = 10
       * which should be sufficient in most (non-exotic) cases.
       *
       * \param[in] args
       * The set of arguments to be forwarded to the original function's constructor.
       */
      template<typename... Args_>
      explicit AutoDerive(Args_&&... args) :
        Function_(std::forward<Args_>(args)...),
        _max_grad_steps(10),
        _max_hess_steps(10),
        _init_grad_h(1E-2),
        _init_hess_h(1E-2)
      {
      }

      /**
       * \brief Configures the Gradient extrapolation scheme.
       *
       * This function configures the Richardson extrapolation scheme for the
       * evaluation of function gradients (first-order derivatives).
       *
       * \param[in] initial_h
       * The initial h for the difference quotent. Must be > 0; default = 1E-2.
       *
       * \param[in] max_steps
       * The maximum number of Richardson extrapolation steps to be performed.
       * Must be > 0; default = 10.
       */
      void config_grad_extrapol(DataType_ initial_h, int max_steps)
      {
        XASSERTM(initial_h > Math::eps<DataType_>(), "Initial h is too small or non-positive!");
        XASSERTM(max_steps > 0, "Invalid maximum extrapolation steps!");
        _init_grad_h = initial_h;
        _max_grad_steps = max_steps;
      }

      /**
       * \brief Configures the Hessian extrapolation scheme.
       *
       * This function configures the Richardson extrapolation scheme for the
       * evaluation of function hessians (second-order derivatives).
       *
       * \param[in] initial_h
       * The initial h for the difference quotent. Must be > 0; default = 1E-2.
       *
       * \param[in] max_steps
       * The maximum number of Richardson extrapolation steps to be performed.
       * Must be > 0; default = 10.
       */
      void config_hess_extrapol(DataType_ initial_h, int max_steps)
      {
        XASSERTM(initial_h > Math::eps<DataType_>(), "Initial h is too small or non-positive!");
        XASSERTM(max_steps > 0, "Invalid maximum extrapolation steps!");
         _init_hess_h = initial_h;
        _max_hess_steps = max_steps;
      }
    }; // class AutoDerive<...>
  } // namespace Analytic
} // namespace FEAT

#endif // KERNEL_ANALYTIC_AUTO_DERIVE_HPP
