#pragma once
#ifndef KERNEL_ASSEMBLY_COMMON_FUNCTIONS_HPP
#define KERNEL_ASSEMBLY_COMMON_FUNCTIONS_HPP 1

// includes, FEAST
#include <kernel/assembly/analytic_function.hpp>
#include <kernel/util/math.hpp>

namespace FEAST
{
  namespace Assembly
  {
    namespace Common
    {
      /**
       * \brief Sine-Bubble Static function
       *
       * This class implements the StaticFunction interface representing the function
       *   - 1D: u(x) = sin(pi*x)
       *   - 2D: u(x) = sin(pi*x) * sin(pi*y)
       *   - 3D: u(x) = sin(pi*x) * sin(pi*y) * sin(pi*z)
       *
       * \author Peter Zajac
       */
      template<typename DataType_>
      class SineBubbleStatic
      {
      public:
        /// returns the constant pi
        static DataType_ pi()
        {
          return Math::Limits<DataType_>::pi();
        }

        /// 1D: function value
        static DataType_ eval(DataType_ x)
        {
          return Math::sin(pi() * x);
        }

        /// 2D: function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return Math::sin(pi() * x) * Math::sin(pi() * y);
        }

        /// 3D: function value
        static DataType_ eval(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sin(pi() * x) * Math::sin(pi() * y) * Math::sin(pi() * z);
        }

        /// 1D: X-derivative
        static DataType_ der_x(DataType_ x)
        {
          return pi() * Math::cos(pi() * x);
        }

        /// 2D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return pi() * Math::cos(pi() * x) * Math::sin(pi() * y);
        }

        /// 2D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return pi() * Math::sin(pi() * x) * Math::cos(pi() * y);
        }

        /// 3D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z)
        {
          return pi() * Math::cos(pi() * x) * Math::sin(pi() * y) * Math::sin(pi() * z);
        }

        /// 3D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z)
        {
          return pi() * Math::sin(pi() * x) * Math::cos(pi() * y) * Math::sin(pi() * z);
        }

        /// 3D: Z-derivative
        static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z)
        {
          return pi() * Math::sin(pi() * x) * Math::sin(pi() * y) * Math::cos(pi() * z);
        }

        /// 1D: XX-derivative
        static DataType_ der_xx(DataType_ x)
        {
          return -Math::sqr(pi()) * Math::sin(pi() * x);
        }

        /// 2D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return -Math::sqr(pi()) * Math::sin(pi() * x) * Math::sin(pi() * y);
        }

        /// 2D: YY-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y)
        {
          return der_xx(x, y);
        }

        /// 2D: XY-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y)
        {
          return Math::sqr(pi()) * Math::cos(pi() * x) * Math::cos(pi() * y);
        }

        /// 2D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x, y);
        }

        /// 3D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y, DataType_ z)
        {
          return -Math::sqr(pi()) * Math::sin(pi() * x) * Math::sin(pi() * y) * Math::sin(pi() * z);
        }

        /// 3D: YY-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xx(x, y, z);
        }

        /// 3D: ZZ-derivative
        static DataType_ der_zz(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xx(x, y, z);
        }

        /// 3D: XY-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(pi()) * Math::cos(pi() * x) * Math::cos(pi() * y) * Math::sin(pi() * z);
        }

        /// 3D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xy(x, y, z);
        }

        /// 3D: XZ-derivative
        static DataType_ der_xz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(pi()) * Math::cos(pi() * x) * Math::sin(pi() * y) * Math::cos(pi() * z);
        }

        /// 3D: ZX-derivative
        static DataType_ der_zx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xz(x, y, z);
        }

        /// 3D: YZ-derivative
        static DataType_ der_yz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(pi()) * Math::sin(pi() * x) * Math::cos(pi() * y) * Math::cos(pi() * z);
        }

        /// 3D: ZY-derivative
        static DataType_ der_zy(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_yz(x, y, z);
        }
      }; // class SineBubbleStatic<...>

      /**
       * \brief Sine-Bubble Analytic function
       *
       * This class implements the AnalyticFunction interface representing the function
       *   - 1D: u(x) = sin(pi*x)
       *   - 2D: u(x) = sin(pi*x) * sin(pi*y)
       *   - 3D: u(x) = sin(pi*x) * sin(pi*y) * sin(pi*z)
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \author Peter Zajac
       */
      typedef StaticWrapperFunction<SineBubbleStatic, true, true, true> SineBubbleFunction;

      /**
       * \brief Constant Analytic function
       *
       * This class implements the AnalyticFunction interface representing a constant function.
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \author Peter Zajac
       */
      class ConstantFunction :
        public AnalyticFunction
      {
      public:
        /** \copydoc AnalyticFunction::FunctionCapabilites */
        enum FunctionCapabilities
        {
          can_value = 1,
          can_grad = 1,
          can_hess = 1
        };

        /** \copydoc AnalyticFunction::ConfigTraits */
        template<typename Config_>
        struct ConfigTraits
        {
          typedef Trafo::ConfigBase TrafoConfig;
        };

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public AnalyticFunction::Evaluator<EvalTraits_>
        {
        public:
          /// trafo evaluator data
          typedef typename EvalTraits_::TrafoEvaluator TrafoEvaluator;
          /// trafo data type
          typedef typename EvalTraits_::TrafoData TrafoData;
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          const ConstantFunction& _function;

        public:
          explicit Evaluator(const ConstantFunction& function) :
            _function(function)
          {
          }

          ValueType value(const TrafoData&) const
          {
            return ValueType(_function._value);
          }

          GradientType gradient(const TrafoData&) const
          {
            return GradientType(DataType(0));
          }

          HessianType hessian(const TrafoData&) const
          {
            return HessianType(DataType(0));
          }
        }; // class ConstantFunction::Evaluator<...>

      private:
        Real _value;

      public:
        explicit ConstantFunction(Real value = Real(0)) :
          _value(value)
        {
        }
      }; // class ConstantFunction
    } // namespace Common
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_COMMON_FUNCTIONS_HPP
