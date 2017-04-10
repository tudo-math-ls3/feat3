#pragma once
#ifndef KERNEL_ANALYTIC_COMMON_HPP
#define KERNEL_ANALYTIC_COMMON_HPP 1

// includes, FEAT
#include <kernel/analytic/static_wrapper.hpp>
#include <kernel/util/math.hpp>

// includes, system
#include <initializer_list>
#include <deque>
#include <vector>

namespace FEAT
{
  namespace Analytic
  {
    /**
     * \brief Analytic Common namespace
     *
     * This namespace encapsulated commonly used functions,
     * which are often used in standard benchmark problems.
     */
    namespace Common
    {
      /**
       * \brief Wrapper class for Tensor-Product scalar static functions
       *
       * This class can be used as a wrapper around another class implementing the
       * StaticFunction interface for 1D. This class will then implement the 2D and 3D
       * by using the tensor-product approach.
       *
       * \tparam Scalar_
       * The name of the class implementing the StaticFunction interface for 1D.
       *
       * \tparam DataType_
       * The data-type used by Scalar_.
       *
       * \author Peter Zajac
       */
      template<typename Scalar_, typename DataType_>
      class TensorStatic
      {
      public:
        /// 1D: function value
        static DataType_ eval(DataType_ x)
        {
          return Scalar_::eval(x);
        }

        /// 1D: X-derivative
        static DataType_ der_x(DataType_ x)
        {
          return Scalar_::der_x(x);
        }

        /// 1D: XX-derivative
        static DataType_ der_xx(DataType_ x)
        {
          return Scalar_::der_xx(x);
        }

        /// 2D: function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return Scalar_::eval(x) * Scalar_::eval(y);
        }

        /// 2D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return Scalar_::der_x(x) * Scalar_::eval(y);
        }

        /// 2D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return Scalar_::eval(x) * Scalar_::der_x(y);
        }

        /// 3D: function value
        static DataType_ eval(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::eval(x) * Scalar_::eval(y) * Scalar_::eval(z);
        }

        /// 3D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::der_x(x) * Scalar_::eval(y) * Scalar_::eval(z);
        }

        /// 3D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::eval(x) * Scalar_::der_x(y) * Scalar_::eval(z);
        }

        /// 3D: Z-derivative
        static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::eval(x) * Scalar_::eval(y) * Scalar_::der_x(z);
        }

        /// 2D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return Scalar_::der_xx(x) * Scalar_::eval(y);
        }

        /// 2D: YY-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y)
        {
          return Scalar_::eval(x) * Scalar_::der_xx(y);
        }

        /// 2D: XY-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y)
        {
          return Scalar_::der_x(x) * Scalar_::der_x(y);
        }

        /// 2D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return Scalar_::der_x(x) * Scalar_::der_x(y);
        }

        /// 3D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::der_xx(x) * Scalar_::eval(y) * Scalar_::eval(z);
        }

        /// 3D: YY-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::eval(x) * Scalar_::der_xx(y) * Scalar_::eval(z);
        }

        /// 3D: ZZ-derivative
        static DataType_ der_zz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::eval(x) * Scalar_::eval(y) * Scalar_::der_xx(z);
        }

        /// 3D: XY-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::der_x(x) * Scalar_::der_x(y) * Scalar_::eval(z);
        }

        /// 3D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::der_x(x) * Scalar_::der_x(y) * Scalar_::eval(z);
        }

        /// 3D: XZ-derivative
        static DataType_ der_xz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::der_x(x) * Scalar_::eval(y) * Scalar_::der_x(z);
        }

        /// 3D: ZX-derivative
        static DataType_ der_zx(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::der_x(x) * Scalar_::eval(y) * Scalar_::der_x(z);
        }

        /// 3D: YZ-derivative
        static DataType_ der_yz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::eval(x) * Scalar_::der_x(y) * Scalar_::der_x(z);
        }

        /// 3D: ZY-derivative
        static DataType_ der_zy(DataType_ x, DataType_ y, DataType_ z)
        {
          return Scalar_::eval(x) * Scalar_::der_x(y) * Scalar_::der_x(z);
        }
      }; // class TensorStatic<...>

      /**
       * \brief Sine-Tensor Static function
       *
       * This class implements the StaticFunction interface representing the function
       *   - 1D: u(x)     = sin(k*pi*x)
       *   - 2D: u(x,y)   = sin(k*pi*x) * sin(k*pi*y)
       *   - 3D: u(x,y,z) = sin(k*pi*x) * sin(k*pi*y) * sin(k*pi*z)
       *
       * For any positive integral \e k, these functions are eigenfunctions of the Laplace operator.
       * The corresponding eigenvalue is \f$ \lambda = -d (k\pi)^2 \f$, where \e d is the dimension of the domain.
       *
       * Moreover, on any rectangular/rectoid domain [x0,x1]x[y0,y1]x[z0,z1] with integral domain
       * boundaries xi,yi,zi, this function fulfills homogeneous Dirichlet boundary conditions.
       *
       * \author Peter Zajac
       */
      template<typename DataType_, int k_ = 1>
      class SineTensorStatic
      {
        static_assert(k_ > 0, "parameter k_ must be a positive integer");

      public:
        /// returns the constant pi
        static DataType_ kpi()
        {
          return DataType_(k_) * Math::pi<DataType_>();
        }

        /// 1D: function value
        static DataType_ eval(DataType_ x)
        {
          return Math::sin(kpi() * x);
        }

        /// 2D: function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return Math::sin(kpi() * x) * Math::sin(kpi() * y);
        }

        /// 3D: function value
        static DataType_ eval(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sin(kpi() * x) * Math::sin(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 1D: X-derivative
        static DataType_ der_x(DataType_ x)
        {
          return kpi() * Math::cos(kpi() * x);
        }

        /// 2D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return kpi() * Math::cos(kpi() * x) * Math::sin(kpi() * y);
        }

        /// 2D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return kpi() * Math::sin(kpi() * x) * Math::cos(kpi() * y);
        }

        /// 3D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z)
        {
          return kpi() * Math::cos(kpi() * x) * Math::sin(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z)
        {
          return kpi() * Math::sin(kpi() * x) * Math::cos(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: Z-derivative
        static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z)
        {
          return kpi() * Math::sin(kpi() * x) * Math::sin(kpi() * y) * Math::cos(kpi() * z);
        }

        /// 1D: XX-derivative
        static DataType_ der_xx(DataType_ x)
        {
          return -Math::sqr(kpi()) * Math::sin(kpi() * x);
        }

        /// 2D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return -Math::sqr(kpi()) * Math::sin(kpi() * x) * Math::sin(kpi() * y);
        }

        /// 2D: YY-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y)
        {
          return der_xx(x, y);
        }

        /// 2D: XY-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y)
        {
          return Math::sqr(kpi()) * Math::cos(kpi() * x) * Math::cos(kpi() * y);
        }

        /// 2D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x, y);
        }

        /// 3D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y, DataType_ z)
        {
          return -Math::sqr(kpi()) * Math::sin(kpi() * x) * Math::sin(kpi() * y) * Math::sin(kpi() * z);
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
          return Math::sqr(kpi()) * Math::cos(kpi() * x) * Math::cos(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xy(x, y, z);
        }

        /// 3D: XZ-derivative
        static DataType_ der_xz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(kpi()) * Math::cos(kpi() * x) * Math::sin(kpi() * y) * Math::cos(kpi() * z);
        }

        /// 3D: ZX-derivative
        static DataType_ der_zx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xz(x, y, z);
        }

        /// 3D: YZ-derivative
        static DataType_ der_yz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(kpi()) * Math::sin(kpi() * x) * Math::cos(kpi() * y) * Math::cos(kpi() * z);
        }

        /// 3D: ZY-derivative
        static DataType_ der_zy(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_yz(x, y, z);
        }
      }; // class SineTensorStatic<...>

      /// \cond internal
      template<typename DataType_>
      using SineBubbleStatic = SineTensorStatic<DataType_, 1>;
      /// \endcond

      /**
       * \brief Sine-Bubble Analytic function
       *
       * This class implements the AnalyticFunction interface representing the function
       *   - 1D: u(x)     = sin(pi*x)
       *   - 2D: u(x,y)   = sin(pi*x) * sin(pi*y)
       *   - 3D: u(x,y,z) = sin(pi*x) * sin(pi*y) * sin(pi*z)
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * This function fulfills homogeneous Dirichlet boundary conditions on the unit-cube domain.
       *
       * \author Peter Zajac
       */
      template<int dim_>
      using SineBubbleFunction = StaticWrapperFunction<dim_, SineBubbleStatic, true, true, true>;

      /**
       * \brief Cosine-Tensor Static function
       *
       * This class implements the StaticFunction interface representing the function
       *   - 1D: u(x)     = cos(k*pi*x)
       *   - 2D: u(x,y)   = cos(k*pi*x) * cos(k*pi*y)
       *   - 3D: u(x,y,z) = cos(k*pi*x) * cos(k*pi*y) * cos(k*pi*z)
       *
       * For any positive integral \e k, these functions are eigenfunctions of the Laplace operator.
       * The corresponding eigenvalue is \f$ \lambda = -d (k\pi)^2 \f$, where \e d is the dimension of the domain.
       *
       * Moreover, on any rectangular/rectoid domain [x0,x1]x[y0,y1]x[z0,z1] with integral domain
       * boundaries xi,yi,zi, this function fulfills homogeneous Neumann boundary conditions including
       * the integral-mean condition \f$ int_\Omega u = 0 \f$.
       *
       * \author Peter Zajac
       */
      template<typename DataType_, int k_ = 1>
      class CosineTensorStatic
      {
        static_assert(k_ > 0, "parameter k_ must be a positive integer");

      public:
        /// Returns k times pi
        static DataType_ kpi()
        {
          return DataType_(k_) * Math::pi<DataType_>();
        }

        /// 1D: function value
        static DataType_ eval(DataType_ x)
        {
          return Math::cos(kpi() * x);
        }

        /// 2D: function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return Math::cos(kpi() * x) * Math::cos(kpi() * y);
        }

        /// 3D: function value
        static DataType_ eval(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::cos(kpi() * x) * Math::cos(kpi() * y) * Math::cos(kpi() * z);
        }

        /// 1D: X-derivative
        static DataType_ der_x(DataType_ x)
        {
          return -kpi() * Math::sin(kpi() * x);
        }

        /// 2D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return -kpi() * Math::sin(kpi() * x) * Math::cos(kpi() * y);
        }

        /// 2D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return -kpi() * Math::cos(kpi() * x) * Math::sin(kpi() * y);
        }

        /// 3D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z)
        {
          return -kpi() * Math::sin(kpi() * x) * Math::cos(kpi() * y) * Math::cos(kpi() * z);
        }

        /// 3D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z)
        {
          return -kpi() * Math::cos(kpi() * x) * Math::sin(kpi() * y) * Math::cos(kpi() * z);
        }

        /// 3D: Z-derivative
        static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z)
        {
          return -kpi() * Math::cos(kpi() * x) * Math::cos(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 1D: XX-derivative
        static DataType_ der_xx(DataType_ x)
        {
          return -Math::sqr(kpi()) * Math::cos(kpi() * x);
        }

        /// 2D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return -Math::sqr(kpi()) * Math::cos(kpi() * x) * Math::cos(kpi() * y);
        }

        /// 2D: YY-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y)
        {
          return der_xx(x, y);
        }

        /// 2D: XY-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y)
        {
          return Math::sqr(kpi()) * Math::sin(kpi() * x) * Math::sin(kpi() * y);
        }

        /// 2D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x, y);
        }

        /// 3D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y, DataType_ z)
        {
          return -Math::sqr(kpi()) * Math::cos(kpi() * x) * Math::cos(kpi() * y) * Math::cos(kpi() * z);
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
          return Math::sqr(kpi()) * Math::sin(kpi() * x) * Math::sin(kpi() * y) * Math::cos(kpi() * z);
        }

        /// 3D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y, DataType_ /*z*/)
        {
          return der_xy(x, y);
        }

        /// 3D: XZ-derivative
        static DataType_ der_xz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(kpi()) * Math::sin(kpi() * x) * Math::cos(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: ZX-derivative
        static DataType_ der_zx(DataType_ x, DataType_ y, DataType_ /*z*/)
        {
          return der_xz(x, y);
        }

        /// 3D: YZ-derivative
        static DataType_ der_yz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(kpi()) * Math::cos(kpi() * x) * Math::sin(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: ZY-derivative
        static DataType_ der_zy(DataType_ x, DataType_ y, DataType_ /*z*/)
        {
          return der_yz(x, y);
        }
      }; // class CosineTensorStatic<...>

      /// \cond internal
      template<typename DataType_>
      using CosineWaveStatic = CosineTensorStatic<DataType_, 1>;
      /// \endcond

      /**
       * \brief Cosine-Wave Analytic function
       *
       * This class implements the AnalyticFunction interface representing the function
       *   - 1D: u(x)     = cos(pi*x)
       *   - 2D: u(x,y)   = cos(pi*x) * cos(pi*y)
       *   - 3D: u(x,y,z) = cos(pi*x) * cos(pi*y) * cos(pi*z)
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * This function fulfills homogeneous Neumann boundary conditions and has vanishing integral
       * mean on the unit-cube domain.
       *
       * \author Peter Zajac
       */
      template<int dim_>
      using CosineWaveFunction = StaticWrapperFunction<dim_, CosineWaveStatic, true, true, true>;

      /**
       * \brief Exponential-Bubble scalar Static function
       *
       * This class implements the StaticFunction interface representing the function
       *   - u(x) = (exp(-(2*x - 1)^2) - exp(-1)) / (exp(-1) - 1)
       *
       * This function fulfills homogeneous Dirichlet boundary conditions on the unit cube domain.
       *
       * \author Peter Zajac
       */
      template<typename DataType_>
      class ExpBubbleScalarStatic
      {
      public:
        /// 1D: function value
        static DataType_ eval(DataType_ x)
        {
          return (Math::exp(-Math::sqr(DataType_(2)*x - DataType_(1))) - Math::exp(-DataType_(1))) / (DataType_(1) - Math::exp(-DataType_(1)));
        }

        /// 1D: X-derivative
        static DataType_ der_x(DataType_ x)
        {
          return (DataType_(8)*x - DataType_(4))*Math::exp(-Math::sqr(DataType_(2)*x - DataType_(1))) / (Math::exp(-DataType_(1)) - DataType_(1));
        }

        /// 1D: XX-derivative
        static DataType_ der_xx(DataType_ x)
        {
          return (DataType_(64)*x*(DataType_(1)-x)-DataType_(8)) * Math::exp(-Math::sqr(DataType_(2)*x - DataType_(1))) / (Math::exp(-DataType_(1)) - DataType_(1));
        }
      }; // class ExpBubbleScalarStatic<...>

      /// \cond internal
      template<typename DataType_>
      using ExpBubbleStatic = TensorStatic<ExpBubbleScalarStatic<DataType_>, DataType_>;
      /// \endcond

      /**
       * \brief Exponential-Bubble Analytic function
       *
       * This class implements the AnalyticFunction interface representing the function
       *   - 1D: u(x)     = (exp(-(2*x - 1)^2) - exp(-1)) / (exp(-1) - 1)
       *   - 2D: u(x,y)   = u(x) * u(y)
       *   - 3D: u(x,y,z) = u(x) * u(y) * u(z)
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * This function fulfills homogeneous Dirichlet boundary conditions on the unit-cube domain.
       *
       * \author Peter Zajac
       */
      template<int dim_>
      using ExpBubbleFunction = StaticWrapperFunction<dim_, ExpBubbleStatic, true, true, true>;

      /**
       * \brief Q2-bubble scalar Static function
       *
       * This class implements the StaticFunction interface representing the function
       *   - u(x) = 4*x*(1-x)
       *
       * This function fulfills homogeneous Dirichlet boundary conditions on the unit-cube domain.
       *
       * \author Peter Zajac
       */
      template<typename DataType_>
      class Q2BubbleScalarStatic
      {
      public:
        /// 1D: function value
        static DataType_ eval(DataType_ x)
        {
          return DataType_(4) * x * (DataType_(1) - x);
        }

        /// 1D: X-derivative
        static DataType_ der_x(DataType_ x)
        {
          return DataType_(4) * (DataType_(1) - DataType_(2) * x);
        }

        /// 1D: XX-derivative
        static DataType_ der_xx(DataType_)
        {
          return -DataType_(8);
        }
      }; // class Q2BubbleScalarStatic<...>

      /// \cond internal
      template<typename DataType_>
      using Q2BubbleStatic = TensorStatic<Q2BubbleScalarStatic<DataType_>, DataType_>;
      /// \endcond

      /**
       * \brief Exponential-Bubble Analytic function
       *
       * This class implements the AnalyticFunction interface representing the function
       *   - 1D: u(x)     =  4*x*(1-x)
       *   - 2D: u(x,y)   = 16*x*(1-x)*y*(1-y)
       *   - 3D: u(x,y,z) = 64*x*(1-x)*y*(1-y)*z*(1-z)
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * This function fulfills homogeneous Dirichlet boundary conditions on the unit-cube domain.
       *
       * \author Peter Zajac
       */
      template<int dim_>
      using Q2BubbleFunction = StaticWrapperFunction<dim_, Q2BubbleStatic, true, true, true>;
      /**
       * \brief Constant Analytic function
       *
       * This class implements the AnalyticFunction interface representing a constant function.
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \author Peter Zajac
       */
      template<int dim_, typename DataType_ = Real>
      class ConstantFunction :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;

        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// Function that is being evaluated
          const DataType _value;

        public:
          /// Constructor
          explicit Evaluator(const ConstantFunction& function) :
            _value(DataType(function._value))
          {
          }


          void value(ValueType& val, const PointType& DOXY(point))
          {
            val = _value;
          }

          void gradient(GradientType& grad, const PointType& DOXY(point))
          {
            grad = DataType(0);
          }

          void hessian(HessianType& hess, const PointType& DOXY(point))
          {
            hess = DataType(0);
          }
        }; // class ConstantFunction::Evaluator<...>

      private:
        /// Value of the constant function
        DataType_ _value;

      public:
        /// Constructor, value defaults to 0
        explicit ConstantFunction(DataType_ value = DataType_(0)) :
          _value(value)
        {
        }
      }; // class ConstantFunction

      /**
       * \brief Analytic distance function
       *
       * This class implements the AnalyticFunction interface representing the distance function
       * \f$ f(x) =\| x - x_0 \|_2 \f$
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \warning Because the function is differentiable everywhere except for x_0, Bad Things (TM) might happen if
       * someone wants to compute the gradient or hessian there, as the functions return 0.
       *
       * \tparam ImgPointType_
       * The type of \f$ x_0 \f$.
       *
       * \author Jordi Paul
       */
      template<int dim_, typename DataType_>
      class DistanceFunction :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// our origin
          const PointType _origin;

        public:
          /// Constructor
          explicit Evaluator(const DistanceFunction& function) :
            _origin(function._origin)
          {
          }

          void value(ValueType& val, const PointType& point) const
          {
            val = (point - _origin).norm_euclid();
          }

          void gradient(GradientType& grad, const PointType& point) const
          {
            DataType norm (DataType(0));
            this->value(norm, point);

            grad.format();
            if(norm > Math::eps<DataType>())
            {
              grad = (DataType(1) / norm) * (point - _origin);
            }
          }

          void hessian(HessianType& hess, const PointType& point) const
          {
            DataType norm (DataType(0));
            this->value(norm, point);

            hess.format();
            if(norm > Math::eps<DataType>())
            {
              norm = DataType(1)/norm;
              DataType denom = Math::sqr(norm)*norm;

              for(Index i(0); i < dim_; ++i)
              {
                hess[i][i] = norm;
                for(Index j(0); j < dim_; ++i)
                {
                  hess[i][j] = ((point[i] - _origin[i]) * (point[j] - _origin[j]))*denom;
                }
              }
            }
          }
        }; // class DistanceFunction::Evaluator<...>

      public:
        /// Point to calculate the distance from
        PointType _origin;

      public:
        /// Constructor
        explicit DistanceFunction(const PointType& origin) :
          _origin(origin)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType& origin)
        {
          _origin = origin;
        }
      }; // class DistanceFunction

      /**
       * \brief Scaled and displaced analytic distance function
       *
       * This class implements the AnalyticFunction interface representing the scaled and displaced distance function
       * \f$ f(x) = a + b \| x - x_0 \|_2 \f$
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \warning Because the function is differentiable everywhere except for x_0, Bad Things (TM) might happen if
       * someone wants to compute the gradient or hessian there, as the functions return 0.
       *
       * \tparam ImgPointType_
       * The type of \f$ x_0 \f$.
       *
       * \author Jordi Paul
       */
      template<int dim_, typename DataType_>
      class DistanceFunctionSD :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;

        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// our origin
          const PointType _origin;
          /// displacement and scaling factor
          const DataType _a, _b;

        public:
          /// Constructor
          explicit Evaluator(const DistanceFunctionSD& function) :
            _origin(function._origin),
            _a(function._a),
            _b(function._b)
          {
          }

          void value(ValueType& val, const PointType& point) const
          {
            val = _a + _b*(point - _origin).norm_euclid();
          }

          void gradient(GradientType& grad, const PointType& point) const
          {
            grad.format();
            DataType norm ((point - _origin).norm_euclid());

            if(norm <= Math::eps<DataType>())
              return;

            grad = (_b / norm) * (point - _origin);
          }

          void hessian(HessianType& hess, const PointType& point) const
          {
            DataType norm (DataType(0));
            this->value(norm, point);
            if(norm <= Math::eps<DataType>())
              return;

            norm = DataType(1)/norm;
            DataType denom = Math::sqr(_b)*Math::sqr(norm)*norm;

            for(Index i(0); i < dim_; ++i)
            {
              hess[i][i] = norm;
              for(Index j(0); j < dim_; ++i)
              {
                hess[i][j] = ((point[i] - _origin[i]) * (point[j] - _origin[j]))*denom;
              }
            }
          }
        }; // class DistanceFunctionSD::Evaluator<...>

      private:
        /// The point to which the distance to is calculated
        PointType _origin;
        /// Displacement of the function
        DataType_ _a;
        /// Scaling factor
        DataType_ _b;

      public:
        /// Constructor
        explicit DistanceFunctionSD(const PointType& origin, const DataType_ a_, const DataType_ b_) :
          _origin(origin),
          _a(a_),
          _b(b_)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType& origin)
        {
          _origin = origin;
        }
      }; // class DistanceFunctionSD

      /**
       * \brief Scaled and displaced analytic distance from i-th coordinate axis/plane function
       *
       * This class implements the AnalyticFunction interface representing the scaled and displaced distance function
       * \f$ f(x) = b | (x - x_0)_i | \f$
       *
       * This class supports function values and gradients for all dimensions.
       *
       * \tparam component_
       * Index of the coordinate axis/plane.
       *
       * \tparam ImgPointType_
       * The type of \f$ x_0 \f$.
       *
       * \author Jordi Paul
       */
      template<int component_, int dim_, typename DataType_>
      class PlaneDistanceFunctionSD :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;

        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;

        typedef Tiny::Vector<DataType_, dim_> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// our origin
          const PointType _origin;
          /// our scaling
          const DataType_ _b;

        public:
          /// Constructor
          explicit Evaluator(const PlaneDistanceFunctionSD& function) :
            _origin(function._origin),
            _b(function._b)
          {
          }

          void value(ValueType& val, const PointType& point) const
          {
            val = _b * (point[component_] - _origin[component_]);
          }

          void gradient(GradientType& grad, const PointType& point) const
          {
            grad.format();
            DataType norm (DataType(0));
            this->value(norm, point);

            if(norm > Math::eps<DataType>())
            {
              grad[component_] = _b;
            }
          }
        }; // class PlaneDistanceFunctionSD::Evaluator<...>

      private:
        /// The point to which the distance to is calculated
        PointType _origin;
        /// Scaling factor
        DataType_ _b;

      public:
        /// Constructor
        explicit PlaneDistanceFunctionSD(const PointType& origin, const DataType_ b_) :
          _origin(origin),
          _b(b_)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType& origin)
        {
          _origin = origin;
        }
      }; // class PlaneDistanceFunctionSD

      /**
       * \brief Function representing the minimum of two analytic functions
       *
       * This is needed i.e. if there are several objects implicitly defined by their zero level sets. The class
       * in general supports function values, gradients and hessians for all dimensions, depending on the two analytic
       * functions supporting these.
       *
       * \warning As min is non differentiable in general, Bad Things(TM) may happen when computing the gradient
       * and/or hessian where the function values are nearly identical.
       *
       * \tparam AnalyticFunctionType1
       * Type for the first AnalyticFunction
       *
       * \tparam AnalyticFunctionType2
       * Type for the second AnalyticFunction
       *
       * \author Jordi Paul
       */
      template<typename AnalyticFunctionType1, typename AnalyticFunctionType2>
      class MinOfTwoFunctions :
        public Analytic::Function
      {
      public:
        /// ensure that the two functions have the same dimension
        static_assert(AnalyticFunctionType1::domain_dim == AnalyticFunctionType2::domain_dim, "domain dimension mismatch");

        /// ensure that both functions are scalar
        static_assert(std::is_same<typename AnalyticFunctionType1::ImageType, Analytic::Image::Scalar>::value, "invalid image type");
        static_assert(std::is_same<typename AnalyticFunctionType2::ImageType, Analytic::Image::Scalar>::value, "invalid image type");

        /// our domain dimension
        static constexpr int domain_dim = AnalyticFunctionType1::domain_dim;
        /// our image type
        typedef Analytic::Image::Scalar ImageType;

        /// Can compute function values if both AnalyticFunctions can do that
        static constexpr bool can_value = (AnalyticFunctionType1::can_value && AnalyticFunctionType2::can_value);
        /// Can compute the function gradient if both AnalyticFunctions can do that
        static constexpr bool can_grad = (AnalyticFunctionType1::can_grad && AnalyticFunctionType2::can_grad);
        /// Can compute the function hessian if both AnalyticFunctions can do that
        static constexpr bool can_hess = (AnalyticFunctionType1::can_hess && AnalyticFunctionType2::can_hess);

      public:
        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// Evaluator for the first AnalyticFunction
          typename AnalyticFunctionType1::template Evaluator<EvalTraits_> _f1_eval;
          /// Evaluator for the second AnalyticFunction
          typename AnalyticFunctionType2::template Evaluator<EvalTraits_> _f2_eval;

        public:
          /// Constructor
          explicit Evaluator(const MinOfTwoFunctions& function) :
            _f1_eval(function._f1),
            _f2_eval(function._f2)
            {
            }

          void value(ValueType& val, const PointType& point) const
          {
            ValueType val1, val2;
            _f1_eval.value(val1, point);
            _f2_eval.value(val2, point);
            val = Math::min(val1, val2);
          }

          void gradient(GradientType& grad, const PointType& point) const
          {
            ValueType val1, val2;
            _f1_eval.value(val1, point);
            _f2_eval.value(val2, point);

            if(Math::abs(val1-val2) < Math::eps<DataType>())
              grad.format();
            else if(val1 < val2)
              _f1_eval.gradient(grad, point);
            else
              _f2_eval.gradient(grad, point);
          }

          void hessian(HessianType& hess, const PointType& point) const
          {
            ValueType val1, val2;
            _f1_eval.value(val1, point);
            _f2_eval.value(val2, point);

            if(Math::abs(val1-val2) < Math::eps<DataType>())
              hess.format();
            else if(val1 < val2)
              _f1_eval.hessian(hess, point);
            else
              _f2_eval.hessian(hess, point);
          }
        }; // class MinOfTwoFunctions::Evaluator<...>

      private:
        /// The first AnalyticFunction
        const AnalyticFunctionType1& _f1;
        /// The second AnalyticFunction
        const AnalyticFunctionType2& _f2;

      public:
        /// Constructor
        explicit MinOfTwoFunctions(const AnalyticFunctionType1& f1_, const AnalyticFunctionType2& f2_) :
          _f1(f1_),
          _f2(f2_)
        {
        }
      }; // class MinOfTwoFunctions

      /// \cond internal
      /**
       * \brief Heaviside static function
       *
       * This class implements the StaticFunction interface representing the function
       * \f[
       *   u(x) = H(x) =
       *   \begin{cases}
       *     0, & x < 0 \\
       *     1, & x \geq 0
       *   \end{cases}
       * \f]
       * in \f$ 1D \f$ and
       * \f[
       *   u(x,y) = H(x) H(y), u(x,y,z) = H(x) H(y) H(z)
       * \f]
       * in \f$ 2D, 3D \f$ respectively.
       *
       * \author Jordi Paul
       */
      template<typename DataType_>
      class HeavisideStatic
      {
      public:

        /// 1D: function value
        static DataType_ eval(DataType_ x)
        {
          return ( x < DataType_(0) ) ? DataType_(0) : DataType_(1);
        }

        /// 2D: function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return eval(x)*eval(y);
        }

        /// 3D: function value
        static DataType_ eval(DataType_ x, DataType_ y, DataType_ z)
        {
          return eval(x)*eval*(y)*eval(z);
        }

      }; // class HeavisideStatic<...>
      /// \endcond

      /**
       * \brief Heaviside function
       *
       * This class implements the AnalyticFunction interface representing the function
       * \f[
       *   u(x) = H(x) =
       *   \begin{cases}
       *     0, & x < 0 \\
       *     1, & x \geq 0
       *   \end{cases}
       * \f]
       * in \f$ 1D \f$ and
       * \f[
       *   u(x,y) = H(x) H(y), u(x,y,z) = H(x) H(y) H(z)
       * \f]
       * in \f$ 2D, 3D \f$ respectively.
       *
       * \author Jordi Paul
       */
      template<int dim_>
      using HeavisideFunction = StaticWrapperFunction<dim_, HeavisideStatic, true, false, false>;

      /// \cond internal
      /**
       * \brief Regularised Heaviside static function
       *
       * This class implements the StaticFunction interface representing the function
       * \f[
       *   u(x) = H(x) =
       *   \begin{cases}
       *     0, & x < 0 \\
       *     1, & 2 (\cosh(x) - 1) \geq 0
       *   \end{cases}
       * \f]
       *
       * in \f$ 1D \f$ and
       * \f[
       *   u(x,y) = H(x) H(y), u(x,y,z) = H(x) H(y) H(z)
       * \f]
       * in \f$ 2D, 3D \f$ respectively.
       *
       * \author Jordi Paul
       */
      template<typename DataType_>
      class HeavisideRegStatic
      {
      public:

        /// 1D: function value
        static DataType_ eval(DataType_ x)
        {
          return ( x < DataType_(0) ) ? DataType_(0) : DataType_(2)*(Math::cosh(x) - DataType_(1));
        }

        /// 2D: function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return eval(x)*eval(y);
        }

        /// 3D: function value
        static DataType_ eval(DataType_ x, DataType_ y, DataType_ z)
        {
          return eval(x)*eval*(y)*eval(z);
        }

        /// 1D: X-derivative
        static DataType_ der_x(DataType_ x)
        {
          return ( x < DataType_(0) ) ? DataType_(0) : DataType_(2)*Math::sinh(x);
        }

        /// 2D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return der_x(x)*eval(y);
        }

        /// 2D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return eval(x)*der_x(y);
        }

        /// 3D: X-derivative
        static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_x(x)*eval(y)*eval(z);
        }

        /// 3D: Y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z)
        {
          return eval(x)*der_x(y)*eval(z);
        }

        /// 3D: Z-derivative
        static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z)
        {
          return eval(x)*eval(y)*der_x(z);
        }

        /// 1D: XX-derivative
        static DataType_ der_xx(DataType_ x)
        {
          return ( x < DataType_(0) ) ? DataType_(0) : DataType_(2)*Math::cosh(x);
        }

        /// 2D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return der_xx(x)*eval(y);
        }

        /// 2D: YY-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y)
        {
          return eval(x)*der_xx(y);
        }

        /// 2D: XY-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y)
        {
          return der_x(x)*der_x(y);
        }

        /// 2D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x, y);
        }

        /// 3D: XX-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xx(x)*eval(y)*eval(z);
        }

        /// 3D: YY-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y, DataType_ z)
        {
          return eval(x)*der_xx(y)*eval(z);
        }

        /// 3D: ZZ-derivative
        static DataType_ der_zz(DataType_ x, DataType_ y, DataType_ z)
        {
          return eval(x)*eval(y)*der_xx(z);
        }

        /// 3D: XY-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_x(x)*der_x(y)*eval(z);
        }

        /// 3D: YX-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_x(x)*der_x(y)*eval(z);
        }

        /// 3D: XZ-derivative
        static DataType_ der_xz(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_x(x)*eval(y)*der_x(z);
        }

        /// 3D: ZX-derivative
        static DataType_ der_zx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xz(x, y, z);
        }

        /// 3D: YZ-derivative
        static DataType_ der_yz(DataType_ x, DataType_ y, DataType_ z)
        {
          return eval(x)*der_x(y)*der_x(z);
        }

        /// 3D: ZY-derivative
        static DataType_ der_zy(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_yz(x, y, z);
        }

      }; // class HeavisideRegStatic<...>
      /// \endcond

      /**
       * \brief Regularised Heaviside static function
       *
       * This class implements the AnalyticFunction interface representing the function
       * \f[
       *   u(x) = H(x) =
       *   \begin{cases}
       *     0, & x < 0 \\
       *     1, & 2 (\cosh(x) - 1) \geq 0
       *   \end{cases}
       * \f]
       *
       * in \f$ 1D \f$ and
       * \f[
       *   u(x,y)   = H(x)*H(y), u(x,y,z) = H(x)*H(y)*H(z)
       * \f]
       * in \f$ 2D, 3D \f$ respectively.
       *
       * \author Jordi Paul
       */
      template<int dim_>
      using HeavisideRegFunction = StaticWrapperFunction<dim_, HeavisideRegStatic, true, false, false>;

      /**
       * \brief 1D Polynomial function class template
       *
       * This class template implements a scalar polynomial implementing the Assembly::AnalyticFunction
       * interface.
       *
       * \tparam DataType_
       * The datatype of the polynomial coefficients. Must be a floating point type.
       *
       * \author Peter Zajac
       */
      template<typename DataType_>
      class PolynomialFunction1D :
        public Analytic::Function
      {
      public:
        /// this is a 1D function
        static constexpr int domain_dim = 1;
        /// this is a scalar function
        typedef Analytic::Image::Scalar ImageType;
        /// we provide function values
        static constexpr bool can_value = true;
        /// we provide function gradients
        static constexpr bool can_grad = true;
        /// we provide function hessians
        static constexpr bool can_hess = true;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          typedef typename EvalTraits_::DataType DataType;
          typedef typename EvalTraits_::PointType PointType;
          typedef typename EvalTraits_::ValueType ValueType;
          typedef typename EvalTraits_::GradientType GradientType;
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// our polynomial coefficients
          std::vector<DataType> _coeff;

        public:
          /// Constructor
          explicit Evaluator(const PolynomialFunction1D& function)
          {
            for(auto it = function._coeff.begin(); it != function._coeff.end(); ++it)
              _coeff.push_back(DataType(*it));
          }

          void value(ValueType& val, const PointType& point) const
          {
            // evaluate polynomial via horner scheme
            DataType x = point[0];
            DataType y = DataType(0);
            for(std::size_t k(_coeff.size()); k > std::size_t(0); )
              y = x * y + _coeff[--k];

            val = y;
          }

          void gradient(GradientType& grad, const PointType& point) const
          {
            grad.format();
            std::size_t k = _coeff.size();
            if(k <= std::size_t(0))
              return;

            // evaluate polynomial via horner scheme
            DataType x = point[0];
            DataType y = DataType(0);
            for( ; (--k) > std::size_t(0); )
              y = x * y + (_coeff[k] * DataType(k));

            grad[0] = y;
          }

          void hessian(HessianType& hess, const PointType& point) const
          {
            hess.format();
            std::size_t k = _coeff.size();
            if(k <= std::size_t(1))
              return;

            // evaluate polynomial via horner scheme
            DataType x = point[0];
            DataType y = DataType(0);
            for( ; (--k) > std::size_t(1); )
              y = x * y + (_coeff[k] * DataType(k*(k-1)));

            hess[0][0] = y;
          }
        }; // class PolynomialFunction1D::Evaluator

      private:
        /// the polynomial coefficient vector
        std::vector<DataType_> _coeff;

      public:
        /// default constructor
        PolynomialFunction1D() :
          _coeff()
        {
        }

        /**
         * \brief Polynomial coefficient constructor
         *
         * \param[in] coeff
         * An initializer list containing the coefficients of the polynomial in ascending monomial degrees.
         */
        explicit PolynomialFunction1D(std::initializer_list<DataType_> coeff) :
          _coeff(coeff)
        {
        }

        /**
         * \brief Sets a monomial coefficient.
         *
         * \param[in] degree
         * The degree of the monomial whose coefficient is to be set.
         *
         * \param[in] coeff
         * The coefficient of the monomial to be set.
         *
         * \returns
         * The degree of the polynomial + 1.
         */
        Index set_coeff(Index degree, DataType_ coeff)
        {
          // check degree
          if(Index(_coeff.size()) <= degree)
            _coeff.resize(std::size_t(degree+1), DataType_(0));

          // set coefficient
          _coeff.at(std::size_t(degree)) = coeff;

          // return degree+1
          return Index(_coeff.size());
        }
      }; // class PolynomialFunction1D

      /// \cond internal
      /**
       * \brief Bazaraa/Shetty function
       *
       * \tparam DataType_
       * Floating point precision
       *
       * This class implements the StaticFunction interface representing the function
       *  \f[
       *    u(x,y) = (x - 2)^4 + (x - 2 y)^2
       *  \f]
       *  This function has a global miminum in \f$ x_0 = (2, 1)^T, u(x_0) = 0 \f$ and is frequently used for
       *  testing optimisation algorithms because the hessian at the minimal point is singular.
       *
       * \author Jordi Paul
       */
      template<typename DataType_>
      class BazaraaShettyStatic
      {
      public:
        /// function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return Math::sqr(Math::sqr(x - DataType_(2))) + Math::sqr(x - DataType_(2)*y);
        }

        /// x-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return DataType_(4)*(x - DataType_(2))*Math::sqr(x - DataType_(2)) + DataType_(2)*(x - DataType_(2)*y);
        }

        /// y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return DataType_(4)*(DataType_(2)*y - x);
        }

        /// xx-derivative
        static DataType_ der_xx(DataType_ x, DataType_ DOXY(y))
        {
          return DataType_(12)*Math::sqr(x - DataType_(2)) + DataType_(2);
        }

        /// yy-derivative
        static DataType_ der_yy(DataType_ DOXY(x), DataType_ DOXY(y))
        {
          return DataType_(8);
        }

        /// xy-derivative
        static DataType_ der_xy(DataType_ DOXY(x), DataType_ DOXY(y))
        {
          return DataType_(4);
        }

        /// yx-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x,y);
        }

      }; // class BazaraaShettyStatic<...>
      /// \endcond

      /// \cond internal
      /**
       * \brief Goldstein/Price function
       *
       * \tparam DataType_
       * Floating point precision
       *
       * This class implements the StaticFunction interface representing the function
       *  \f[
       *    u(x,y) =(1+(1+x+y)^2 * (19 - 14*x+3*x^2-14*y+6*x*y+3*y^2))
       *    *(30 + (2*x-3*y)^2 * (18 - 32*x +12*x^2 + 48*y - 36*x*y + 27*y^2))
       *  \f]
       *  This function has a global miminum in \f$ x_0 = (0, -0.5)^T, u(x_0) = 3 \f$, three more local minima
       *  \f[
       *    x_1 = (-0.6, -0.4)^T, x_2 = (1.2, 0.8)^T, x_3 = (1.8, 0.2),
       *  \f]
       *  one saddle point \f$ x_4 = (-0.4, 0.6)^T \f$ and a local maximum in \f$ x_5 = (1.2, -0.2)^T.
       *
       *
       * \author Jordi Paul
       */
      template<typename DataType_>
      class GoldsteinPriceStatic
      {
      public:
        /// function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (30 + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)));
        }

        /// x-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (30 + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (DataType_(4) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(24) * x - DataType_(36) * y - DataType_(32)));
        }

        /// y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (30 + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (-DataType_(6) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48)));
        }

        /// xx-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return DataType_((DataType_(6) * Math::pow(x, DataType_(2)) + DataType_(12) * x * y + DataType_(6) * Math::pow(y, DataType_(2)) - DataType_(28) * x - DataType_(28) * y + DataType_(38) + DataType_(4) * (DataType_(1) + x + y) * (DataType_(6) * x + DataType_(6) * y - DataType_(14)) + DataType_(6) * Math::pow(DataType_(1) + x + y, DataType_(2))) * (30 + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + DataType_(2) * (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (DataType_(4) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(24) * x - DataType_(36) * y - DataType_(32))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (DataType_(96) * Math::pow(x, DataType_(2)) - DataType_(288) * x * y + DataType_(216) * Math::pow(y, DataType_(2)) - DataType_(256) * x + DataType_(384) * y + DataType_(144) + DataType_(8) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(24) * x - DataType_(36) * y - DataType_(32)) + DataType_(24) * Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2))));
        }

        /// yy-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y)
        {
          return (DataType_(6) * Math::pow(x, DataType_(2)) + DataType_(12) * x * y + DataType_(6) * Math::pow(y, DataType_(2)) - DataType_(28) * x - DataType_(28) * y + DataType_(38) + DataType_(4) * (DataType_(1) + x + y) * (DataType_(6) * x + DataType_(6) * y - DataType_(14)) + DataType_(6) * Math::pow(DataType_(1) + x + y, DataType_(2))) * (30 + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + DataType_(2) * (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (-DataType_(6) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (DataType_(216) * Math::pow(x, DataType_(2)) - DataType_(648) * x * y + DataType_(486) * Math::pow(y, DataType_(2)) - DataType_(576) * x + DataType_(864) * y + DataType_(324) - DataType_(12) * (DataType_(2) * x - DataType_(3) * y) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48)) + DataType_(54) * Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)));
        }

        /// xy-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y)
        {
          return (DataType_(6) * Math::pow(x, DataType_(2)) + DataType_(12) * x * y + DataType_(6) * Math::pow(y, DataType_(2)) - DataType_(28) * x - DataType_(28) * y + DataType_(38) + DataType_(4) * (DataType_(1) + x + y) * (DataType_(6) * x + DataType_(6) * y - DataType_(14)) + DataType_(6) * Math::pow(DataType_(1) + x + y, DataType_(2))) * (30 + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (DataType_(4) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(24) * x - DataType_(36) * y - DataType_(32))) + (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (-DataType_(6) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (-DataType_(144) * Math::pow(x, DataType_(2)) + DataType_(432) * x * y - DataType_(324) * Math::pow(y, DataType_(2)) + DataType_(384) * x - DataType_(576) * y - DataType_(216) - DataType_(6) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(24) * x - DataType_(36) * y - DataType_(32)) + DataType_(4) * (DataType_(2) * x - DataType_(3) * y) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48)) - DataType_(36) * Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)));
        }

        /// xy-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x, y);
        }

      }; // class GoldsteinPriceStatic<...>
      /// \endcond

      /**
       * \brief Goldstein/Price function
       *
       * \tparam DataType_
       * Floating point precision
       *
       * This class implements the StaticFunction interface representing the function
       *  \f[
       *    u(x,y) =(1+(1+x+y)^2 * (19 - 14*x+3*x^2-14*y+6*x*y+3*y^2))
       *    *(30 + (2*x-3*y)^2 * (18 - 32*x +12*x^2 + 48*y - 36*x*y + 27*y^2))
       *  \f]
       *  This function has a global miminum in \f$ x_0 = (0, -0.5)^T, u(x_0) = 3 \f$.
       *
       * \author Jordi Paul
       */
      using GoldsteinPriceFunction = StaticWrapperFunction<2, GoldsteinPriceStatic, true, true, true>;

      /**
       * \brief Bazaraa/Shetty function
       *
       * This class implements the AnalyticFunction interface representing the function
       *  \f[
       *    u(x,y) = (x - 2)^4 + (x - 2 y)^2
       *  \f]
       *  This function has a global miminum in \f$ x_0 = (2, 1)^T, u(x_0) = 0 \f$ and is frequently used for
       *  testing optimisation algorithms because the hessian at the minimal point is singular.
       *
       * \author Jordi Paul
       */
      using BazaraaShettyFunction = StaticWrapperFunction<2, BazaraaShettyStatic, true, true, true>;

      /// \cond internal
      /**
       * \brief Himmelblau function
       *
       * \tparam DataType_
       * Floating point precision
       *
       * This class implements the StaticFunction interface representing the function
       *  \f[
       *    u(x,y) = (x^2 + y^2 - 11)^2 + (x + y^2 - 7)^2
       *  \f]
       *  The function is nonconvex and has 4 local minima at
       *  \f{align*}{
       *    x_0 & \approx (-3.77931025337774689189076584129, -3.28318599128616941226600051437)^T \\
       *    x_1 & \approx (-2.80511808695274485305357239809, 3.13131251825057296580430072341)^T \\
       *    x_2 & = (3, 2)^T \\
       *    x_3 & \approx (3.58442834033049174494433823938, -1.84812652696440355353830020904)^T.
       *  \f}
       *  with \f$ \forall i = 0, \dots, 3: u(x_i) = 0 \f$.
       *
       *  It is often used for testing optimisation algorithms because of the nonconvexity and existence of a local
       *  maximum in \f$ x_4 \approx ( -0.270845, -0.93039)^T \f$ and several saddle points.
       *
       * \author Jordi Paul
       */
      template<typename DataType_>
      class HimmelblauStatic
      {
      public:
        /// function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return Math::sqr(Math::sqr(x) + y - DataType_(11))
            + Math::sqr( x + Math::sqr(y) - DataType_(7));
        }

        /// x-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return DataType_(4)*x*(Math::sqr(x) + y - DataType_(11)) + DataType_(2)*(x + Math::sqr(y) - DataType_(7));
        }

        /// y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return DataType_(2)* (Math::sqr(x) + y - DataType_(11)) + DataType_(4)*y*(x + Math::sqr(y) - DataType_(7));
        }

        /// xx-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return DataType_(12)*Math::sqr(x) + DataType_(4)*y - DataType_(42);
        }

        /// yy-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y)
        {
          return DataType_(4)*x + DataType_(12)*Math::sqr(y) - DataType_(26);
        }

        /// xy-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y)
        {
          return DataType_(4)*(x + y);
        }

        /// yx-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x,y);
        }

      }; // class HimmelblauStatic<...>
      /// \endcond

      /**
       * \brief Himmelblau function
       *
       * This class implements the AnalyticFunction interface representing the function
       *  \f[
       *    u(x,y) = (x^2 + y^2 - 11)^2 + (x + y^2 - 7)^2
       *  \f]
       *  The function is nonconvex and has 4 local minima at
       *  \f{align*}{
       *    x_0 & \approx (-3.77931025337774689189076584129, -3.28318599128616941226600051437)^T \\
       *    x_1 & \approx (-2.80511808695274485305357239809, 3.13131251825057296580430072341)^T \\
       *    x_2 & = (3, 2)^T \\
       *    x_3 & \approx (3.58442834033049174494433823938, -1.84812652696440355353830020904)^T.
       *  \f}
       *  with \f$ \forall i = 0, \dots, 3: u(x_i) = 0 \f$.
       *
       *  It is often used for testing optimisation algorithms because of the nonconvexity and existence of a local
       *  maximum in \f$ x_4 \approx ( -0.270845, -0.93039)^T \f$ and several saddle points.
       *
       * \author Jordi Paul
       */
      using HimmelblauFunction = StaticWrapperFunction<2, HimmelblauStatic, true, true, true>;

      /// \cond internal
      template<typename DataType_>
      class RosenbrockStatic
      {
      public:
        /// function value
        static DataType_ eval(DataType_ x, DataType_ y)
        {
          return DataType_(100)*Math::sqr( y - Math::sqr(x)) + Math::sqr(DataType_(1) - x);
        }

        /// x-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return -DataType_(400)*x*(y - Math::sqr(x) ) + DataType_(2)*(x - DataType_(1));
        }

        /// y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return DataType_(200)*(y - Math::sqr(x));
        }

        /// xx-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return DataType_(1200)*Math::sqr(x) - DataType_(400)*y + DataType_(2);
        }

        /// yy-derivative
        static DataType_ der_yy(DataType_ DOXY(x), DataType_ DOXY(y))
        {
          return DataType_(200);
        }

        /// xy-derivative
        static DataType_ der_xy(DataType_ x, DataType_ DOXY(y))
        {
          return -DataType_(400)*x;
        }

        /// yx-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x,y);
        }

      }; // class RosenbrockStatic<...>
      /// \endcond

      /**
       * \brief Rosenbrock function
       *
       * This class implements the AnalyticFunction interface representing the function
       * \f[
       *   u(x,y) = 100(y-x^2)^2 + (1-x)^2.
       * \f]
       *
       * The function has a global minimum in \f$ x_0 = (1, 1)^T\f$ and a "steep valley" along the parabola
       * \f$ y = x^2 \f$. This is a great challenge to descend-based optimisation algorithms like steepest descent or
       * nonlinear CG and the reason it is frequently used as a target function to test such algorithms.
       *
       * \author Jordi Paul
       */
      using RosenbrockFunction = StaticWrapperFunction<2, RosenbrockStatic, true, true, true>;

      /**
       * \brief 2D Parabolic Profile function base-class
       *
       * This class represents a parabolic profile along a 2D line segment, which is given
       * by the coordinates (x_0, y_0) and (x_1, y_1).
       *
       * More precisely: For a given point (x,y) let \f$ t\in(0,1) \f$ denote the interpolation
       * parameter of the orthonal projection of (x,y) onto the line segment, then the
       * parabolic profile function value is given by \f$4 v_{max} t (t-1)\f$.
       *
       * \author Peter Zajac
       */
      class ParProfileBase :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;

      protected:
        // coordinates of line segment
        Real _x0, _y0, _x1, _y1;
        // maximum value
        Real _vmax;

      public:
        /**
         * \brief Default Constructor
         *
         * This constructor initialises a parabolic profile along the segment (0,0)-(0,1)
         * with maximum value 1.
         */
        ParProfileBase() :
          _x0(0.0), _y0(0.0), _x1(0.0), _y1(1.0), _vmax(1.0)
        {
        }

        /**
         * \brief Constructor
         *
         * \param[in] x0, y0
         * Coordinates of the first segment point.
         *
         * \param[in] x1, y1
         * Coordinates of the second segment point.
         *
         * \param[in] vmax
         * Maximum value of the parabolic profile.
         */
        explicit ParProfileBase(Real x0, Real y0, Real x1, Real y1, Real vmax = Real(1.0)) :
          _x0(x0), _y0(y0), _x1(x1), _y1(y1), _vmax(vmax)
        {
        }

        /**
         * \brief Parses the profile parameters from a string.
         *
         * This function can be used to parse the configuration for the parabolic profile
         * from a string, which may e.g. be read from the command line or a text file.
         *
         * The supported syntax for the string is
         * \f[(x_0 y_0 , x_1 y_1)\f]
         * or
         * \f[(x_0 y_0 , x_1 y_1 , v_{max})\f]
         *
         * \param[in] sbc
         * The string to be parsed.
         *
         * \returns
         * \c true, if the string was parsed successfully, otherwise \c false.
         */
        bool parse(const String& sbc)
        {
          auto li = sbc.find_first_of('(');
          auto ri = sbc.find_last_of(')');
          if((li == sbc.npos) || (ri == sbc.npos))
            return false;

          std::deque<String> sv, sv0, sv1;
          sbc.substr(li+1, ri-li-1).split_by_charset(sv, ",");
          if((sv.size() < std::size_t(2)) || (sv.size() > std::size_t(3)))
            return false;

          sv[0].trim().split_by_charset(sv0);
          sv[1].trim().split_by_charset(sv1);
          if(sv0.size() != std::size_t(2)) return false;
          if(sv1.size() != std::size_t(2)) return false;

          if(!sv0[0].parse(_x0)) return false;
          if(!sv0[1].parse(_y0)) return false;
          if(!sv1[0].parse(_x1)) return false;
          if(!sv1[1].parse(_y1)) return false;

          if((sv.size() == std::size_t(3)) && !sv.back().parse(_vmax))
            return false;

          return true;
        }
      }; // class ParProfileBase

      /**
       * \brief 2D Scalar Parabolic Profile function
       *
       * This class represents a scalar parabolic profile along a 2D line segment,
       * which is given by the coordinates (x_0, y_0) and (x_1, y_1).
       *
       * \author Peter Zajac
       */
      class ParProfileScalar :
        public ParProfileBase
      {
      public:
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;

        using ParProfileBase::ParProfileBase;

        template<typename Traits_>
        class Evaluator :
          public Analytic::Function::Evaluator<Traits_>
        {
        protected:
          typedef typename Traits_::DataType DataType;
          typedef typename Traits_::PointType PointType;
          typedef typename Traits_::ValueType ValueType;

          PointType _vo, _ve;
          DataType _den, _vmax;

        public:
          explicit Evaluator(const ParProfileScalar& function)
          {
            _vo[0] = DataType(function._x0);
            _vo[1] = DataType(function._y0);
            _ve[0] = DataType(function._x1 - function._x0);
            _ve[1] = DataType(function._y1 - function._y0);
            _den = DataType(1) / Tiny::dot(_ve, _ve);
            _vmax = DataType(function._vmax);
          }

          void value(ValueType& val, const PointType& point)
          {
            // project point onto line segment
            DataType x = Math::clamp(Tiny::dot(point - _vo, _ve) * _den, DataType(0), DataType(1));
            // compute function value
            val = _vmax * DataType(4) * x * (DataType(1) - x);
          }
        };
      }; // class ParProfileScalar

      /**
       * \brief 2D Vector-Valued Parabolic Profile function
       *
       * This class represents a parabolic profile vector field along a 2D line segment,
       * which is given by the coordinates (x_0, y_0) and (x_1, y_1).
       *
       * \author Peter Zajac
       */
      class ParProfileVector :
        public ParProfileBase
      {
      public:
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Vector<2> ImageType;
        static constexpr bool can_value = true;

        using ParProfileBase::ParProfileBase;

        template<typename Traits_>
        class Evaluator :
          public Analytic::Function::Evaluator<Traits_>
        {
        protected:
          typedef typename Traits_::DataType DataType;
          typedef typename Traits_::PointType PointType;
          typedef typename Traits_::ValueType ValueType;

          PointType _vo, _ve, _vn;
          DataType _den, _vmax;

        public:
          explicit Evaluator(const ParProfileVector& function)
          {
            _vo[0] = DataType(function._x0);
            _vo[1] = DataType(function._y0);
            _ve[0] = DataType(function._x1 - function._x0);
            _ve[1] = DataType(function._y1 - function._y0);
            _vn[0] =  _ve[1];
            _vn[1] = -_ve[0];
            _vn.normalise();
            _den = DataType(1) / Tiny::dot(_ve, _ve);
            _vmax = DataType(function._vmax);
          }

          void value(ValueType& val, const PointType& point)
          {
            // project point onto line segment
            DataType x = Math::clamp(Tiny::dot(point - _vo, _ve) * _den, DataType(0), DataType(1));
            // compute function value
            DataType v = _vmax * DataType(4) * x * (DataType(1) - x);
            val[0] = _vn[0] * v;
            val[1] = _vn[1] * v;
          }
        };
      }; // class ParProfileVector

      template<typename DT_>
        class ExpScalarStatic
        {
          public:
            static constexpr DT_ p = DT_(10);

            static DT_ eval(DT_ x)
            {
              return (Math::exp(p) - Math::exp(p*x*x)) / (Math::exp(p) - DT_(1));
            }

            static DT_ der_x(DT_ x)
            {
              return -DT_(2)*p*x*Math::exp(p*x*x)/(Math::exp(p)-DT_(1));
            }

            static DT_ der_xx(DT_ x)
            {
              return -DT_(2)*p*Math::exp(p*x*x)*(DT_(2)*p*x*x+DT_(1))/(Math::exp(p)-DT_(1));
            }
        };

      template<typename DataType_>
        using ExpStatic = Analytic::Common::TensorStatic<ExpScalarStatic<DataType_>, DataType_>;

      template<int dim_>
        using ExpFunction = Analytic::StaticWrapperFunction<dim_, ExpStatic, true, true, true>;

      /**
       * \brief Velocity field for a rigid body rotation in the x,y plane
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * This is the mapping
       * \f[
       *    f: \mathbb{R}^d \to \mathbb{R}^d, (x, y, z)^T \mapsto \omega (-(y - y_0), x - x_0, 0)
       * \f]
       * for a given angular velocity \f$ \omega \f$ and and origin \f$ (x_0, y_0, z_0)^T \f$.
       *
       * \author Jordi Paul
       */
      template<typename DT_, int dim_>
      class XYPlaneRotation : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<dim_> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = dim_;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can compute the Hessian
        static constexpr bool can_hess = true;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// The angular velocity
          const DataType _angular_velocity;
          /// our origin
          const PointType& _origin;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] function
           * The x,y plane velocity field function
           */
          explicit Evaluator(const XYPlaneRotation& function) :
            _angular_velocity(function._angular_velocity),
            _origin(function._origin)
          {
          }

          /**
           * \brief Computes the value
           *
           * \param[out] val
           * The (vector valued) function value
           *
           * \param[in] point
           * The domain point
           */
          void value(ValueType& val, const PointType& point)
          {
            val.format();
            val[0] = _angular_velocity*(-(point[1] - _origin[1]));
            val[1] = _angular_velocity*( (point[0] - _origin[0]));
          }

          /**
           * \brief Computes the value
           *
           * \param[out] grad
           * The (matrix valued) gradient
           *
           * \param[in] point
           * The domain point
           */
          void gradient(GradientType& grad, const PointType& DOXY(point)) const
          {
            grad.format();
            grad[0][1] = -DataType(1);
            grad[1][0] = DataType(1);
          }

          /**
           * \brief Computes the Hessian
           *
           * \param[out] hess
           * The (tensor valued) Hessian
           *
           * \param[in] point
           * The domain point
           */
          void hessian(HessianType& hess, const PointType& DOXY(point)) const
          {
            hess.format();
          }

        }; // class XYPlaneRotation::Evaluator<...>

      public:
        /// The angular velocity
        const DataType _angular_velocity;
        /// Point to rotate around
        const PointType _origin;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] angular_velocity
         * The angular velocity to use
         *
         * \param[in] origin
         * The origin around which to rotate
         *
         */
        explicit XYPlaneRotation(const DataType angular_velocity, const PointType& origin) :
          _angular_velocity(angular_velocity),
          _origin(origin)
        {
        }
      }; // class XYPlaneRotation

      /**
       * \brief Velocity field for a parabolic profile in the y,z plane, constant in x
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * This is the mapping
       * \f[
       *    f: \mathbb{R}^d \to \mathbb{R}^d, f_0(x_0, x_1, x_2)^T =
       *    \alpha \frac{2^d}{\prod_{i=0}^{d-1}(b_i - a_i)}\prod_{i=0}^{d-1}( (x_i - a_i)(b_i - x_i) ),
       *    f_i \equiv 0, i=1,\dots,d
       * \f]
       *
       * for a given amplitude \f$ \alpha \f$. This means the function's \f$ x \f$-component is \f$ \alpha \f$ at
       * the midpoint of the rectangle \f$ \{x\}\times [a_1, b_1] \times [a_2, b_2] \f$ and zero on its boundary.
       *
       * \author Jordi Paul
       */
      template<typename DT_, int dim_>
      class YZPlaneParabolic : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<dim_> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = dim_;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// The scaling factor according to the amplitude
          DataType _fac;
          /// The points where the function becomes zero
          const std::vector<Tiny::Vector<DataType, 2>>& _zeros;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] function
           * The x,y plane velocity field function
           */
          explicit Evaluator(const YZPlaneParabolic& function) :
            _fac(function._amplitude),
            _zeros(function._zeros)

          {
            XASSERT(_zeros.size() == size_t(domain_dim-1));

            for(int d(1); d < dim_; ++d)
            {
              _fac *= DataType(2)/Math::sqr(_zeros.at(d-1)[1]-_zeros.at(d-1)[0]);
            }
          }

          /**
           * \brief Computes the value
           *
           * \param[out] val
           * The (vector valued) function value
           *
           * \param[in] point
           * The domain point
           */
          void value(ValueType& val, const PointType& point)
          {
            val.format();
            val(0) = _fac;

            for(Index d(1); d < Index(PointType::n); ++d)
            {
              val(0) *= (point[d] - _zeros.at(d-1)[0])*(_zeros.at(d-1)[1] - point[d]);
            }
          }

          /**
           * \brief Computes the gradient
           *
           * \param[out] grad
           * The (matrix valued) gradient
           *
           * \param[in] point
           * The domain point
           */
          void gradient(GradientType& grad, const PointType& point) const
          {
            grad.format();
            for(int d(1); d < domain_dim; ++d)
            {
              grad[d][0] = _fac[d]*(_zeros.at(d-1)[0] + _zeros.at(d-1)[1] - DataType(2)*point(d));
            }
          }

          /**
           * \brief Computes the Hessian
           *
           * \param[out] hess
           * The (tensor valued) Hessian
           *
           * \param[in] point
           * The domain point
           */
          void hessian(HessianType& hess, const PointType& DOXY(point)) const
          {
            hess.format();
          }

        }; // class YZPlaneParabolic::Evaluator<...>

      public:
        /// The maximum value of the parabolic profile
        const DataType _amplitude;
        /// The points where the function becomes zero
        std::vector<Tiny::Vector<DataType, 2>> _zeros;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] amplitude
         * The amplitude to use
         *
         * \param[in] zeros_y
         * The roots for the y part
         *
         * \param[in] zeros_z
         * The roots for the z part
         *
         */
        explicit YZPlaneParabolic(const DataType amplitude, const PointType& zeros_y) :
          _amplitude(amplitude),
          _zeros(1)
        {
          _zeros.at(0) = zeros_y;
        }

        /**
         * \brief Constructor
         *
         * \param[in] amplitude
         * The amplitude to use
         *
         * \param[in] zeros_y
         * The roots for the y part
         *
         * \param[in] zeros_z
         * The roots for the z part
         *
         */
        explicit YZPlaneParabolic(const DataType amplitude, const PointType& zeros_y, const PointType& zeros_z) :
          _amplitude(amplitude),
          _zeros(2)
        {
          _zeros.at(0) = zeros_y;
          _zeros.at(1) = zeros_z;
        }
      }; // class YZPlaneParabolic

      /**
       * \brief Time dependent divergence free velocity field
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * \author Jordi Paul
       */
      template<typename DT_, int dim_>
      class SinYT0 : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<dim_> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = dim_;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// The scaling factor according to the amplitude
          const DataType _t;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] function
           * The x,y plane velocity field function
           */
          explicit Evaluator(const SinYT0& function) :
            _t(function._t)

          {
            XASSERT(_t >= DataType(0));
          }

          /**
           * \brief Computes the value
           *
           * \param[out] val
           * The (vector valued) function value
           *
           * \param[in] point
           * The domain point
           */
          void value(ValueType& val, const PointType& point)
          {
            val.format();
            val(0) = Math::sin(_t*point[1]);
          }

          /**
           * \brief Computes the gradient
           *
           * \param[out] grad
           * The (matrix valued) gradient
           *
           * \param[in] point
           * The domain point
           */
          void gradient(GradientType& grad, const PointType& point) const
          {
            grad.format();
            grad[0][1] = _t*Math::cos(_t*point[1]);
          }

          /**
           * \brief Computes the Hessian
           *
           * \param[out] hess
           * The (tensor valued) Hessian
           *
           * \param[in] point
           * The domain point
           */
          void hessian(HessianType& hess, const PointType& DOXY(point)) const
          {
            hess.format();
          }

        }; // class SinYT0::Evaluator<...>

      private:
        /// The time
        DataType _t;

      public:
        /**
         * \brief Constructor
         *
         */
        explicit SinYT0(DataType t = DataType(0)) :
          _t(t)
        {
        }

        void set_time(const DataType t)
        {
          _t = t;
        }

      };

      /**
       * \brief Time dependent divergence free velocity field
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * \author Jordi Paul
       */
      template<typename DT_, int dim_>
      class SinYT0StokesRhs : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<dim_> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = dim_;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = false;
        /// We can compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          const DataType _t;
          const DataType _fac;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] function
           * The x,y plane velocity field function
           */
          explicit Evaluator(const SinYT0StokesRhs& function) :
            _t(function._t),
            _fac(DataType(1)/function._reynolds)

          {
            XASSERT(_t >= DataType(0));
          }

          /**
           * \brief Computes the value
           *
           * \param[out] val
           * The (vector valued) function value
           *
           * \param[in] point
           * The domain point
           */
          void value(ValueType& val, const PointType& point)
          {
            val.format();
            val(0) = point[1]*Math::cos(_t*point[1]) + _fac*Math::sqr(_t)*Math::sin(point[1]*_t);
          }

          /**
           * \brief Computes the gradient
           *
           * \param[out] grad
           * The (matrix valued) gradient
           *
           * \param[in] point
           * The domain point
           */
          void gradient(GradientType& grad, const PointType& DOXY(point)) const
          {
            grad.format();
          }

          /**
           * \brief Computes the Hessian
           *
           * \param[out] hess
           * The (tensor valued) Hessian
           *
           * \param[in] point
           * The domain point
           */
          void hessian(HessianType& hess, const PointType& DOXY(point)) const
          {
            hess.format();
          }

        }; // class SinYT0StokesRhs::Evaluator<...>

      private:
        /// The Reynolds number of the associated flow with the solution SinYT0
        const DataType _reynolds;
        /// The time
        DataType _t;

      public:
        /**
         * \brief Constructor
         *
         */
        explicit SinYT0StokesRhs(DataType reynolds, DataType t = DataType(0)) :
          _reynolds(reynolds),
          _t(t)
        {
          XASSERT(reynolds > DataType(0));
        }

        void set_time(const DataType t)
        {
          _t = t;
        }
      };

      /**
       * \brief Pressure to GuermondStokesSol
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * \author Jordi Paul
       */
      template<typename DT_, int dim_>
      class GuermondStokesSolPressure : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Scalar ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = dim_;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// The scaling factor according to the amplitude
          const DataType _t;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] function
           * The x,y plane velocity field function
           */
          explicit Evaluator(const GuermondStokesSolPressure& function) :
            _t(function._t)

          {
            XASSERT(_t >= DataType(0));
          }

          /**
           * \brief Computes the value
           *
           * \param[out] val
           * The (vector valued) function value
           *
           * \param[in] point
           * The domain point
           */
          void value(ValueType& val, const PointType& point)
          {
            val = Math::sin(point[0] - point[1] + _t);
            // To make it mean value 0 on [0,1]x[0,1]
            val -= (DataType(2)*Math::sin(_t) - Math::sin(DataType(1)+_t) - Math::sin(-DataType(1) + _t));
          }

          /**
           * \brief Computes the gradient
           *
           * \param[out] grad
           * The (matrix valued) gradient
           *
           * \param[in] point
           * The domain point
           */
          void gradient(GradientType& grad, const PointType& point) const
          {
            grad.format();

            grad[0] = Math::cos(point[0] - point[1] + _t);
            grad[1] = -grad[0];
          }

          /**
           * \brief Computes the Hessian
           *
           * \param[out] hess
           * The (tensor valued) Hessian
           *
           * \param[in] point
           * The domain point
           */
          void hessian(HessianType& hess, const PointType& DOXY(point)) const
          {
            hess.format();
          }

        }; // class GuermondStokesSolPressure::Evaluator<...>

      private:
        /// The time
        DataType _t;

      public:
        /**
         * \brief Constructor
         *
         */
        explicit GuermondStokesSolPressure(DataType t = DataType(0)) :
          _t(t)
        {
        }

        void set_time(const DataType t)
        {
          _t = t;
        }

      };

      /**
       * \brief Time dependent divergence free velocity field
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * \author Jordi Paul
       */
      template<typename DT_, int dim_>
      class GuermondStokesSol : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<dim_> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = dim_;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          /// The scaling factor according to the amplitude
          const DataType _t;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] function
           * The x,y plane velocity field function
           */
          explicit Evaluator(const GuermondStokesSol& function) :
            _t(function._t)

          {
            XASSERT(_t >= DataType(0));
          }

          /**
           * \brief Computes the value
           *
           * \param[out] val
           * The (vector valued) function value
           *
           * \param[in] point
           * The domain point
           */
          void value(ValueType& val, const PointType& point)
          {
            val.format();
            val(0) = Math::sin(point[0] + _t)*Math::sin(point[1] + _t);
            val(1) = Math::cos(point[0] + _t)*Math::cos(point[1] + _t);
          }

          /**
           * \brief Computes the gradient
           *
           * \param[out] grad
           * The (matrix valued) gradient
           *
           * \param[in] point
           * The domain point
           */
          void gradient(GradientType& grad, const PointType& point) const
          {
            grad.format();

            grad[0][0] = Math::cos(point[0] + _t)*Math::sin(point[1] + _t);
            grad[0][1] = Math::sin(point[0] + _t)*Math::cos(point[1] + _t);

            grad[1][0] = -Math::sin(point[0] + _t)*Math::cos(point[1] + _t);
            grad[1][1] = -Math::cos(point[0] + _t)*Math::sin(point[1] + _t);
          }

          /**
           * \brief Computes the Hessian
           *
           * \param[out] hess
           * The (tensor valued) Hessian
           *
           * \param[in] point
           * The domain point
           */
          void hessian(HessianType& hess, const PointType& DOXY(point)) const
          {
            hess.format();
          }

        }; // class GuermondStokesSol::Evaluator<...>

      private:
        /// The time
        DataType _t;

      public:
        /**
         * \brief Constructor
         *
         */
        explicit GuermondStokesSol(DataType t = DataType(0)) :
          _t(t)
        {
        }

        void set_time(const DataType t)
        {
          _t = t;
        }

      };

      /**
       * \brief Time dependent divergence free velocity field
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * \author Jordi Paul
       */
      template<typename DT_, int dim_>
      class GuermondStokesSolRhs : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<dim_> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = dim_;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = false;
        /// We can compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /** \copydoc AnalyticFunction::Evaluator */
        template<typename EvalTraits_>
        class Evaluator :
          public Analytic::Function::Evaluator<EvalTraits_>
        {
        public:
          /// coefficient data type
          typedef typename EvalTraits_::DataType DataType;
          /// evaluation point type
          typedef typename EvalTraits_::PointType PointType;
          /// value type
          typedef typename EvalTraits_::ValueType ValueType;
          /// gradient type
          typedef typename EvalTraits_::GradientType GradientType;
          /// hessian type
          typedef typename EvalTraits_::HessianType HessianType;

        private:
          const DataType _t;
          const DataType _fac;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] function
           * The x,y plane velocity field function
           */
          explicit Evaluator(const GuermondStokesSolRhs& function) :
            _t(function._t),
            _fac(DataType(1)/function._reynolds)

          {
            XASSERT(_t >= DataType(0));
          }

          /**
           * \brief Computes the value
           *
           * \param[out] val
           * The (vector valued) function value
           *
           * \param[in] point
           * The domain point
           */
          void value(ValueType& val, const PointType& point)
          {
            val.format();

            // Stationary
            //val(0) = _fac*DataType(2)*Math::sin(point[0]+_t)*Math::sin(point[1]+_t) + Math::cos(point[0]-point[1]+_t);
            //val(1) = _fac*DataType(2)*Math::cos(point[0]+_t)*Math::cos(point[1]+_t) - Math::cos(point[0]-point[1]+_t);

            val(0) =  Math::cos(point[0]+_t)*Math::sin(point[1]+_t) + Math::sin(point[0]+_t)*Math::cos(point[1]+_t)
              + _fac*DataType(2)*Math::sin(point[0]+_t)*Math::sin(point[1]+_t) + Math::cos(point[0]-point[1]+_t);
            val(1) = -Math::sin(point[0]+_t)*Math::cos(point[1]+_t) - Math::cos(point[0]+_t)*Math::sin(point[1]+_t)
              + _fac*DataType(2)*Math::cos(point[0]+_t)*Math::cos(point[1]+_t) - Math::cos(point[0]-point[1]+_t);
          }

          /**
           * \brief Computes the gradient
           *
           * \param[out] grad
           * The (matrix valued) gradient
           *
           * \param[in] point
           * The domain point
           */
          void gradient(GradientType& grad, const PointType& DOXY(point)) const
          {
            grad.format();
          }

          /**
           * \brief Computes the Hessian
           *
           * \param[out] hess
           * The (tensor valued) Hessian
           *
           * \param[in] point
           * The domain point
           */
          void hessian(HessianType& hess, const PointType& DOXY(point)) const
          {
            hess.format();
          }

        }; // class GuermondStokesSolRhs::Evaluator<...>

      private:
        /// The Reynolds number of the associated flow with the solution GuermondStokesSol
        const DataType _reynolds;
        /// The time
        DataType _t;

      public:
        /**
         * \brief Constructor
         *
         */
        explicit GuermondStokesSolRhs(DataType reynolds, DataType t = DataType(0)) :
          _reynolds(reynolds),
          _t(t)
        {
          XASSERT(reynolds > DataType(0));
        }

        void set_time(const DataType t)
        {
          _t = t;
        }
      };

    } // namespace Common
  } // namespace Analytic
} // namespace FEAT

#endif // KERNEL_ANALYTIC_COMMON_HPP
