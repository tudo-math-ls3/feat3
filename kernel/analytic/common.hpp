#pragma once
#ifndef KERNEL_ANALYTIC_COMMON_HPP
#define KERNEL_ANALYTIC_COMMON_HPP 1

// includes, FEAST
#include <kernel/analytic/static_wrapper.hpp>
#include <kernel/util/math.hpp>

// includes, system
#include <initializer_list>
#include <vector>

namespace FEAST
{
  namespace Analytic
  {
    /**
     * \brief Assembly Common namespace
     *
     * This namespace encapsulated commonly used functions, operators and functionals for use
     * with the various assembly classes, which are often used in standard benchmark problems.
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
        explicit DistanceFunction(const PointType origin) :
          _origin(origin)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType origin)
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
        explicit PlaneDistanceFunctionSD(const PointType origin, const DataType_ b_) :
          _origin(origin),
          _b(b_)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType origin)
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

      /**
       * \brief Heaviside static function
       *
       * This class implements the StaticFunction interface representing the function
       *                           ( 0, x  < 0
       *   - 1D: u(x)     = H(x) = (
       *                           ( 1, x >= 0
       *   - 2D: u(x,y)   = H(x)*H(y)
       *   - 3D: u(x,y,z) = H(x)*H(y)*H(z)
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

      /**
       * \brief Heaviside function
       *
       * This class implements the StaticFunction interface representing the function
       *                           ( 0, x  < 0
       *   - 1D: u(x)     = H(x) = (
       *                           ( 1, x >= 0
       *   - 2D: u(x,y)   = H(x)*H(y)
       *   - 3D: u(x,y,z) = H(x)*H(y)*H(z)
       *
       *   This class supports only values.
       *
       * \author Jordi Paul
       */
      template<int dim_>
      using HeavisideFunction = StaticWrapperFunction<dim_, HeavisideStatic, true, false, false>;

      /**
       * \brief Regularised Heaviside static function
       *
       * This class implements the StaticFunction interface representing the function
       *                           ( 0,               x  < 0
       *   - 1D: u(x)     = H(x) = (
       *                           ( 2 (cosh(x) - 1), x >= 0
       *   - 2D: u(x,y)   = H(x)*H(y)
       *   - 3D: u(x,y,z) = H(x)*H(y)*H(z)
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

      /**
       * \brief Regularised Heaviside function
       *
       * This class implements the StaticFunction interface representing the function
       *                           ( 0,               x  < 0
       *   - 1D: u(x)     = H(x) = (
       *                           ( 2 (cosh(x) - 1), x >= 0
       *   - 2D: u(x,y)   = H(x)*H(y)
       *   - 3D: u(x,y,z) = H(x)*H(y)*H(z)
       *
       *   This class supports function values, gradients and hessians for all dimensions.
       *
       * \author Jordi Paul
       */
      template<int dim_>
      using HeavisideRefFunction = StaticWrapperFunction<dim_, HeavisideRegStatic, true, false, false>;

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
    } // namespace Common
  } // namespace Analytic
} // namespace FEAST

#endif // KERNEL_ANALYTIC_COMMON_HPP
