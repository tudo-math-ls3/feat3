// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ANALYTIC_COMMON_HPP
#define KERNEL_ANALYTIC_COMMON_HPP 1

// includes, FEAT
#include <kernel/analytic/static_wrapper.hpp>
#include <kernel/analytic/wrappers.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/cgal.hpp>

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
        static DataType_ der_yx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xy(x, y, z);
        }

        /// 3D: XZ-derivative
        static DataType_ der_xz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(kpi()) * Math::sin(kpi() * x) * Math::cos(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: ZX-derivative
        static DataType_ der_zx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xz(x, y, z);
        }

        /// 3D: YZ-derivative
        static DataType_ der_yz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(kpi()) * Math::cos(kpi() * x) * Math::sin(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: ZY-derivative
        static DataType_ der_zy(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_yz(x, y, z);
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
       *   - u(x) = (exp(-(2*x - 1)^2) - exp(-1)) / (1 - exp(-1))
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
       *   - 1D: u(x)     = (exp(-(2*x - 1)^2) - exp(-1)) / (1 - exp(-1))
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
       * \brief Q2-Bubble Analytic function
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

          ValueType value(const PointType& DOXY(point))
          {
            return _value;
          }

          GradientType gradient(const PointType& DOXY(point))
          {
            return GradientType::null();
          }

          HessianType hessian(const PointType& DOXY(point))
          {
            return HessianType::null();
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
       * \brief Constant vector valued Analytic function
       *
       * This class implements the AnalyticFunction interface representing a constant vector valued function.
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       * \author Maximilian Esser
       */
      template<int dim_, typename DataType_ = Real>
      class ConstantVectorFunction :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Vector<domain_dim> ImageType;

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
          const ValueType _value_vec;

        public:
          /// Constructor
          explicit Evaluator(const ConstantVectorFunction& function) :
            _value_vec(function._value_vector)
          {
          }

          ValueType value(const PointType& DOXY(point))
          {
            return _value_vec;
          }

          GradientType gradient(const PointType& DOXY(point))
          {
            return GradientType::null();
          }

          HessianType hessian(const PointType& DOXY(point))
          {
            return HessianType::null();
          }
        }; // class ConstantVectorFunction::Evaluator<...>

      private:
        /// Tiny Vector representing the values of the constant function
        Tiny::Vector<DataType_, domain_dim> _value_vector;

      public:
        /// Constructor, value defaults to 0
        /// This constructor sets all values to the same given value
        explicit ConstantVectorFunction(DataType_ value = DataType_(0)) :
          _value_vector(Tiny::Vector<DataType_, domain_dim>(value))
        {
        }

        /// This constructor sets _value_vec to the given Tiny::Vector vec
        explicit ConstantVectorFunction(Tiny::Vector<DataType_, domain_dim> vec) :
          _value_vector(vec)
        {
        }
      }; // class ConstantVectorFunction

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

          const DataType _eps;

        public:
          /// Constructor
          explicit Evaluator(const MinOfTwoFunctions& function) :
            _f1_eval(function._f1),
            _f2_eval(function._f2),
            _eps(Math::eps<DataType>())
          {
          }

          ValueType value(const PointType& point)
          {
            return Math::min(_f1_eval.value(point), _f2_eval.value(point));
          }

          GradientType gradient(const PointType& point)
          {
            ValueType val1 = _f1_eval.value(point);
            ValueType val2 = _f2_eval.value(point);

            if(Math::abs(val1-val2) < _eps)
              return GradientType::null();
            else if(val1 < val2)
              return _f1_eval.gradient(point);
            else
              return _f2_eval.gradient(point);
          }

          HessianType hessian(const PointType& point)
          {
            ValueType val1 = _f1_eval.value(point);
            ValueType val2 = _f2_eval.value(point);

            if(Math::abs(val1-val2) < _eps)
              return HessianType::null();
            else if(val1 < val2)
              return _f1_eval.hessian(point);
            else
              return _f2_eval.hessian(point);
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
          return eval(x)*eval(y)*eval(z);
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
       * \brief Regularized Heaviside static function
       *
       * This class implements the StaticFunction interface representing the function
       * \f[
       *   u(x) = H(x) =
       *   \begin{cases}
       *     0, & x < 0 \\
       *     2 (\cosh(x) - 1), & x \geq 0
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
          return eval(x)*eval(y)*eval(z);
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
       * \brief Regularized Heaviside static function
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

          ValueType value(const PointType& point) const
          {
            // evaluate polynomial via horner scheme
            DataType x = point[0];
            DataType y = DataType(0);
            for(std::size_t k(_coeff.size()); k > std::size_t(0); )
              y = x * y + _coeff[--k];

            return y;
          }

          GradientType gradient(const PointType& point) const
          {
            std::size_t k = _coeff.size();
            if(k <= std::size_t(0))
              return GradientType::null();

            // evaluate polynomial via horner scheme
            DataType x = point[0];
            DataType y = DataType(0);
            for( ; (--k) > std::size_t(0); )
              y = x * y + (_coeff[k] * DataType(k));

            return GradientType(y);
          }

          HessianType hessian(const PointType& point) const
          {
            std::size_t k = _coeff.size();
            if(k <= std::size_t(1))
              return HessianType::null();

            // evaluate polynomial via horner scheme
            DataType x = point[0];
            DataType y = DataType(0);
            for( ; (--k) > std::size_t(1); )
              y = x * y + (_coeff[k] * DataType(k*(k-1)));

            return HessianType(y);
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
       *  testing optimization algorithms because the hessian at the minimal point is singular.
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
          return DataType_(-4);
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
       * \brief Goldstein-Price function
       *
       * \tparam DataType_
       * Floating point precision
       *
       * This class implements the StaticFunction interface representing the function
       *  \f[
       *    u(x,y) =(1+(1+x+y)^2 * (19 - 14*x+3*x^2-14*y+6*x*y+3*y^2))
       *    *(30 + (2*x-3*y)^2 * (18 - 32*x +12*x^2 + 48*y - 36*x*y + 27*y^2))
       *  \f]
       *  This function has a global miminum in \f$ x_0 = (0, -1)^T, u(x_0) = 3 \f$.
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
          return (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (DataType_(30) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)));
        }

        /// x-derivative
        static DataType_ der_x(DataType_ x, DataType_ y)
        {
          return (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (DataType_(30) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (DataType_(4) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(24) * x - DataType_(36) * y - DataType_(32)));
        }

        /// y-derivative
        static DataType_ der_y(DataType_ x, DataType_ y)
        {
          return (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (DataType_(30) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (-DataType_(6) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48)));
        }

        /// xx-derivative
        static DataType_ der_xx(DataType_ x, DataType_ y)
        {
          return DataType_((DataType_(6) * Math::pow(x, DataType_(2)) + DataType_(12) * x * y + DataType_(6) * Math::pow(y, DataType_(2)) - DataType_(28) * x - DataType_(28) * y + DataType_(38) + DataType_(4) * (DataType_(1) + x + y) * (DataType_(6) * x + DataType_(6) * y - DataType_(14)) + DataType_(6) * Math::pow(DataType_(1) + x + y, DataType_(2))) * (DataType_(30) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + DataType_(2) * (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (DataType_(4) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(24) * x - DataType_(36) * y - DataType_(32))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (DataType_(96) * Math::pow(x, DataType_(2)) - DataType_(288) * x * y + DataType_(216) * Math::pow(y, DataType_(2)) - DataType_(256) * x + DataType_(384) * y + DataType_(144) + DataType_(8) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(24) * x - DataType_(36) * y - DataType_(32)) + DataType_(24) * Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2))));
        }

        /// yy-derivative
        static DataType_ der_yy(DataType_ x, DataType_ y)
        {
          return (DataType_(6) * Math::pow(x, DataType_(2)) + DataType_(12) * x * y + DataType_(6) * Math::pow(y, DataType_(2)) - DataType_(28) * x - DataType_(28) * y + DataType_(38) + DataType_(4) * (DataType_(1) + x + y) * (DataType_(6) * x + DataType_(6) * y - DataType_(14)) + DataType_(6) * Math::pow(DataType_(1) + x + y, DataType_(2))) * (DataType_(30) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + DataType_(2) * (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (-DataType_(6) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (DataType_(216) * Math::pow(x, DataType_(2)) - DataType_(648) * x * y + DataType_(486) * Math::pow(y, DataType_(2)) - DataType_(576) * x + DataType_(864) * y + DataType_(324) - DataType_(12) * (DataType_(2) * x - DataType_(3) * y) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48)) + DataType_(54) * Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)));
        }

        /// xy-derivative
        static DataType_ der_xy(DataType_ x, DataType_ y)
        {
          return (DataType_(6) * Math::pow(x, DataType_(2)) + DataType_(12) * x * y + DataType_(6) * Math::pow(y, DataType_(2)) - DataType_(28) * x - DataType_(28) * y + DataType_(38) + DataType_(4) * (DataType_(1) + x + y) * (DataType_(6) * x + DataType_(6) * y - DataType_(14)) + DataType_(6) * Math::pow(DataType_(1) + x + y, DataType_(2))) * (DataType_(30) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18))) + (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (DataType_(4) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (DataType_(24) * x - DataType_(36) * y - DataType_(32))) + (DataType_(2) * (DataType_(1) + x + y) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19)) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(6) * x + DataType_(6) * y - DataType_(14))) * (-DataType_(6) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(12) * Math::pow(x, DataType_(2)) - DataType_(36) * x * y + DataType_(27) * Math::pow(y, DataType_(2)) - DataType_(32) * x + DataType_(48) * y + DataType_(18)) + Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48))) + (DataType_(1) + Math::pow(DataType_(1) + x + y, DataType_(2)) * (DataType_(3) * Math::pow(x, DataType_(2)) + DataType_(6) * x * y + DataType_(3) * Math::pow(y, DataType_(2)) - DataType_(14) * x - DataType_(14) * y + DataType_(19))) * (-DataType_(144) * Math::pow(x, DataType_(2)) + DataType_(432) * x * y - DataType_(324) * Math::pow(y, DataType_(2)) + DataType_(384) * x - DataType_(576) * y - DataType_(216) - DataType_(6) * (DataType_(2) * x - DataType_(3) * y) * (DataType_(24) * x - DataType_(36) * y - DataType_(32)) + DataType_(4) * (DataType_(2) * x - DataType_(3) * y) * (-DataType_(36) * x + DataType_(54) * y + DataType_(48)) - DataType_(36) * Math::pow(DataType_(2) * x - DataType_(3) * y, DataType_(2)));
        }

        /// xy-derivative
        static DataType_ der_yx(DataType_ x, DataType_ y)
        {
          return der_xy(x, y);
        }

      }; // class GoldsteinPriceStatic<...>
      /// \endcond

      /**
       * \brief Goldstein-Price function
       *
       * \tparam DataType_
       * Floating point precision
       *
       * This class implements the StaticFunction interface representing the function
       *  \f[
       *    u(x,y) =(1+(1+x+y)^2 * (19 - 14*x+3*x^2-14*y+6*x*y+3*y^2))
       *    *(30 + (2*x-3*y)^2 * (18 - 32*x +12*x^2 + 48*y - 36*x*y + 27*y^2))
       *  \f]
       *  This function has a global miminum in \f$ x_0 = (0, -1)^T, u(x_0) = 3 \f$.
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
       *  testing optimization algorithms because the hessian at the minimal point is singular.
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
       *    u(x,y) = (x^2 + y - 11)^2 + (x + y^2 - 7)^2
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
       *  It is often used for testing optimization algorithms because of the nonconvexity and existence of a local
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
       *  It is often used for testing optimization algorithms because of the nonconvexity and existence of a local
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
       * \f$ y = x^2 \f$. This is a great challenge to descend-based optimization algorithms like steepest descent or
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
       * parameter of the orthogonal projection of (x,y) onto the line segment, then the
       * parabolic profile function value is given by \f$4 v_{max} t (t-1)\f$.
       *
       * \author Peter Zajac
       */
      template<typename DataType_>
      class ParProfileBase :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;

      protected:
        // coordinates of line segment
        DataType_ _x0, _y0, _x1, _y1;
        // maximum value
        DataType_ _vmax;

      public:
        /**
         * \brief Default Constructor
         *
         * This constructor initializes a parabolic profile along the segment (0,0)-(0,1)
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
        explicit ParProfileBase(DataType_ x0, DataType_ y0, DataType_ x1, DataType_ y1, DataType_ vmax = DataType_(1.0)) :
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

          std::deque<String> sv = sbc.substr(li+1, ri-li-1).split_by_string(",");
          if((sv.size() < std::size_t(2)) || (sv.size() > std::size_t(3)))
            return false;

          std::deque<String> sv0 = sv[0].trim().split_by_whitespaces();
          std::deque<String> sv1 = sv[1].trim().split_by_whitespaces();
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
      template<typename DataType_>
      class ParProfileScalar :
        public ParProfileBase<DataType_>
      {
      public:
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        using ParProfileBase<DataType_>::ParProfileBase;

        template<typename Traits_>
        class Evaluator :
          public Analytic::Function::Evaluator<Traits_>
        {
        protected:
          typedef typename Traits_::DataType DataType;
          typedef typename Traits_::PointType PointType;
          typedef typename Traits_::ValueType ValueType;
          typedef typename Traits_::GradientType GradientType;
          typedef typename Traits_::HessianType HessianType;

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

          ValueType value(const PointType& point)
          {
            // project point onto line segment
            const DataType x = Math::clamp(Tiny::dot(point - _vo, _ve) * _den, DataType(0), DataType(1));
            // compute function value
            return _vmax * DataType(4) * x * (DataType(1) - x);
          }

          GradientType gradient(const PointType& point)
          {
            // project point onto line segment
            const DataType x = Tiny::dot(point - _vo, _ve) * _den;

            // Note: the gradient has a singularity at x=0 and x=1
            if((x < DataType(0)) || (x > DataType(1)))
              return GradientType::null();

            return (_vmax * _den * DataType(4) * (DataType(1) - DataType(2)*x)) * _ve;
          }

          HessianType hessian(const PointType& point)
          {
            // project point onto line segment
            const DataType x = Tiny::dot(point - _vo, _ve) * _den;

            // Note: the hessian has a singularity at x=0 and x=1
            if((x < DataType(0)) || (x > DataType(1)))
              return HessianType::null();

            HessianType hess;
            const DataType v = -DataType(8) * _vmax * _den * _den;
            hess[0][0] = v * _ve[0] * _ve[0];
            hess[0][1] = hess[1][0] = v * _ve[0] * _ve[1];
            hess[1][1] = v * _ve[1] * _ve[1];
            return hess;
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
      template<typename DataType_>
      class ParProfileVector :
        public ParProfileBase<DataType_>
      {
      public:
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Vector<2> ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        using ParProfileBase<DataType_>::ParProfileBase;

        template<typename Traits_>
        class Evaluator :
          public Analytic::Function::Evaluator<Traits_>
        {
        protected:
          typedef typename Traits_::DataType DataType;
          typedef typename Traits_::PointType PointType;
          typedef typename Traits_::ValueType ValueType;
          typedef typename Traits_::GradientType GradientType;
          typedef typename Traits_::HessianType HessianType;

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
            _vn.normalize();
            _den = DataType(1) / Tiny::dot(_ve, _ve);
            _vmax = DataType(function._vmax);
          }

          ValueType value(const PointType& point)
          {
            // project point onto line segment
            const DataType x = Math::clamp(Tiny::dot(point - _vo, _ve) * _den, DataType(0), DataType(1));
            // compute function value
            const DataType v = _vmax * DataType(4) * x * (DataType(1) - x);
            ValueType val;
            val[0] = _vn[0] * v;
            val[1] = _vn[1] * v;
            return val;
          }

          GradientType gradient(const PointType& point)
          {
            // project point onto line segment
            const DataType x = Tiny::dot(point - _vo, _ve) * _den;

            // Note: the gradient has a singularity at x=0 and x=1
            if((x < DataType(0)) || (x > DataType(1)))
              return GradientType::null();

            const DataType v = _vmax * _den * DataType(4) * (DataType(1) - DataType(2)*x);
            GradientType grad;
            grad(0,0) =  v * _vn[0] * _ve[0];
            grad(0,1) =  v * _vn[0] * _ve[1];
            grad(1,0) =  v * _vn[1] * _ve[0];
            grad(1,1) =  v * _vn[1] * _ve[1];
            return grad;
          }

          HessianType hessian(const PointType& point)
          {
            // project point onto line segment
            const DataType x = Tiny::dot(point - _vo, _ve) * _den;

            // Note: the hessian has a singularity at x=0 and x=1
            if((x < DataType(0)) || (x > DataType(1)))
              return HessianType::null();

            const DataType v = -DataType(8) * _vmax * _den * _den;
            HessianType hess;
            hess(0,0,0) = v * _vn[0] * _ve[0] * _ve[0];
            hess(0,1,1) = v * _vn[0] * _ve[1] * _ve[1];
            hess(0,1,0) = hess(0,0,1) = v * _vn[0] * _ve[0] * _ve[1];
            hess(1,1,0) = hess(1,0,1) = v * _vn[1] * _ve[0] * _ve[1];
            hess(1,0,0) = v * _vn[1] * _ve[0] * _ve[0];
            hess(1,1,1) = v * _vn[1] * _ve[1] * _ve[1];
            return hess;
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

      #ifdef FEAT_HAVE_HALFMATH
      template<>
      class ExpScalarStatic<Half>
      {
      public:
        typedef Half DT_;
        static DT_ eval(DT_ x)
        {
          return (Math::exp(DT_(10)) - Math::exp(DT_(10)*x*x)) / (Math::exp(DT_(10)) - DT_(1));
        }

        static DT_ der_x(DT_ x)
        {
          return -DT_(2)*DT_(10)*x*Math::exp(DT_(10)*x*x)/(Math::exp(DT_(10))-DT_(1));
        }

        static DT_ der_xx(DT_ x)
        {
          return -DT_(2)*DT_(10)*Math::exp(DT_(10)*x*x)*(DT_(2)*DT_(10)*x*x+DT_(1))/(Math::exp(DT_(10))-DT_(1));
        }
      };
      #endif

      /// \cond internal
      template<typename DataType_>
      using ExpStatic = Analytic::Common::TensorStatic<ExpScalarStatic<DataType_>, DataType_>;
      /// \endcond

      /**
       * \brief Exponential Analytic function
       *
       * This class implements the AnalyticFunction interface representing the function
       *   - 1D: u(x)     =  (exp(10) - exp(10*x^2)) / (exp(10) - 1)
       *   - 2D: u(x,y)   =  u(x)*u(y)
       *   - 3D: u(x,y,z) =  u(x)*u(y)*u(z)
       *
       * This class supports function values, gradient and hessians for all dimensions.
       *
       */
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
          explicit Evaluator(const XYPlaneRotation& function) :
            _angular_velocity(function._angular_velocity),
            _origin(function._origin)
          {
          }

          ValueType value(const PointType& point)
          {
            ValueType val;
            val.format();
            val[0] = _angular_velocity*(-(point[1] - _origin[1]));
            val[1] = _angular_velocity*( (point[0] - _origin[0]));
            return val;
          }

          GradientType gradient(const PointType& DOXY(point)) const
          {
            GradientType grad;
            grad.format();
            grad[0][1] = -_angular_velocity;
            grad[1][0] = _angular_velocity;
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
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
       *    f: \mathbb{R}^d \to \mathbb{R}^d, f_0(x_0, ..., x_d)^T =
       *    \alpha \frac{4^{d-1}}{\prod_{i=1}^{d-1}(b_i - a_i)^2}\prod_{i=1}^{d-1}( (x_i - a_i)(b_i - x_i) ),
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
        typedef Tiny::Vector<DT_, 2> RangeType;

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
          const std::vector<Tiny::Vector<DataType, 2>>& _range;

        public:
          explicit Evaluator(const YZPlaneParabolic& function) :
            _fac(function._amplitude),
            _range(function._range)

          {
            XASSERT(_range.size() == size_t(domain_dim-1));

            for (int d(1); d < domain_dim; ++d)
            {
              _fac *= DataType(4)/Math::sqr(_range.at(std::size_t(d-1))[1]-_range.at(std::size_t(d-1))[0]);
            }
          }

          ValueType value(const PointType& point)
          {
            ValueType val;
            val.format();
            val(0) = _fac;

            for (int d(1); d < domain_dim; ++d)
            {
              val(0) *= (point[d] - _range.at(std::size_t(d-1))[0])*(_range.at(std::size_t(d-1))[1] - point[d]);
            }

            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad;
            grad.format();
            for(int d(1); d < domain_dim; ++d)
            {
              grad[d][0] = _fac*(_range.at(std::size_t(d-1))[0] + _range.at(std::size_t(d-1))[1] - DataType(2)*point(d));

              for (int q(1); q < domain_dim; ++q)
              {
                if (q != d)
                  grad[d][0] *= (point[q] - _range.at(std::size_t(q - 1))[0]) * (_range.at(std::size_t(q - 1))[1] - point[q]);
              }
            }
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
          }

        }; // class YZPlaneParabolic::Evaluator<...>

      public:
        /// The maximum value of the parabolic profile
        const DataType _amplitude;
        /// The points where the function becomes zero
        std::vector<Tiny::Vector<DataType, 2>> _range;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] amplitude
         * The amplitude to use
         *
         * \param[in] range_y
         * The roots for the y part
         *
         */
        explicit YZPlaneParabolic(const DataType amplitude, const RangeType& range_y) :
          _amplitude(amplitude),
          _range(1)
        {
          _range.at(0) = range_y;
        }

        /**
         * \brief Constructor
         *
         * \param[in] amplitude
         * The amplitude to use
         *
         * \param[in] range_y
         * The roots for the y part
         *
         * \param[in] range_z
         * The roots for the z part
         *
         */
        explicit YZPlaneParabolic(const DataType amplitude, const RangeType& range_y, const RangeType& range_z) :
          _amplitude(amplitude),
          _range(2)
        {
          _range.at(0) = range_y;
          _range.at(1) = range_z;
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
          explicit Evaluator(const SinYT0& function) :
            _t(function._t)

          {
            XASSERT(_t >= DataType(0));
          }

          ValueType value(const PointType& point)
          {
            ValueType val;
            val.format();
            val(0) = Math::sin(_t*point[1]);
            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad;
            grad.format();
            grad[0][1] = _t*Math::cos(_t*point[1]);
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
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
          explicit Evaluator(const SinYT0StokesRhs& function) :
            _t(function._t),
            _fac(DataType(1)/function._reynolds)

          {
            XASSERT(_t >= DataType(0));
          }

          ValueType value(const PointType& point)
          {
            ValueType val;
            val.format();
            val(0) = point[1]*Math::cos(_t*point[1]) + _fac*Math::sqr(_t)*Math::sin(point[1]*_t);
            return val;
          }

          GradientType gradient(const PointType& DOXY(point)) const
          {
            return GradientType::null();
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
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
          explicit Evaluator(const GuermondStokesSolPressure& function) :
            _t(function._t)

          {
            XASSERT(_t >= DataType(0));
          }

          ValueType value( const PointType& point)
          {
            ValueType val;
            val = Math::sin(point[0] - point[1] + _t);
            // To make it mean value 0 on [0,1]x[0,1]
            val -= (DataType(2)*Math::sin(_t) - Math::sin(DataType(1)+_t) - Math::sin(-DataType(1) + _t));
            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad;
            grad.format();
            grad[0] = Math::cos(point[0] - point[1] + _t);
            grad[1] = -grad[0];
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
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
          explicit Evaluator(const GuermondStokesSol& function) :
            _t(function._t)
          {
            XASSERT(_t >= DataType(0));
          }

          ValueType value(const PointType& point)
          {
            ValueType val;
            val.format();
            val(0) = Math::sin(point[0] + _t)*Math::sin(point[1] + _t);
            val(1) = Math::cos(point[0] + _t)*Math::cos(point[1] + _t);
            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad;
            grad.format();

            grad[0][0] = Math::cos(point[0] + _t)*Math::sin(point[1] + _t);
            grad[0][1] = Math::sin(point[0] + _t)*Math::cos(point[1] + _t);

            grad[1][0] = -Math::sin(point[0] + _t)*Math::cos(point[1] + _t);
            grad[1][1] = -Math::cos(point[0] + _t)*Math::sin(point[1] + _t);
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
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
          explicit Evaluator(const GuermondStokesSolRhs& function) :
            _t(function._t),
            _fac(DataType(1)/function._reynolds)

          {
            XASSERT(_t >= DataType(0));
          }

          ValueType value(const PointType& point)
          {
            ValueType val;
            val.format();

            // Stationary
            //val(0) = _fac*DataType(2)*Math::sin(point[0]+_t)*Math::sin(point[1]+_t) + Math::cos(point[0]-point[1]+_t);
            //val(1) = _fac*DataType(2)*Math::cos(point[0]+_t)*Math::cos(point[1]+_t) - Math::cos(point[0]-point[1]+_t);

            val(0) =  Math::cos(point[0]+_t)*Math::sin(point[1]+_t) + Math::sin(point[0]+_t)*Math::cos(point[1]+_t)
              + _fac*DataType(2)*Math::sin(point[0]+_t)*Math::sin(point[1]+_t) + Math::cos(point[0]-point[1]+_t);
            val(1) = -Math::sin(point[0]+_t)*Math::cos(point[1]+_t) - Math::cos(point[0]+_t)*Math::sin(point[1]+_t)
              + _fac*DataType(2)*Math::cos(point[0]+_t)*Math::cos(point[1]+_t) - Math::cos(point[0]-point[1]+_t);
            return val;
          }

          GradientType gradient(const PointType& DOXY(point)) const
          {
            return GradientType::null();
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
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


      /**
       * \brief Sphere-Normalized Sine-Bubble Function
       *
       * This function represents the 3D sine-bubble function, whose input
       * point has been normalized, i.e.:
       *
       * \f[u(x,y,z) = \sin\bigg(\frac{\pi x}{\sqrt{x^2+y^2+z^2}}\bigg)\cdot
          \sin\bigg(\frac{\pi y}{\sqrt{x^2+y^2+z^2}}\bigg)\cdot
          \sin\bigg(\frac{\pi z}{\sqrt{x^2+y^2+z^2}}\bigg)\f]
       *
       * \author Peter Zajac
       */
      class SphereSinBubbleFunction : public Analytic::Function
      {
      public:
        /// What type this mapping maps to
        typedef Analytic::Image::Scalar ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 3;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can compute the Hessian
        static constexpr bool can_hess = true;

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
          const DataType pi;

        public:
          explicit Evaluator(const SphereSinBubbleFunction&) :
            pi(Math::pi<DataType>())
          {
          }

          ValueType value(const PointType& point)
          {
            const DataType re = pi / point.norm_euclid();
            return Math::sin(re*point[0]) * Math::sin(re*point[1]) * Math::sin(re*point[2]);
          }

          GradientType gradient(const PointType& point) const
          {
            const DataType x = point[0];
            const DataType y = point[1];
            const DataType z = point[2];
            const DataType re = DataType(1) / point.norm_euclid();
            const DataType cx = Math::cos(pi * re * x);
            const DataType cy = Math::cos(pi * re * y);
            const DataType cz = Math::cos(pi * re * z);
            const DataType sx = Math::sin(pi * re * x);
            const DataType sy = Math::sin(pi * re * y);
            const DataType sz = Math::sin(pi * re * z);
            GradientType grad;
            grad[0] = +pi * Math::cub(re) * (cx*sz*sy*y*y + cx*sz*sy*z*z - sx*x*y*cy*sz - sx*sy*x*z*cz);
            grad[1] = -pi * Math::cub(re) * (x*y*cx*sy*sz - sz*sx*cy*x*x - sz*sx*cy*z*z + sx*sy*z*y*cz);
            grad[2] = -pi * Math::cub(re) * (x*z*cx*sy*sz + sx*z*y*cy*sz - sx*cz*sy*x*x - sx*cz*sy*y*y);
            return grad;
          }

          HessianType hessian(const PointType& point) const
          {
            const DataType x = point[0];
            const DataType y = point[1];
            const DataType z = point[2];
            const DataType di = Math::sqrt(x*x + y*y + z*z);
            const DataType re = DataType(1) / di;
            const DataType qo = DataType(1) / Math::cub(x*x + y*y + z*z);
            const DataType cx = Math::cos(pi * re * x);
            const DataType cy = Math::cos(pi * re * y);
            const DataType cz = Math::cos(pi * re * z);
            const DataType sx = Math::sin(pi * re * x);
            const DataType sy = Math::sin(pi * re * y);
            const DataType sz = Math::sin(pi * re * z);
            HessianType hess;
            hess[0][0] = -pi*qo*(
              DataType(2)*cx*sy*cz*pi*x*y*y*z +
              DataType(2)*cx*sy*cz*pi*x*z*z*z +
              DataType(2)*cx*sz*cy*pi*x*y*y*y +
              DataType(2)*cx*sz*cy*pi*x*y*z*z +
              DataType(2)*sy*sz*sx*pi*y*y*z*z +
              -DataType(2)*sx*y*x*x*cy*sz*di +
              -DataType(2)*sx*y*pi*x*x*cy*z*cz +
              -DataType(2)*sx*sy*z*x*x*cz*di +
              DataType(3)*cx*sy*sz*x*di*y*y +
              DataType(3)*cx*sy*sz*x*di*z*z +
              sx*y*y*pi*x*x*sy*sz +
              sx*sy*z*z*pi*x*x*sz +
              sy*sz*sx*pi*y*y*y*y +
              sy*sz*sx*pi*z*z*z*z +
              sy*cz*sx*di*y*y*z +
              sy*cz*sx*di*z*z*z +
              sz*cy*sx*di*y*y*y +
              sz*cy*sx*di*y*z*z);
            hess[1][1] = -pi*qo*(
              DataType(2)*cx*sz*cy*pi*x*x*x*y +
              DataType(2)*cx*sz*cy*pi*x*y*z*z +
              DataType(2)*sx*sy*z*z*pi*x*x*sz +
              DataType(2)*sx*y*pi*x*x*cy*z*cz +
              DataType(2)*sx*cy*cz*pi*y*z*z*z +
              -DataType(2)*cx*sy*cz*pi*x*y*y*z +
              -DataType(2)*cx*sy*sz*x*di*y*y +
              -DataType(2)*sy*cz*sx*di*y*y*z +
              DataType(3)*sx*y*x*x*cy*sz*di +
              DataType(3)*sz*cy*sx*di*y*z*z +
              sy*sz*sx*pi*x*x*x*x +
              sx*y*y*pi*x*x*sy*sz +
              sy*sz*sx*pi*y*y*z*z +
              sy*sz*sx*pi*z*z*z*z +
              di*cx*sy*sz*x*x*x +
              cx*sy*sz*x*di*z*z +
              sx*sy*z*x*x*cz*di +
              sy*cz*sx*di*z*z*z);
            hess[2][2] = -pi*qo*(
              DataType(2)*cx*sy*cz*pi*x*x*x*z +
              DataType(2)*cx*sy*cz*pi*x*y*y*z +
              DataType(2)*sx*y*y*pi*x*x*sy*sz +
              DataType(2)*sx*y*pi*x*x*cy*z*cz +
              DataType(2)*sx*cy*cz*pi*y*y*y*z +
              -DataType(2)*cx*sy*sz*x*di*z*z +
              -DataType(2)*cx*sz*cy*pi*x*y*z*z +
              -DataType(2)*sz*cy*sx*di*y*z*z +
              DataType(3)*sx*sy*z*x*x*cz*di +
              DataType(3)*sy*cz*sx*di*y*y*z +
              sy*sz*sx*pi*x*x*x*x +
              sx*sy*z*z*pi*x*x*sz +
              sy*sz*sx*pi*y*y*y*y +
              sy*sz*sx*pi*y*y*z*z +
              di*cx*sy*sz*x*x*x +
              cx*sy*sz*x*di*y*y +
              sx*y*x*x*cy*sz*di +
              sz*cy*sx*di*y*y*y);
            hess[0][1] = hess[1][0] = pi*qo*(
              DataType(2)*y*y*pi*x*x*cx*cy*sz +
              DataType(2)*di*cx*sy*sz*x*x*y +
              DataType(2)*sx*y*y*x*cy*sz*di +
              DataType(3)*sx*sy*z*x*cz*y*di +
              y*pi*x*x*cx*sy*z*cz +
              -cx*sy*cz*pi*y*y*y*z +
              -cx*sy*cz*pi*y*z*z*z +
              cx*sz*cy*pi*x*x*z*z +
              cx*sz*cy*pi*y*y*z*z +
              cx*sz*cy*pi*z*z*z*z +
              sy*sz*sx*pi*x*x*x*y +
              sy*sz*sx*pi*x*y*y*y +
              sx*sy*z*z*pi*x*y*sz +
              -sx*cy*cz*pi*x*x*x*z +
              sx*y*y*pi*x*cy*z*cz +
              -sx*cy*cz*pi*x*z*z*z +
              -cx*sy*sz*di*y*y*y +
              -cx*sy*sz*di*y*z*z +
              -cy*sz*sx*x*x*x*di +
              -cy*sz*sx*x*di*z*z);
            hess[0][2] = hess[2][0] = pi*qo*(
              DataType(2)*z*z*pi*x*x*cx*sy*cz +
              DataType(2)*sx*sy*z*z*x*cz*di +
              DataType(2)*di*cx*sy*sz*x*x*z +
              DataType(3)*sx*y*x*cy*sz*z*di +
              cx*sy*cz*pi*x*x*y*y +
              cx*sy*cz*pi*y*y*y*y +
              cx*sy*cz*pi*y*y*z*z +
              sy*sz*sx*pi*x*x*x*z +
              sx*y*y*pi*x*z*sy*sz +
              sy*sz*sx*pi*x*z*z*z +
              -sx*cy*cz*pi*x*x*x*y +
              -sx*cy*cz*pi*x*y*y*y +
              sx*y*pi*z*z*cy*x*cz +
              z*pi*x*x*cx*y*cy*sz +
              -cx*sz*cy*pi*y*y*y*z +
              -cx*sz*cy*pi*y*z*z*z +
              -sy*cz*sx*x*x*x*di +
              -sy*cz*sx*x*di*y*y +
              -sy*sz*cx*di*y*y*z +
              -sy*sz*cx*di*z*z*z);
            hess[1][2] = hess[2][1] = pi*qo*(
              DataType(2)*sx*y*y*pi*z*z*cy*cz +
              DataType(2)*sx*sy*z*z*y*cz*di +
              DataType(2)*di*sz*sx*cy*y*y*z +
              DataType(3)*y*x*cx*sy*sz*z*di +
              -cx*sy*cz*pi*x*x*x*y +
              -cx*sy*cz*pi*x*y*y*y +
              z*z*pi*x*cx*sy*y*cz +
              -cx*sz*cy*pi*x*x*x*z +
              y*y*pi*x*cx*z*cy*sz +
              -cx*sz*cy*pi*x*z*z*z +
              y*pi*x*x*z*sx*sy*sz +
              sy*sz*sx*pi*y*y*y*z +
              sy*sz*sx*pi*y*z*z*z +
              sx*cy*cz*pi*x*x*x*x +
              sx*cy*cz*pi*x*x*y*y +
              sx*cy*cz*pi*x*x*z*z +
              -sy*cz*sx*x*x*di*y +
              -sy*cz*sx*di*y*y*y +
              -sz*cy*sx*x*x*di*z +
              -sz*cy*sx*di*z*z*z);
            return hess;
          }
        }; // class SphereSinBubbleFunction::Evaluator<...>
      }; // class SphereSinBubbleFunction

      /**
       * \brief Standing-Vortex function in 2D
       *
       * This function represents the 2D standing vortex velocity field, which is a common benchmark
       * velocity field for testing the stability of Navier-Stokes discretizations and stabilizations.
       *
       * The field is given with \f$r := \sqrt{(x-0.5)^2 \cdot(y-0.5)^2}\f$ as
       *
       * \f[u(x,y,z) = \begin{cases}\begin{bmatrix}5(y-0.5)\\ -5(x-0.5)\end{bmatrix},& r < 0.2\\
         \begin{bmatrix}(2/r - 5)(y-0.5) \\ (-2/r + 5)(x-0.5)\end{bmatrix},& 0.2 \leq r < 0.4\\0, & \text{else}\end{cases}\f]
       *
       * \author Peter Zajac
       */
      class StandingVortexFunction2D :
        public Analytic::Function
      {
      public:
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Vector<2> ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = false;

        template<typename Traits_>
        class Evaluator :
          public Analytic::Function::Evaluator<Traits_>
        {
        protected:
          typedef typename Traits_::DataType DataType;
          typedef typename Traits_::PointType PointType;
          typedef typename Traits_::ValueType ValueType;
          typedef typename Traits_::GradientType GradientType;

        public:
          explicit Evaluator(const StandingVortexFunction2D&)
          {
          }

          ValueType value(const PointType& point) const
          {
            ValueType val(DataType(0));
            const DataType x = point[0] - DataType(0.5);
            const DataType y = point[1] - DataType(0.5);
            const DataType r = Math::sqrt(x*x + y*y);
            if (r < DataType(0.2))
            {
              val[0] =  DataType(5)*y;
              val[1] = -DataType(5)*x;
            }
            else if (r < DataType(0.4))
            {
              val[0] = ( DataType(2)/r - DataType(5))*y;
              val[1] = (-DataType(2)/r + DataType(5))*x;
            }
            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad(DataType(0));
            const DataType x = point[0] - DataType(0.5);
            const DataType y = point[1] - DataType(0.5);
            const DataType r = Math::sqrt(x*x + y*y);
            if (r < DataType(0.2))
            {
              grad[0][1] =  DataType(5);
              grad[1][0] = -DataType(5);
            }
            else if (r < DataType(0.4))
            {
              const DataType s = DataType(1) / (r*r*r);
              grad[0][0] = -DataType(2)*x*y*s;
              grad[0][1] = -DataType(2)*y*y*s + DataType(2)/r - DataType(5);
              grad[1][0] =  DataType(2)*x*x*s - DataType(2)/r + DataType(5);
              grad[1][1] =  DataType(2)*x*y*s;
            }
            return grad;
          }
        };
      }; // class StandingVortexFunction2D

      /**
       * \brief Taylor-Green Vortex velocity field
       *
       * This class implements the Taylor-Green vortex velocity field,
       * which is the velocity field of an analytical solution to the
       * incompressible 2D Navier-Stokes equations on the unit-square.
       * The compatible pressure function is TaylorGreenVortexPres2D.
       *
       * \f[v(t,x,y) := \begin{bmatrix}\sin(\pi x)\cos(\pi y)\exp(-2\pi^2 \nu t)\\
         -\cos(\pi x)\sin(\pi y)\exp(-2\pi^2 \nu t)\end{bmatrix}\f]
       *
       * \author Peter Zajac
       */
      template<typename DT_>
      class TaylorGreenVortexVelo2D : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<2> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 2;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can't compute the Hessian
        static constexpr bool can_hess = true;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /// The viscosity parameter
        DataType nu;
        /// The current simulation time
        DataType cur_t;

        /**
         * \brief Constructor
         *
         * \param[in] nu_
         * The viscosity parameter nu.
         *
         * \param[in] t_
         * The current simulation time t.
         */
        explicit TaylorGreenVortexVelo2D(DataType nu_ = DataType(1), DataType t_ = DataType(0)) :
          nu(nu_),
          cur_t(t_)
        {
        }

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
          const DataType _pi;
          const DataType _fac;

        public:
          explicit Evaluator(const TaylorGreenVortexVelo2D& function) :
            _pi(Math::pi<DataType>()),
            _fac(Math::exp(-DataType(2)*_pi*_pi*function.nu*function.cur_t))
          {
          }

          ValueType value(const PointType& point)
          {
            ValueType val = ValueType::null();
            val[0] = +_fac * Math::sin(_pi*point[0]) * Math::cos(_pi*point[1]);
            val[1] = -_fac * Math::cos(_pi*point[0]) * Math::sin(_pi*point[1]);
            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad = GradientType::null();
            grad[0][0] = +_fac * _pi * Math::cos(_pi*point[0]) * Math::cos(_pi*point[1]);
            grad[0][1] = -_fac * _pi * Math::sin(_pi*point[0]) * Math::sin(_pi*point[1]);
            grad[1][0] = +_fac * _pi * Math::sin(_pi*point[0]) * Math::sin(_pi*point[1]);
            grad[1][1] = -_fac * _pi * Math::cos(_pi*point[0]) * Math::cos(_pi*point[1]);
            return grad;
          }

          HessianType hessian(const PointType& point) const
          {
            HessianType hess = HessianType::null();
            hess[0][0][0] = -_fac * _pi*_pi * Math::sin(_pi*point[0]) * Math::cos(_pi*point[1]);
            hess[0][0][1] = -_fac * _pi*_pi * Math::cos(_pi*point[0]) * Math::sin(_pi*point[1]);
            hess[0][1][0] = -_fac * _pi*_pi * Math::cos(_pi*point[0]) * Math::sin(_pi*point[1]);
            hess[0][1][1] = -_fac * _pi*_pi * Math::sin(_pi*point[0]) * Math::cos(_pi*point[1]);
            hess[1][0][0] = +_fac * _pi*_pi * Math::cos(_pi*point[0]) * Math::sin(_pi*point[1]);
            hess[1][0][1] = +_fac * _pi*_pi * Math::sin(_pi*point[0]) * Math::cos(_pi*point[1]);
            hess[1][1][0] = +_fac * _pi*_pi * Math::sin(_pi*point[0]) * Math::cos(_pi*point[1]);
            hess[1][1][1] = +_fac * _pi*_pi * Math::cos(_pi*point[0]) * Math::sin(_pi*point[1]);
            return hess;
          }
        }; // class TaylorGreenVortexVelo::Evaluator<...>
      }; // class TaylorGreenVortexVelo<...>

      /**
       * \brief Taylor-Green Vortex pressure function
       *
       * This class implements the Taylor-Green vortex pressure function,
       * which is the pressure function of an analytical solution to the
       * incompressible 2D Navier-Stokes equations on the unit-square.
       * The compatible pressure function is TaylorGreenVortexVelo2D.
       *
       * \f[p(t,x,y) := \frac{1}{2}\Big(\cos(\pi x)^2 + \cos(\pi y)^2 -1\Big)\exp(-4\pi^2 \nu t)\f]
       *
       * \author Peter Zajac
       */
      template<typename DT_>
      class TaylorGreenVortexPres2D : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Scalar ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 2;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /// The viscosity parameter
        DataType nu;
        /// The current simulation time
        DataType cur_t;

        /**
         * \brief Constructor
         *
         * \param[in] nu_
         * The viscosity parameter nu.
         *
         * \param[in] t_
         * The current simulation time t.
         */
        explicit TaylorGreenVortexPres2D(DataType nu_ = DataType(1), DataType t_ = DataType(0)) :
          nu(nu_),
          cur_t(t_)
        {
        }

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
          const DataType _pi;
          const DataType _fac;

        public:
          explicit Evaluator(const TaylorGreenVortexPres2D& function) :
            _pi(Math::pi<DataType>()),
            _fac(Math::exp(-DataType(4)*_pi*_pi*function.nu*function.cur_t) * DataType(0.5))
          {
          }

          ValueType value(const PointType& point)
          {
            return ValueType(_fac * (Math::sqr(Math::cos(_pi*point[0])) + Math::sqr(Math::cos(_pi*point[1])) - DataType(1)));
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad = GradientType::null();
            grad[0] = -DataType(2) * _pi * _fac * Math::sin(_pi*point[0]) * Math::cos(_pi*point[0]);
            grad[1] = -DataType(2) * _pi * _fac * Math::sin(_pi*point[1]) * Math::cos(_pi*point[1]);
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
          }
        }; // class TaylorGreenVortexPres::Evaluator<...>
      }; // class TaylorGreenVortexPres<...>

      /**
       * \brief Parabolic Poiseuille Pipe-Flow velocity field
       *
       * This function is effectively a 3D equivalent of the ParProfileVector function, which represents a
       * parabolic flow inside a 3D cylinder. This function is an analytic steady state solution to the
       * 3D incompressible Navier-Stokes equations in a cylinder domain (with constant radius) and thus
       * represent a Poiseuille flow in this domain.
       *
       * This function is parameterized in the rotational axis and the radius of the cylinder domain as
       * well as the maximum flow velocity in the center of the cylinder.
       *
       * \author Peter Zajac
       */
      template<typename DataType_, int dim_>
      class PoiseuillePipeFlow :
        public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DataType_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<dim_> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = dim_;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can't compute the Hessian
        static constexpr bool can_hess = true;
        /// Type to map from
        typedef Tiny::Vector<DataType, domain_dim> PointType;

        /// The flow origin
        PointType origin;
        /// The flow axis
        PointType axis;
        /// The pipe radius
        DataType radius;
        /// the maximum flow velocity
        DataType v_max;

        /**
         * \brief Creates this function
         *
         * \param[in] origin_
         * Specifies the origin of the cylinder domain, i.e. some arbitrary point along the rotational axis
         *
         * \param[in] axis_
         * Specifies the rotational axis of the cylinder domain. Must not be the null vector.
         *
         * param[in] radius_
         * Specifies the radius of the cylinder domain. Must be > 0.
         *
         * \param[in] v_max_
         * Specifies the maximum flow velocity in the center of the cylinder. Should be > 0.
         */
        explicit PoiseuillePipeFlow(PointType origin_, PointType axis_, DataType radius_ = DataType(1), DataType v_max_ = DataType(1)) :
          origin(origin_),
          axis(axis_),
          radius(radius_),
          v_max(v_max_)
        {
        }

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
          const PointType _origin;
          const PointType _axis;
          const DataType _v_max;
          const DataType _radius_s;

        public:
          explicit Evaluator(const PoiseuillePipeFlow& function) :
            _origin(function.origin),
            _axis(PointType(function.axis).normalize()),
            _v_max(function.v_max),
            _radius_s(DataType(1) / Math::sqr(function.radius))
          {
          }

          ValueType value(const PointType& point)
          {
            // compute ray to project point onto axis
            const PointType q = _origin + Tiny::dot(point - _origin, _axis) * _axis - point;
            // compute velocity based on distance to axis and in direction of axis
            return _v_max * (DataType(1) - q.norm_euclid_sqr()*_radius_s) * _axis;
          }

          GradientType gradient(const PointType& point) const
          {
            // compute ray to project point onto axis
            const PointType q = _origin + Tiny::dot(point - _origin, _axis) * _axis - point;
            // compute gradient
            GradientType grad = GradientType::null();
            //grad.add_outer_product(_axis, Tiny::dot(_axis, q)*_axis - q, -DataType(2) * _radius_s * _v_max);
            grad.set_outer_product(_axis, q);
            grad.add_outer_product(_axis, _axis, -Tiny::dot(_axis, q));
            return DataType(2) * _radius_s * _v_max * grad;
          }

          HessianType hessian(const PointType&) const
          {
            // seriously, getting this mumbo jumbo below right was mostly a lucky guess...
            GradientType grad_q = GradientType::null();
            grad_q.set_outer_product(_axis, _axis);
            grad_q.add_scalar_main_diag(-DataType(1));
            GradientType at_x_a;
            at_x_a.set_outer_product(_axis, _axis);
            HessianType hess = HessianType::null();
            hess.add_vec_mat_outer_product(_axis, grad_q);
            hess.add_vec_mat_outer_product(_axis * grad_q, at_x_a, -DataType(1));
            return DataType(2) * _radius_s * _v_max * hess;
          }
        }; // class PoiseuillePipeFlow::Evaluator<...>
      }; // class PoiseuillePipeFlow

      /**
       * \brief Rigid-Body Vortex velocity field
       *
       * This class implements the Rigid-Body vortex velocity field,
       * which is the velocity field of an analytical solution to the
       * incompressible 2D Navier-Stokes equations on the unit-circle.
       * The compatible pressure function is RigidBodyVortexPres2D.
       *
       * \f[v(x,y) := \begin{bmatrix}-y\\x\end{bmatrix}\f]
       *
       * \author Peter Zajac
       */
      template<typename DT_>
      class RigidBodyVortexVelo2D : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<2> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 2;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can't compute the Hessian
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

        public:
          explicit Evaluator(const RigidBodyVortexVelo2D&)
          {
          }

          ValueType value(const PointType& point)
          {
            ValueType val = ValueType::null();
            val[0] = -point[1];
            val[1] =  point[0];
            return val;
          }

          GradientType gradient(const PointType& DOXY(point)) const
          {
            GradientType grad = GradientType::null();
            grad[0][1] = -DataType(1);
            grad[1][0] =  DataType(1);
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
          }
        }; // class RigidBodyVortexVelo2D::Evaluator<...>
      }; // class RigidBodyVortexVelo2D<...>

      /**
       * \brief Rigid-Body Vortex pressure function
       *
       * This class implements the Rigid-Body vortex pressure function,
       * which is the pressure function of an analytical solution to the
       * incompressible 2D Navier-Stokes equations on the unit-square.
       * The compatible pressure function is RigidBodyVortexVelo2D.
       *
       * \f[p(x,y) := \frac{1}{2}\Big(x^2 + y^2\Big) - \frac{1}{4}\f]
       *
       * \author Peter Zajac
       */
      template<typename DT_>
      class RigidBodyVortexPres2D : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Scalar ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 2;
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

        public:
          explicit Evaluator(const RigidBodyVortexPres2D&)
          {
          }

          ValueType value(const PointType& point)
          {
            return DataType(0.5) * (point[0]*point[0] + point[1]*point[1]) - DataType(0.25);
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad = GradientType::null();
            grad[0] = point[0];
            grad[1] = point[1];
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            HessianType hess = HessianType::null();
            hess[0][0] = DataType(1);
            hess[1][1] = DataType(1);
            return hess;
          }
        }; // class RigidBodyVortexPres2D::Evaluator<...>
      }; // class RigidBodyVortexPres2D<...>

      /**
       * \brief Sine-Vortex velocity field on the 2D ring domain with distance [1/2,1] around (0,0)
       *
       * This function can be used as an analytic solution to the steady-state Stokes or Navier-Stokes
       * equations on the 2D ring domain with distance [1/2,1] around the origin (0,0) when used in
       * combination with the SineRingVortexPres2D and SineRingVortexRHS2D functions for pressure and
       * right-hand-side, respectively.
       *
       * The velocity field is defined as
       *
       * \f[v(x,y) := \frac{1}{d(x,y)} \begin{bmatrix}\hphantom{-}y\sin(2\pi d(x,y))\\ -x\sin(2\pi d(x,y))\end{bmatrix}\f]
       *
       * where \f$d(x,y) := \sqrt{x^2+y^2}\f$ is the euclidean distance to the origin.
       *
       * \author Peter Zajac
       */
      template<typename DT_>
      class SineRingVortexVelo2D : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<2> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 2;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can't compute the Hessian
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
          /// 2*pi
          const DataType _pi2;

        public:
          explicit Evaluator(const SineRingVortexVelo2D&) :
            _pi2(DataType(2) * Math::pi<DataType>())
          {
          }

          ValueType value(const PointType& point)
          {
            ValueType val = ValueType::null();
            const DataType d = point.norm_euclid();
            if(d < DataType(1E-3))
              return val;

            const DataType s2pdi = Math::sin(_pi2 * d) / d;

            val[0] =  point[1]*s2pdi;
            val[1] = -point[0]*s2pdi;
            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad = GradientType::null();
            const DataType d = point.norm_euclid();
            if(d < DataType(1E-3))
              return grad;

            const DataType x(point[0]), y(point[1]);
            const DataType di2 = DataType(1) / (d*d);
            const DataType s2pdi = Math::sin(_pi2 * d) / d;
            const DataType c2pdi2 = Math::cos(_pi2 * d) * di2;

            grad[0][0] = -y*x*s2pdi*di2 + _pi2*y*x*c2pdi2;
            grad[0][1] = -y*y*s2pdi*di2 + _pi2*y*y*c2pdi2 + s2pdi;
            grad[1][0] =  x*x*s2pdi*di2 - _pi2*x*x*c2pdi2 - s2pdi;
            grad[1][1] =  x*y*s2pdi*di2 - _pi2*y*x*c2pdi2;
            return grad;
          }

          HessianType hessian(const PointType& point) const
          {
            HessianType hess = HessianType::null();
            const DataType d = point.norm_euclid();
            if(d < DataType(1E-3))
              return hess;

            const DataType x(point[0]), y(point[1]);
            const DataType di2 = DataType(1) / (d*d);
            const DataType s2pdi3 = Math::sin(_pi2 * d) * di2 / d;
            const DataType c2pdi2 = Math::cos(_pi2 * d) * di2;
            const DataType t3 = DataType(3);

            hess[0][0][0] =  t3*y*x*x*s2pdi3*di2 - t3*_pi2*y*x*x*c2pdi2*di2 - y*s2pdi3 + _pi2*y*c2pdi2 - _pi2*_pi2*y*x*x*s2pdi3;
            hess[0][0][1] =
            hess[0][1][0] =  t3*y*y*x*s2pdi3*di2 - x*s2pdi3 - t3*_pi2*y*y*x*c2pdi2*di2 + _pi2*x*c2pdi2 - _pi2*_pi2*y*y*x*s2pdi3;
            hess[0][1][1] =  t3*y*y*y*s2pdi3*di2 - t3*y*s2pdi3 - t3*_pi2*y*y*y*c2pdi2*di2 + t3*_pi2*y*c2pdi2 - _pi2*_pi2*y*y*y*s2pdi3;
            hess[1][0][0] = -t3*x*x*x*s2pdi3*di2 + t3*x*s2pdi3 + t3*_pi2*x*x*x*c2pdi2*di2 - t3*_pi2*x*c2pdi2 + _pi2*_pi2*x*x*x*s2pdi3;
            hess[1][0][1] =
            hess[1][1][0] = -t3*y*x*x*s2pdi3*di2 + t3*_pi2*y*x*x*c2pdi2*di2 + y*s2pdi3 - _pi2*y*c2pdi2 + _pi2*_pi2*y*x*x*s2pdi3;
            hess[1][1][1] = -t3*y*y*x*s2pdi3*di2 + x*s2pdi3 + t3*_pi2*y*y*x*c2pdi2*di2 - _pi2*x*c2pdi2 + _pi2*_pi2*y*y*x*s2pdi3;
            return hess;
          }
        }; // class SineRingVortexVelo2D::Evaluator<...>
      }; // class SineRingVortexVelo2D<...>

      /**
       * \brief Sine-Vortex pressure function on the 2D ring domain with distance [1/2,1] around (0,0)
       *
       * This function can be used as an analytic solution to the steady-state Stokes or Navier-Stokes
       * equations on the 2D ring domain with distance [1/2,1] around the origin (0,0) when used in
       * combination with the SineRingVortexVelo2D and SineRingVortexRHS2D functions for velocity and
       * right-hand-side, respectively.
       *
       * The pressure function is defined as
       *
       * \f[p(x,y) := \frac{\sin(2\pi d(x,y))}{2\pi}\f]
       *
       * where \f$d(x,y) := \sqrt{x^2+y^2}\f$ is the euclidean distance to the origin.
       *
       * \author Peter Zajac
       */
      template<typename DT_>
      class SineRingVortexPres2D : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Scalar ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 2;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = true;
        /// We can compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /// the shift of the mean value; defaults to that the pressure has integral mean equal to 0
        DataType mean_shift = DataType(1) / Math::sqr(Math::pi<DataType>());

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
          /// 2*pi and 1/pi^2
          const DataType _pi2, _spi;

        public:
          explicit Evaluator(const SineRingVortexPres2D& func) :
            _pi2(DataType(2) * Math::pi<DataType>()),
            _spi(func.mean_shift)
          {
          }

          ValueType value(const PointType& point)
          {
            const DataType d = point.norm_euclid();
            return -Math::sin(_pi2*d) / _pi2 - _spi;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad = GradientType::null();
            const DataType d = point.norm_euclid();
            if(d < DataType(1E-3))
              return grad;

            const DataType c2pdi = Math::cos(_pi2 * d)  / d;

            grad[0] = -point[0]*c2pdi;
            grad[1] = -point[1]*c2pdi;
            return grad;
          }

          HessianType hessian(const PointType& DOXY(point)) const
          {
            return HessianType::null();
          }
        }; // class SineRingVortexPres2D::Evaluator<...>
      }; // class SineRingVortexPres2D<...>

      /**
       * \brief Sine-Vortex RHS field on the 2D ring domain with distance [1/2,1] around (0,0)
       *
       * This function can be used as the right-hand-side for the momentum equation to the steady-state Stokes or
       * Navier-Stokes equations on the 2D ring domain with distance [1/2,1] around the origin (0,0) which results
       * in a system whose analytical solution is given by the SineRingVortexVelo2D and SineRingVortexPres2D functions,
       * respectively.
       *
       * \todo check whether this is still the correct RHS if using deformation tensor formulation (should be the case)
       *
       * \author Peter Zajac
       */
      template<typename DT_>
      class SineRingVortexRHS2D : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef Analytic::Image::Vector<2> ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 2;
        /// We can compute the value
        static constexpr bool can_value = true;
        /// We can compute the gradient
        static constexpr bool can_grad = false;
        /// We can't compute the Hessian
        static constexpr bool can_hess = false;
        /// Type to map from
        typedef Tiny::Vector<DT_, domain_dim> PointType;

        /// coefficients
        DataType nu, beta, theta, sigma;

        /**
         * \brief Constructor
         *
         * \param[in] nu_
         * The viscosity parameter for the velocity diffusion term
         *
         * \param[in] beta_
         * The convection coefficient; 0 for Stokes and 1 for Navier-Stokes
         *
         * \param[in] theta_
         * The reaction coefficient; 0 for steady-state, but may be set > 0 for a quasi-Stokes RHS
         *
         * \param[in] sigma_
         * The factor for the pressure gradient; should be 1
         */
        explicit SineRingVortexRHS2D(DataType nu_ = DataType(1), DataType beta_ = DataType(0), DataType theta_ = DataType(0), DataType sigma_ = DataType(1)) :
          nu(nu_),
          beta(beta_),
          theta(theta_),
          sigma(sigma_)
        {
        }

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
          const DataType _pi2, _nu, _beta, _theta, _sigma;

        public:
          explicit Evaluator(const SineRingVortexRHS2D& function) :
            _pi2(DataType(2) * Math::pi<DataType>()),
            _nu(function.nu),
            _beta(function.beta),
            _theta(function.theta),
            _sigma(function.sigma)
          {
          }

          ValueType value(const PointType& point)
          {
            ValueType val = ValueType::null();
            const DataType d = point.norm_euclid();
            if(d < DataType(1E-3))
              return val;

            const DataType x(point[0]), y(point[1]);
            const DataType q = DataType(1) / d;
            const DataType s2pdq = Math::sin(_pi2 * d) * q;
            const DataType c2pdq = Math::cos(_pi2 * d) * q;

            val[0] = -_nu*q*y*( _pi2*c2pdq - (_pi2*_pi2*d + q)*s2pdq) - _beta*x*s2pdq*s2pdq + _theta*y*s2pdq - _sigma*x*c2pdq;
            val[1] = -_nu*q*x*(-_pi2*c2pdq + (_pi2*_pi2*d + q)*s2pdq) - _beta*y*s2pdq*s2pdq - _theta*x*s2pdq - _sigma*y*c2pdq;
            return val;
          }
        }; // class SineRingVortexRHS2D::Evaluator<...>
      }; // class SineRingVortexRHS2D<...>

      /**
       * \brief Sphere-Cap Function on the 2D unit circle / 3D unit ball
       *
       * This function is given by the formula \f$ (\sqrt{a - \|x\|_2^2} - \sqrt{a-1}) \f$,
       * which represents the upper half-sphere for a_ = 1. For a_ > 1, this function represents
       * just a part of the upper half-sphere, which is regular even for \f$\|x\|_2=1\f$.
       *
       * It is recommended to use this function with a_ = 2, which yields a sufficiently regular
       * function on the 2D/3D unit circle/ball.
       *
       * \author Peter Zajac
       */
      template<int dim_, int a_ = 2>
      class SphereCapFunction :
        public Analytic::Function
      {
      public:
        static_assert(a_ >= 1, "parameter a_ must be >= 1");

        static constexpr int domain_dim = dim_;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

        template<typename Traits_>
        class Evaluator :
          public Analytic::Function::Evaluator<Traits_>
        {
        public:
          typedef typename Traits_::DataType DataType;
          typedef typename Traits_::PointType PointType;
          typedef typename Traits_::ValueType ValueType;
          typedef typename Traits_::GradientType GradientType;
          typedef typename Traits_::HessianType HessianType;

          static constexpr DataType a = DataType(a_);
          const DataType sa;

          explicit Evaluator(const SphereCapFunction&) :
            sa(Math::sqrt(a - DataType(1)))
          {
          }

          ValueType value(const PointType& point) const
          {
            return Math::sqrt(a - point.norm_euclid_sqr()) - sa;
          }

          GradientType gradient(const PointType& point) const
          {
            return (DataType(-1) / Math::sqrt(a - point.norm_euclid_sqr())) * point;
          }

          HessianType hessian(const PointType& point) const
          {
            const DataType s = DataType(1) / Math::cub(Math::sqrt(a - point.norm_euclid_sqr()));
            //const DataType s = DataType(1) / Math::pow(a - point.norm_euclid_sqr(), DataType(1.5));
            HessianType hess;
            hess.format();
            for(int i = 0; i < dim_; ++i)
            {
              for(int j = 0; j < dim_; ++j)
              {
                if(i == j)
                {
                  hess[i][j] = DataType(0);
                  for(int k = 0; k < dim_; ++k)
                    hess[i][j] += (k == i ? -a : Math::sqr(point[k])) * s;
                }
                else
                {
                  hess[i][j] = -point[i] * point[j] * s;
                }
              }
            }
            return hess;
          }
        }; // class SphereCapFunction::Evaluator<...>
      }; // class SphereCapFunction

      /**
       * \brief Singularity function for a 2D reentry corner with homogeneous boundary conditions.
       *
       * This function represents the lowest order singular function that occurs for the poisson
       * equation with a re-entry corner, which has zero boundary conditions in the corner itself.
       *
       * The natural space to describe this function is in radial basis and is given by
       *         \f$ f(r,\varphi) = r^{\frac{\pi}{\theta}}\sin(\frac{\pi}{\theta}\varphi) \f$
       * where \theta is the (right-handed) opening angle of the re-entry corner.
       * This function is in \f$ H^1 \f$ but for degrees larger \f$ \frac{\pi}{2} \f$ not in
       * \f$ H^2 \f$ despite the corresponding Poisson-Equation having homogeneous right hand side
       * and arbitrarily smooth boundary conditions.
       *
       * The actual function space is provided in standard euclidean coordinates and requires a center point p
       * as well as the two directional vectors v_1, v_2 of the re-entry corner provided in a right-hand oriented manner:
       *
       *                  \ v_2
       *                   \
       *                    \  theta
       *                     + -----------------
       *                   p                  v_1
       *
       * \note This version directly calculates the values, gradients and hessians in euclidean basis, allowing for a few
       *       simplifications. There is also a variant based on radial basis function, which has to be wrapped in
       *       the PolarCoordinate wrapper.
       * \author Maximilian Esser
       */
      template<typename DT_>
      class CornerSingularity2D :
        public Analytic::Function
      {
      public:
        typedef DT_ DataType;
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

      private:
        Tiny::Vector<DataType, 2> base_point;
        DataType alpha, k, offset;

      public:
        CornerSingularity2D() = delete;

        CornerSingularity2D(const Tiny::Vector<DT_, 2>& center, Tiny::Vector<DT_, 2> v_1, Tiny::Vector<DT_, 2> v_2, DT_ alpha_ = DT_(1.)) :
          base_point(center),
          alpha(alpha_)
        {
          ASSERTM((v_1.norm_euclid_sqr() > Math::eps<DataType>()) && (v_2.norm_euclid_sqr() > Math::eps<DataType>()), "Spanvectors are zero length");
          k = Math::pi<DataType>() / Tiny::calculate_opening_angle(v_1,v_2);
          offset = Tiny::calculate_opening_angle(Tiny::Vector<DT_, 2>{DataType(1.0), DataType(0.)}, v_1);
        }

        template<typename Traits_>
        class Evaluator :
          public Analytic::Function::Evaluator<Traits_>
        {
        public:
          typedef typename Traits_::DataType DataType;
          typedef typename Traits_::PointType PointType;
          typedef typename Traits_::ValueType ValueType;
          typedef typename Traits_::GradientType GradientType;
          typedef typename Traits_::HessianType HessianType;


          Tiny::Vector<DataType, 2> base_point;
          const DataType alpha, k, offset;
          const DataType eps;

          explicit Evaluator(const CornerSingularity2D& func) :
            base_point(func.base_point),
            alpha(func.alpha),
            k(func.k),
            offset(func.offset),
            eps(Math::eps<DataType>())
          {
          }

          ValueType value(const PointType& point) const
          {
            //first calculate radius
            const DataType x = point[0] - base_point[0];
            const DataType y = point[1] - base_point[1];
            const DataType r = Math::sqrt(Math::sqr(x) + Math::sqr(y));
            //TODO: Set x_rel to zero if r is very small?? <- Better conditioning?
            const DataType x_rel = r > DataType(10.)* eps ? x/r : DataType(0);
            DataType phi = Math::acos(x_rel);
            phi = (y >= 0) ? phi : DataType(2)*Math::pi<DataType>()-phi;
            return alpha * Math::pow(r, k) * Math::sin(k * (phi - offset));
          }

          GradientType gradient(const PointType& point) const
          {
            const DataType x = point[0] - base_point[0];
            const DataType y = point[1] - base_point[1];
            const DataType r = Math::sqrt(Math::sqr(x) + Math::sqr(y));
            const DataType x_rel = r > DataType(10.)* eps ? x/r : DataType(0);
            const DataType y_rel = r > DataType(10.)* eps ? y/r : DataType(0);
            //Note: While in practice we use a rotated coordinate system to express our polar-coordinate
            //based function, we have to use the "real" cos(phi) and sin(phi) values
            //for our base transformation, which are conveniently x_rel and y_rel
            //TODO: This could be badly conditioned for r -> 0, have to test this...
            DataType phi = Math::acos(x_rel);
            phi = (y >= 0) ? phi : DataType(2)*Math::pi<DataType>()-phi;
            const DataType ex = k * Math::pow(r, k-1);
            const DataType sval = Math::sin(k * (phi - offset));
            const DataType cval = Math::cos(k * (phi - offset));
            GradientType grad;
            grad[0] = alpha * ex * (x_rel * sval - y_rel * cval);
            grad[1] = alpha * ex * (y_rel * sval + x_rel * cval);
            return grad;
          }

          HessianType hessian(const PointType& point) const
          {
            const DataType x = point[0] - base_point[0];
            const DataType y = point[1] - base_point[1];
            const DataType r = Math::sqrt(Math::sqr(x) + Math::sqr(y));
            const DataType x_rel = r > DataType(10.)* eps ? x/r : DataType(0);
            const DataType y_rel = r > DataType(10.)* eps ? y/r : DataType(0);
            DataType phi = Math::acos(x_rel);
            phi = (y >= 0) ? phi : DataType(2)*Math::pi<DataType>()-phi;
            const DataType ex = k * (k - DataType(1)) * Math::pow(r, k-2);
            const DataType sval = Math::sin(k * (phi - offset));
            const DataType cval = Math::cos(k * (phi - offset));
            HessianType hess;
            hess[0][0] = alpha * ex * ((Math::sqr(x_rel) - Math::sqr(y_rel)) * sval - DataType(2) * x_rel * y_rel * cval);
            hess[1][1] = -hess[0][0];
            hess[0][1] = hess[1][0] = alpha * ex * (DataType(2) * y_rel * x_rel * sval + (Math::sqr(x_rel) - Math::sqr(y_rel)) * cval);
            return hess;
          }
        }; // class CornerSingularity2D::Evaluator<...>
      }; // CornerSingularity2D

      /// radial version of the CornerSingularity function... should be used with the PolarCoordinateWrapper
      template<typename DT_>
      class CornerSingularity2DRadial :
        public Analytic::Function
      {
      public:
        typedef DT_ DataType;
        static constexpr int domain_dim = 2;
        typedef Analytic::Image::Scalar ImageType;
        static constexpr bool can_value = true;
        static constexpr bool can_grad = true;
        static constexpr bool can_hess = true;

      private:
        DataType alpha, k;

      public:
        CornerSingularity2DRadial() = delete;

        CornerSingularity2DRadial(DT_ theta_, DT_ alpha_ = DT_(1.)) :
          alpha(alpha_),
          k(Math::pi<DataType>() / theta_)
        {}

        /// ctor, theta is the angle required to rotate vec_r into vec_l in a counterclockwise manner
        CornerSingularity2DRadial(const Tiny::Vector<DT_, 2>& vec_r, const Tiny::Vector<DT_, 2>& vec_l, DT_ alpha_ = DT_(1.)) :
          alpha(alpha_),
          k(Math::pi<DataType>() / Math::calc_opening_angle(vec_r[0], vec_r[1], vec_l[0], vec_l[1]))
        {}

        template<typename Traits_>
        class Evaluator :
          public Analytic::Function::Evaluator<Traits_>
        {
        public:
          typedef typename Traits_::DataType DataType;
          typedef typename Traits_::PointType PointType;
          typedef typename Traits_::ValueType ValueType;
          typedef typename Traits_::GradientType GradientType;
          typedef typename Traits_::HessianType HessianType;


          const DataType alpha, k;

          explicit Evaluator(const CornerSingularity2DRadial& func) :
            alpha(func.alpha),
            k(func.k)
          {
          }

          ValueType value(const PointType& point) const
          {
            return alpha * Math::pow(point[0], k) * Math::sin(k * point[1]);
          }

          GradientType gradient(const PointType& point) const
          {
            const DataType ex = k * Math::pow(point[0], k-1);
            const DataType sval = Math::sin(k * point[1]);
            const DataType cval = Math::cos(k * point[1]);
            GradientType grad;
            grad[0] = alpha * ex * sval;
            grad[1] = alpha * ex * cval * point[0];
            return grad;
          }

          HessianType hessian(const PointType& point) const
          {
            const DataType ex = k * Math::pow(point[0], k-2);
            const DataType sval = Math::sin(k * point[1]);
            const DataType cval = Math::cos(k * point[1]);
            HessianType hess;
            hess[0][0] = alpha * ex * (k-1) * sval;
            hess[1][1] = - alpha * ex * Math::sqr(point[0]) * k * sval;
            hess[0][1] = hess[1][0] = alpha * ex * k * point[0] * cval;
            return hess;
          }
        }; // class CornerSingularity2DRadial::Evaluator<...>

      }; // class CornerSingularity2DRadial

      template<typename DT_>
      using CornerSingluratity2DSimple = PolarCoordinate<CornerSingularity2DRadial<DT_>>;

      /**
       * \brief Saddle-Point-Gauss type interpolation test function
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * For details, see Franke, R. (1979). A critical comparison of some methods for interpolation of scattered data
       *
       * \author Maximilian Esser
       */
      template<typename DT_>
      class FrankesFunction : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef typename Analytic::Image::Scalar ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 2;
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

        public:
          explicit Evaluator(const FrankesFunction&)
          {
          }

          ValueType value(const PointType& point)
          {
            const DataType x = point[0];
            const DataType y = point[1];

            DataType val = (DataType(3) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)))) / DataType(4) - Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) / DataType(5) +(DataType(3) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) -
                  Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)))) / DataType(4) + Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) / DataType(2);

            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad;
            const DataType x = point[0];
            const DataType y = point[1];
            grad.format();
            grad[0] = (Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) *(DataType(162) * x - DataType(72))) / DataType(5) -(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(9))) / DataType(4) -
                      (Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(63) / DataType(2))) / DataType(2) -(DataType(3) * Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) /
                      DataType(10)) *((DataType(162) * x) / DataType(49) + DataType(18) / DataType(49))) / DataType(4);
            grad[1] = (Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) *(DataType(162) * y - DataType(126))) / DataType(5) -(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)) *((DataType(81) * y) / DataType(2) - DataType(9))) / DataType(4) -
                      (Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) *((DataType(81) * y) / DataType(2) - DataType(27) / DataType(2))) / DataType(2) -(DataType(27) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) /
                      DataType(10)))) / DataType(40);
            return grad;
          }

          HessianType hessian(const PointType& point) const
          {
            HessianType hess;
            const DataType x = point[0];
            const DataType y = point[1];
            hess.format();
            hess[0][0] = (DataType(162) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))))) / DataType(5) -(DataType(243) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)))) / DataType(98) -(DataType(243) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) /
                          DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)))) / DataType(8) -(DataType(81) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)))) / DataType(4) +(DataType(3) * Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) /
                          DataType(49) - DataType(1) / DataType(10)) * Math::sqr((DataType(162) * x) / DataType(49) + DataType(18) / DataType(49))) / DataType(4) +(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)) * Math::sqr((DataType(81) * x) / DataType(2) - DataType(9))) / DataType(4) +(Math::exp(-
                          Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) * Math::sqr((DataType(81) * x) / DataType(2) - DataType(63) / DataType(2))) / DataType(2) -(Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) * Math::sqr(DataType(162) * x - DataType(72))) / DataType(5);
            hess[0][1] = hess[1][0] = (DataType(27) * Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)) *((DataType(162) * x) / DataType(49) + DataType(18) / DataType(49))) / DataType(40) +(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) /
                          DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(9)) *((DataType(81) * y) / DataType(2) - DataType(9))) / DataType(4) +(Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(63) / DataType(2)) *((DataType(81) *
                          y) / DataType(2) - DataType(27) / DataType(2))) / DataType(2) -(Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) *(DataType(162) * x - DataType(72)) *(DataType(162) * y - DataType(126))) / DataType(5);
            hess[1][1] = (DataType(243) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)))) / DataType(400) +(DataType(162) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))))) / DataType(5) -(DataType(243) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) /
                          DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)))) / DataType(8) -(DataType(81) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)))) / DataType(4) +(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y -
                          DataType(2)) / DataType(4)) * Math::sqr((DataType(81) * y) / DataType(2) - DataType(9))) / DataType(4) +(Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) * Math::sqr((DataType(81) * y) / DataType(2) - DataType(27) / DataType(2))) / DataType(2) -(Math::exp(- Math::sqr(DataType(9) * x -
                          DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) * Math::sqr(DataType(162) * y - DataType(126))) / DataType(5);
            return hess;
          }

        }; // class FrankesFunction::Evaluator<...>

      public:
        explicit FrankesFunction() = default;
      }; // class FrankesFunction

      /**
       * \brief 3D Version of Franke's Function
       *
       * \tparam DT_
       * The floating point type
       *
       * \tparam dim_ The space dimension
       *
       * For Franke's function f(x,y), this function represents \f$ g(x,y,z) = \exp(-(10*z-2)^2/25) \cdot f(x,y) $\f
       *
       * \author Maximilian Esser
       */
      template<typename DT_>
      class Frankes3DVariantFunction : public Analytic::Function
      {
      public:
        /// The floating point type
        typedef DT_ DataType;
        /// What type this mapping maps to
        typedef typename Analytic::Image::Scalar ImageType;
        /// The dimension to map from
        static constexpr int domain_dim = 3;
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

        public:
          explicit Evaluator(const Frankes3DVariantFunction&)
          {
          }

          ValueType value(const PointType& point)
          {
            const DataType x = point[0];
            const DataType y = point[1];
            const DataType z = point[2];

            DataType val = (DataType(3) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)))) / DataType(4) - Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) / DataType(5) +(DataType(3) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) -
                  Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)))) / DataType(4) + Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) / DataType(2);
            val *= Math::exp(-Math::sqr(DataType(2)*z-DataType(0.4)));

            return val;
          }

          GradientType gradient(const PointType& point) const
          {
            GradientType grad;
            const DataType x = point[0];
            const DataType y = point[1];
            const DataType z = point[2];

            const DataType franks_val = (DataType(3) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)))) / DataType(4) - Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) / DataType(5) +(DataType(3) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) -
                  Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)))) / DataType(4) + Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) / DataType(2);


            const DataType zexp = Math::exp(-Math::sqr(DataType(2)*z-DataType(0.4)));
            const DataType dzexp = -DataType(4)*(DataType(2)*z - DataType(0.4)) * zexp;

            grad.format();
            grad[0] = zexp * ((Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) *(DataType(162) * x - DataType(72))) / DataType(5) -(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(9))) / DataType(4) -
                      (Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(63) / DataType(2))) / DataType(2) -(DataType(3) * Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) /
                      DataType(10)) *((DataType(162) * x) / DataType(49) + DataType(18) / DataType(49))) / DataType(4));
            grad[1] = zexp * ((Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) *(DataType(162) * y - DataType(126))) / DataType(5) -(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)) *((DataType(81) * y) / DataType(2) - DataType(9))) / DataType(4) -
                      (Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) *((DataType(81) * y) / DataType(2) - DataType(27) / DataType(2))) / DataType(2) -(DataType(27) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) /
                      DataType(10)))) / DataType(40));
            grad[2] = dzexp * franks_val;
            return grad;
          }

          HessianType hessian(const PointType& point) const
          {
            HessianType hess;
            const DataType x = point[0];
            const DataType y = point[1];
            const DataType z = point[2];
            const DataType zexp = Math::exp(-Math::sqr(DataType(2)*z-DataType(0.4)));
            const DataType dzexp = -DataType(4)*(DataType(2)*z - DataType(0.4)) * zexp;
            const DataType franks_val = (DataType(3) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)))) / DataType(4) - Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) / DataType(5) +(DataType(3) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) -
                  Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)))) / DataType(4) + Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) / DataType(2);
            hess.format();
            hess[0][0] = zexp * ((DataType(162) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))))) / DataType(5) -(DataType(243) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)))) / DataType(98) -(DataType(243) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) /
                          DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)))) / DataType(8) -(DataType(81) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)))) / DataType(4) +(DataType(3) * Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) /
                          DataType(49) - DataType(1) / DataType(10)) * Math::sqr((DataType(162) * x) / DataType(49) + DataType(18) / DataType(49))) / DataType(4) +(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)) * Math::sqr((DataType(81) * x) / DataType(2) - DataType(9))) / DataType(4) +(Math::exp(-
                          Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) * Math::sqr((DataType(81) * x) / DataType(2) - DataType(63) / DataType(2))) / DataType(2) -(Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) * Math::sqr(DataType(162) * x - DataType(72))) / DataType(5));
            hess[0][1] = hess[1][0] = zexp * ((DataType(27) * Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)) *((DataType(162) * x) / DataType(49) + DataType(18) / DataType(49))) / DataType(40) +(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) /
                          DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(9)) *((DataType(81) * y) / DataType(2) - DataType(9))) / DataType(4) +(Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(63) / DataType(2)) *((DataType(81) *
                          y) / DataType(2) - DataType(27) / DataType(2))) / DataType(2) -(Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) *(DataType(162) * x - DataType(72)) *(DataType(162) * y - DataType(126))) / DataType(5));
            hess[1][1] = zexp * ((DataType(243) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) / DataType(10)))) / DataType(400) +(DataType(162) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))))) / DataType(5) -(DataType(243) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) /
                          DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)))) / DataType(8) -(DataType(81) *(Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)))) / DataType(4) +(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y -
                          DataType(2)) / DataType(4)) * Math::sqr((DataType(81) * y) / DataType(2) - DataType(9))) / DataType(4) +(Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) * Math::sqr((DataType(81) * y) / DataType(2) - DataType(27) / DataType(2))) / DataType(2) -(Math::exp(- Math::sqr(DataType(9) * x -
                          DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) * Math::sqr(DataType(162) * y - DataType(126))) / DataType(5));
            hess[0][2] = hess[2][0] = dzexp * ((Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) *(DataType(162) * x - DataType(72))) / DataType(5) -(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(9))) / DataType(4) -
                          (Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) *((DataType(81) * x) / DataType(2) - DataType(63) / DataType(2))) / DataType(2) -(DataType(3) * Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) /
                          DataType(10)) *((DataType(162) * x) / DataType(49) + DataType(18) / DataType(49))) / DataType(4));
            hess[1][2] = hess[2][1] = dzexp * ((Math::exp(- Math::sqr(DataType(9) * x - DataType(4)) - Math::sqr(DataType(9) * y - DataType(7))) *(DataType(162) * y - DataType(126))) / DataType(5) -(DataType(3) * Math::exp(- Math::sqr(DataType(9) * x - DataType(2)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(2)) / DataType(4)) *((DataType(81) * y) / DataType(2) - DataType(9))) / DataType(4) -
                      (Math::exp(- Math::sqr(DataType(9) * x - DataType(7)) / DataType(4) - Math::sqr(DataType(9) * y - DataType(3)) / DataType(4)) *((DataType(81) * y) / DataType(2) - DataType(27) / DataType(2))) / DataType(2) -(DataType(27) *(Math::exp(-(DataType(9) * y) / DataType(10) - Math::sqr(DataType(9) * x + DataType(1)) / DataType(49) - DataType(1) /
                      DataType(10)))) / DataType(40));
            hess[2][2] = franks_val * (Math::sqr(DataType(4)*(DataType(2)*z - DataType(0.4))) - DataType(8)) * zexp;
            return hess;
          }

        }; // class Frankes3DVariantFunction::Evaluator<...>

      public:
        explicit Frankes3DVariantFunction() = default;
      }; // class FrankesFunction

      /**
       * \brief Harmonic shell function
       *
       * This class implements a harmonic function u, which fulfills u(x) = 1 for |x - x_0|_2 = 1 and u(x) = 2 for |x - x_0|_2 = 1/2.
       *
       * This function is given by the following formulae:
       * - 1D: u(x) := 3 - 2*|x-x_0|
       * - 2D: u(x,y) := 1 - ln((x-x_0)^2 + (y-y_0)^2)/(2*ln(2))
       * - 3D: u(x,y,z) := 1/sqrt((x-x_0)^2 + (y-y_0)^2 + (z-z0)^2)
       *
       * \note This function is harmonic and it is basically a scaled Green's function for the Laplace operator,
       * therefore it can be used as an analytical solution for the Laplace equation on the unit ring/shell domain
       * with Dirichlet boundary conditions u=1 for the outer boundary and u=2 for the inner boundary.
       *
       * \warning Because the function is differentiable everywhere except for x_0, Bad Things (TM) might happen if
       * someone wants to compute the gradient or hessian there, as the functions return 0.
       *
       * \tparam ImgPointType_
       * The type of \f$ x_0 \f$.
       *
       * \author Peter Zajac
       */
      template<int dim_, typename DataType_>
      class HarmonicShellFunction :
        public Analytic::Function
      {
      public:
        static_assert((1 <= dim_) && (dim_ <= 3), "harmonic shell function is only implemented in 1D, 2D and 3D");

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
          const DataType _eps;
          const DataType _ln2; // = ln(2)
          const DataType _i2ln2; // = 1 /(2*ln(2))

        public:
          /// Constructor
          explicit Evaluator(const HarmonicShellFunction& function) :
            _origin(function._origin),
            _eps(Math::eps<DataType>()),
            _ln2(Math::log(DataType(2))),
            _i2ln2(DataType(0.5)/_ln2)
          {
          }

          ValueType value(const PointType& point) const
          {
            if constexpr (dim_ == 1)
              return DataType(3)  - DataType(2) * Math::abs(point[0] - _origin[0]);
            else if constexpr(dim_ == 2)
              return DataType(1) - _i2ln2 * Math::log((point - _origin).norm_euclid_sqr());
            else if constexpr(dim_ == 3)
              return DataType(1) / (point - _origin).norm_euclid();
            else
              return DataType(0);
          }

          GradientType gradient(const PointType& point) const
          {
            const DataType norm_sqr = (point - _origin).norm_euclid_sqr();
            if(norm_sqr < _eps)
              return GradientType::null();

            GradientType grad(DataType(0));
            if constexpr (dim_ == 1)
              grad[0] = DataType(-2) * Math::signum(point[0] - _origin[0]);
            else if constexpr(dim_ == 2)
              grad = (DataType(-1) / (_ln2 * norm_sqr)) * (point - _origin);
            else if constexpr(dim_ == 3)
              grad = (DataType(1) / Math::pow(norm_sqr, DataType(1.5))) * (_origin - point);

            return grad;
          }

          HessianType hessian(const PointType& point) const
          {
            const DataType norm_sqr = (point - _origin).norm_euclid_sqr();
            if(norm_sqr < _eps)
              return HessianType::null();

            HessianType hess(DataType(0));
            if constexpr(dim_ == 1)
            {
              // 1D hessian is null
            }
            else if constexpr(dim_ == 2)
            {
              const DataType denom = DataType(1) / (_ln2 * Math::sqr(norm_sqr));
              hess[0][0] = denom * (Math::sqr(point[0] - _origin[0]) - Math::sqr(point[1] - _origin[1]));
              hess[1][1] = denom * (Math::sqr(point[1] - _origin[1]) - Math::sqr(point[0] - _origin[0]));
              hess[0][1] = hess[1][0] = denom * DataType(2) * (point[0] - _origin[0]) * (point[1] - _origin[1]);
            }
            else if constexpr(dim_ == 3)
            {
              const DataType denom = DataType(1) / Math::pow(norm_sqr, DataType(2.5));
              for(int i(0); i < dim_; ++i)
              {
                for(int j(0); j < dim_; ++j)
                {
                  if(i == j)
                  {
                    for(int k(0); k < dim_; ++k)
                      hess[i][j] += DataType(k == i ? 2 : -1) * Math::sqr(point[k] - _origin[k]) * denom;
                  }
                  else
                  {
                    hess[i][j] = DataType(3) * (point[i] - _origin[i]) * (point[j] - _origin[j]) * denom;
                  }
                }
              }
            }
            return hess;
          }
        }; // class HarmonicShellFunction::Evaluator<...>

      public:
        /// Point to calculate the distance from
        PointType _origin;

      public:
        /// default constructor
        HarmonicShellFunction()
        {
          _origin.format();
        }

        /// Constructor
        explicit HarmonicShellFunction(const PointType& origin) :
          _origin(origin)
        {
        }

        /// Sets _point to x0_
        void set_point(const PointType& origin)
        {
          _origin = origin;
        }
      }; // class HarmonicShellFunction
    } // namespace Common
  } // namespace Analytic
} // namespace FEAT

#endif // KERNEL_ANALYTIC_COMMON_HPP
