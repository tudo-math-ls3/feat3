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
       * boundaries xi,yi,zi, this function fulfills homogene Dirichlet boundary conditions.
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
       * This function fulfills homogene Dirichlet boundary conditions on the unit-cube domain.
       *
       * \author Peter Zajac
       */
      typedef StaticWrapperFunction<SineBubbleStatic, true, true, true> SineBubbleFunction;

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
       * boundaries xi,yi,zi, this function fulfills homogene Neumann boundary conditions including
       * the integral-mean condition \f$ int_\Omega u = 0 \f$.
       *
       * \author Peter Zajac
       */
      template<typename DataType_, int k_ = 1>
      class CosineTensorStatic
      {
        static_assert(k_ > 0, "parameter k_ must be a positive integer");

      public:
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
          return der_xy(x, y);
        }

        /// 3D: XZ-derivative
        static DataType_ der_xz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(kpi()) * Math::sin(kpi() * x) * Math::cos(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: ZX-derivative
        static DataType_ der_zx(DataType_ x, DataType_ y, DataType_ z)
        {
          return der_xz(x, y);
        }

        /// 3D: YZ-derivative
        static DataType_ der_yz(DataType_ x, DataType_ y, DataType_ z)
        {
          return Math::sqr(kpi()) * Math::cos(kpi() * x) * Math::sin(kpi() * y) * Math::sin(kpi() * z);
        }

        /// 3D: ZY-derivative
        static DataType_ der_zy(DataType_ x, DataType_ y, DataType_ z)
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
       * This function fulfills homogene Neumann boundary conditions and has vanishing integral
       * mean on the unit-cube domain.
       *
       * \author Peter Zajac
       */
      typedef StaticWrapperFunction<CosineWaveStatic, true, true, true> CosineWaveFunction;

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
      template<typename ImgPointType_>
      class DistanceFunction :
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
          struct TrafoConfig :
            public Trafo::ConfigBase
          {
            enum
            {
              need_img_point = 1
            };
          };
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
          /// type for the points the analytic function is evaluated at
          typedef typename EvalTraits_::ImagePointType ImgPointType;

        private:
          const DistanceFunction& _function;

        public:
          explicit Evaluator(const DistanceFunction& function) :
            _function(function)
          {
          }

          ValueType value(const TrafoData& tau) const
          {
            ImgPointType tmp (tau.img_point - _function._point);
            return tmp.norm_euclid();
          }

          GradientType gradient(const TrafoData& tau) const
          {
            GradientType grad(DataType(0));
            DataType norm(value(tau));

            if(norm <= Math::eps<DataType>())
            {
              return grad;
            }
            else
            {
              for(Index d(0); d < TrafoData::image_dim; ++d)
                grad(d,tau.img_point[d] - _function._point(d));

              return grad/norm;
            }
          }

          HessianType hessian(const TrafoData& tau) const
          {
            HessianType hess(DataType(0));
            DataType norm(value(tau));
            if(norm <= Math::eps<DataType>())
            {
              return hess;
            }
            else
            {
              norm = DataType(1)/norm;
              DataType denom = Math::sqr(norm)*norm;

              for(Index i(0); i < TrafoData::image_dim; ++i)
              {
                hess(i,i, DataType(1)*norm);
                for(Index j(0); j < TrafoData::image_dim; ++i)
                  hess(i,j,
                  ((tau.img_point[i] - _function._point(i)) * (tau.img_point[j] - _function._point(j)) )*denom);
              }


              return hess;
            }
          }
        }; // class DistanceFunction::Evaluator<...>

      public:
        ImgPointType_ _point;

      public:
        /// Constructor
        explicit DistanceFunction(const ImgPointType_ x0_) :
          _point(x0_)
        {
        }

        /// Sets _point to x0_
        void set_point(const ImgPointType_ x0_)
        {
          _point = x0_;
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
      template<typename ImgPointType_>
      class DistanceFunctionSD :
        public AnalyticFunction
      {
      public:
        typedef typename ImgPointType_::DataType DataType;
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
          struct TrafoConfig :
            public Trafo::ConfigBase
          {
            enum
            {
              need_img_point = 1
            };
          };
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
          /// type for the points the analytic function is evaluated at
          typedef typename EvalTraits_::ImagePointType ImgPointType;

        private:
          const DistanceFunctionSD& _function;

        public:
          explicit Evaluator(const DistanceFunctionSD& function) :
            _function(function)
          {
          }

          ValueType value(const TrafoData& tau) const
          {
            ImgPointType tmp (tau.img_point - _function._point);
            return _function._a + _function._b*tmp.norm_euclid();
          }

          GradientType gradient(const TrafoData& tau) const
          {
            GradientType grad(DataType(0));
            DataType norm(value(tau));
            if(norm <= Math::eps<DataType>())
            {
              return grad;
            }
            else
            {

              for(Index d(0); d < TrafoData::image_dim; ++d)
                grad(d,tau.img_point[d] - _function._point(d));

              return _function._b*grad/norm;
            }
          }

          HessianType hessian(const TrafoData& tau) const
          {
            HessianType hess(DataType(0));
            DataType norm(value(tau));
            if(norm <= Math::eps<DataType>())
            {
              return hess;
            }
            else
            {
              norm = DataType(1)/norm;
              DataType denom = Math::sqr(_function._b)*Math::sqr(norm)*norm;

              for(Index i(0); i < TrafoData::image_dim; ++i)
              {
                hess(i,i, DataType(1)*norm);
                for(Index j(0); j < TrafoData::image_dim; ++i)
                  hess(i,j,
                  ((tau.img_point[i] - _function._point(i)) * (tau.img_point[j] - _function._point(j)) )*denom);
              }

              return hess;
            }
          }
        }; // class DistanceFunctionSD::Evaluator<...>

      private:
        /// The point to which the distance to is calculated
        ImgPointType_ _point;
        /// Displacement of the function
        DataType _a;
        /// Scaling factor
        DataType _b;

      public:
        /// Constructor
        explicit DistanceFunctionSD(const ImgPointType_ x0_, const DataType a_, const DataType b_) :
          _point(x0_),
          _a(a_),
          _b(b_)
        {
        }

        /// Sets _point to x0_
        void set_point(const ImgPointType_ x0_)
        {
          _point = x0_;
        }
      }; // class DistanceFunctionSD

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
      typedef StaticWrapperFunction<HeavisideStatic, true, false, false> HeavisideFunction;

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
      typedef StaticWrapperFunction<HeavisideRegStatic, true, true, true> HeavisideRegFunction;
    } // namespace Common
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_COMMON_FUNCTIONS_HPP
