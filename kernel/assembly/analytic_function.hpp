#pragma once
#ifndef KERNEL_ASSEMBLY_ANALYTIC_FUNCTION_HPP
#define KERNEL_ASSEMBLY_ANALYTIC_FUNCTION_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/trafo/base.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Analytic function interface class
     *
     * This class acts as a base class and interface documentation for basic analytic functors, which are used
     * by various assemblies.
     *
     * \author Peter Zajac
     */
    class AnalyticFunction
    {
    public:
      /**
       * \brief Function Value capability
       *
       * This entry specifies whether the functor is capable of computing function values.\n
       * If this value is non-zero, the Evaluator class template offers the value() function.
       */
      static constexpr bool can_value = false;

      /**
       * \brief Gradient capability
       *
       * This entry specifies whether the functor is capable of computing gradients.\n
       * If this value is non-zero, the Evaluator class template offers the gradient() function.
       */
      static constexpr bool can_grad = false;

      /**
       * \brief Hessian capability
       *
       * This entry specifies whether the functor is capable of computing hessians.\n
       * If this value is non-zero, the Evaluator class template offers the hessian() function.
       */
      static constexpr bool can_hess = false;

      /**
       * \brief Configuration Traits
       *
       * \tparam Config_
       * A function configuration tag class. See Trafo::AnalyticConfigBase for details.
       */
      template<typename Config_>
      struct ConfigTraits
      {
        /**
         * \brief Trafo configuration tag class
         *
         * \see Trafo::ConfigBase
         */
        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          /// Let's assume that the function needs image point coordinates
          static constexpr bool need_img_point = true;
        };
      };

      /**
       * \brief Analytic function evaluator class template
       *
       * \tparam EvalTraits_
       * Specifies the traits for the evaluation.
       */
      template<typename EvalTraits_>
      class Evaluator
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

        /**
         * \brief Prepares the evaluator for a given cell
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& DOXY(trafo_eval))
        {
          // do nothing
        }

        /**
         * \brief Releases the evaluator from the current cell.
         */
        void finish()
        {
          // do nothing
        }

#ifdef DOXYGEN
        /**
         * \brief Value evaluation operator
         *
         * This operator evaluates the analytic function for a given point.
         *
         * \param[in] tau
         * The transformation data in the current evaluation point.
         *
         * \returns
         * The value of the function.
         */
        ValueType value(const TrafoData& tau) const;

        /**
         * \brief Gradient evaluation operator
         *
         * This operator evaluates the analytic function's gradient for a given point.
         *
         * \param[in] tau
         * The transformation data in the current evaluation point.
         *
         * \returns
         * The gradient of the function.
         */
        GradientType gradient(const TrafoData& tau) const;

        /**
         * \brief Hessian evaluation operator
         *
         * This operator evaluates the analytic function's hessian for a given point.
         *
         * \param[in] tau
         * The transformation data in the current evaluation point.
         *
         * \returns
         * The hessian of the function.
         */
        HessianType hessian(const TrafoData& tau) const;
#endif // DOXYGEN
      }; // class Functor::Evaluator<...>
    }; // class Functor

    /**
     * \brief Analytic static function interface class template
     *
     * \tparam DataType_
     * The data-type that is to be used for the evaluation of the function.
     *
     * \author Peter Zajac
     */
    template<typename DataType_>
    class StaticFunction
    {
    public:
#ifdef DOXYGEN
      /**
       * \brief Evaluates the function value.
       *
       * \param[in] x
       * The coordinates of the evaluation point.
       *
       * \returns
       * The function value
       */
      static DataType_ eval(DataType_ x);

      /**
       * \brief Evaluates the function value.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The function value
       */
      static DataType_ eval(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the function value.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The function value
       */
      static DataType_ eval(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the first order x-derivative.
       *
       * \param[in] x
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_x(DataType_ x);

      /**
       * \brief Evaluates the first order x-derivative.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_x(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the first order x-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the first order y-derivative.
       *
       * \param[in] x
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_y(DataType_ x);

      /**
       * \brief Evaluates the first order y-derivative.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_y(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the first order y-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the first order z-derivative.
       *
       * \param[in] x
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_z(DataType_ x);

      /**
       * \brief Evaluates the first order z-derivative.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_z(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the first order z-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order xx-derivative.
       *
       * \param[in] x
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_xx(DataType_ x);

      /**
       * \brief Evaluates the second order xx-derivative.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_xx(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the second order xx-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_xx(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order yy-derivative.
       *
       * \param[in] x
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_yy(DataType_ x);

      /**
       * \brief Evaluates the second order yy-derivative.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_yy(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the second order yy-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_yy(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order zz-derivative.
       *
       * \param[in] x
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_zz(DataType_ x);

      /**
       * \brief Evaluates the second order zz-derivative.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_zz(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the second order zz-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_zz(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order xy-derivative.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_xy(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the second order xy-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_xy(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order yx-derivative.
       *
       * \param[in] x,y
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_yx(DataType_ x, DataType_ y);

      /**
       * \brief Evaluates the second order yx-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_yx(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order xz-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_xz(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order zx-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_zx(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order yz-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_yz(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the second order zy-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       *
       * \returns
       * The value of the derivative
       */
      static DataType_ der_zy(DataType_ x, DataType_ y, DataType_ z);
#endif // DOXYGEN
    }; // class StaticFunction

    /// \cond internal
    namespace Intern
    {
      template<
        template<typename> class Function_,
        typename Traits_,
        bool enable_,
        int dim_ = Traits_::image_dim>
      struct StaticFunctionWrapper;

      // specialisation for 1D functions
      template<template<typename> class Function_, typename Traits_>
      struct StaticFunctionWrapper<Function_, Traits_, true, 1>
      {
        static typename Traits_::ValueType eval(const typename Traits_::TrafoData& tau)
        {
          return Function_<typename Traits_::DataType>::eval(tau.img_point[0]);
        }

        static typename Traits_::GradientType grad(const typename Traits_::TrafoData& tau)
        {
          typename Traits_::GradientType value;
          value(0) = Function_<typename Traits_::DataType>::der_x(tau.img_point[0]);
          return value;
        }

        static typename Traits_::HessianType hess(const typename Traits_::TrafoData& tau)
        {
          typename Traits_::HessianType value;
          value(0,0) = Function_<typename Traits_::DataType>::der_xx(tau.img_point[0]);
          return value;
        }
      };

      // specialisation for 2D functions
      template<template<typename> class Function_, typename Traits_>
      struct StaticFunctionWrapper<Function_, Traits_, true, 2>
      {
        static typename Traits_::ValueType eval(const typename Traits_::TrafoData& tau)
        {
          typedef Function_<typename Traits_::DataType> FuncType;
          return FuncType::eval(tau.img_point[0], tau.img_point[1]);
        }

        static typename Traits_::GradientType grad(const typename Traits_::TrafoData& tau)
        {
          typedef Function_<typename Traits_::DataType> FuncType;
          typename Traits_::GradientType value;
          value(0) = FuncType::der_x(tau.img_point[0], tau.img_point[1]);
          value(1) = FuncType::der_y(tau.img_point[0], tau.img_point[1]);
          return value;
        }

        static typename Traits_::HessianType hess(const typename Traits_::TrafoData& tau)
        {
          typedef Function_<typename Traits_::DataType> FuncType;
          typename Traits_::HessianType value;
          value(0,0) = FuncType::der_xx(tau.img_point[0], tau.img_point[1]);
          value(0,1) = FuncType::der_xy(tau.img_point[0], tau.img_point[1]);
          value(1,0) = FuncType::der_yx(tau.img_point[0], tau.img_point[1]);
          value(1,1) = FuncType::der_yy(tau.img_point[0], tau.img_point[1]);
          return value;
        }
      };

      // specialisation for 3D functions
      template<template<typename> class Function_, typename Traits_>
      struct StaticFunctionWrapper<Function_, Traits_, true, 3>
      {
        static typename Traits_::ValueType eval(const typename Traits_::TrafoData& tau)
        {
          typedef Function_<typename Traits_::DataType> FuncType;
          return FuncType::eval(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
        }

        static typename Traits_::GradientType grad(const typename Traits_::TrafoData& tau)
        {
          typedef Function_<typename Traits_::DataType> FuncType;
          typename Traits_::GradientType value;
          value(0) = FuncType::der_x(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(1) = FuncType::der_y(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(2) = FuncType::der_z(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          return value;
        }

        static typename Traits_::HessianType hess(const typename Traits_::TrafoData& tau)
        {
          typedef Function_<typename Traits_::DataType> FuncType;
          typename Traits_::HessianType value;
          value(0,0) = FuncType::der_xx(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(0,1) = FuncType::der_xy(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(0,2) = FuncType::der_xz(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(1,0) = FuncType::der_yx(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(1,1) = FuncType::der_yy(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(1,2) = FuncType::der_yz(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(2,0) = FuncType::der_zx(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(2,1) = FuncType::der_zy(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value(2,2) = FuncType::der_zz(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          return value;
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief StaticFunction wrapper class template for AnalyticFunction interface
     *
     * This class template implements the AnalyticFunction interface by wrapping around a StaticFunction
     * class template.
     *
     * \tparam Function_
     * A class template implementing the StaticFunction interface that is to be wrapped.
     *
     * \tparam can_value_
     * Specifies whether the static function class template supports the evaluation of function values.
     *
     * \tparam can_grad_
     * Specifies whether the static function class template supports the evaluation of first order derivatives.
     *
     * \tparam can_hess_
     * Specifies whether the static function class template supports the evaluation of second order derivatives.
     *
     * \author Peter Zajac
     */
    template<
      template<typename> class Function_,
      bool can_value_ = true,
      bool can_grad_ = false,
      bool can_hess_ = false>
    class StaticWrapperFunction :
      public AnalyticFunction
    {
    public:
      static constexpr bool can_value = can_value_;
      static constexpr bool can_grad = can_grad_;
      static constexpr bool can_hess = can_hess_;

      /** \copydoc AnalyticFunction::ConfigTraits */
      template<typename Config_>
      struct ConfigTraits
      {
        // ensure that we have everything available
        static_assert(can_value_ || (!Config_::need_value), "static function can't compute function values");
        static_assert(can_grad_  || (!Config_::need_grad),  "static function can't compute function gradients");
        static_assert(can_hess_  || (!Config_::need_hess),  "static function can't compute function hessians");

        /**
         * \brief Trafo configuration tag class
         *
         * \see Trafo::ConfigBase
         */
        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          /// Let's assume that the function needs image point coordinates
          static constexpr bool need_img_point = true;
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

      public:
        /// Constructor
        explicit Evaluator(const StaticWrapperFunction&)
        {
        }

        ValueType value(const TrafoData& tau) const
        {
          return Intern::StaticFunctionWrapper<Function_, EvalTraits_, can_value_>::eval(tau);
        }

        GradientType gradient(const TrafoData& tau) const
        {
          return Intern::StaticFunctionWrapper<Function_, EvalTraits_, can_grad_>::grad(tau);
        }

        HessianType hessian(const TrafoData& tau) const
        {
          return Intern::StaticFunctionWrapper<Function_, EvalTraits_, can_hess_>::hess(tau);
        }
      }; // class StaticWrapperFunction::Evaluator<...>
    }; // class StaticWrapperFunction
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_ANALYTIC_FUNCTION_HPP
