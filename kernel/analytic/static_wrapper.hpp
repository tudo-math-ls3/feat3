// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ANALYTIC_STATIC_WRAPPER_FUNCTION_HPP
#define KERNEL_ANALYTIC_STATIC_WRAPPER_FUNCTION_HPP 1

// includes, FEAT
#include <kernel/analytic/function.hpp>

namespace FEAT
{
  namespace Analytic
  {
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
        typename DataType_,
        int dim_,
        bool enable_>
      struct StaticFunctionWrapper;

      // specialisation for 1D functions
      template<template<typename> class Function_, typename DataType_>
      struct StaticFunctionWrapper<Function_, DataType_, 1, true>
      {
        template<typename Value_, typename Point_>
        static void eval(Value_& v, const Point_& x)
        {
          v = Function_<DataType_>::eval(x[0]);
        }

        template<typename Value_, typename Point_>
        static void grad(Value_& v, const Point_& x)
        {
          v[0] = Function_<DataType_>::der_x(x[0]);
        }

        template<typename Value_, typename Point_>
        static void hess(Value_& v, const Point_& x)
        {
          v[0][0] = Function_<DataType_>::der_xx(x[0]);
        }
      };

      // specialisation for 2D functions
      template<template<typename> class Function_, typename DataType_>
      struct StaticFunctionWrapper<Function_, DataType_, 2, true>
      {
        template<typename Value_, typename Point_>
        static void eval(Value_& v, const Point_& x)
        {
          v = Function_<DataType_>::eval(x[0], x[1]);
        }

        template<typename Value_, typename Point_>
        static void grad(Value_& v, const Point_& x)
        {
          v[0] = Function_<DataType_>::der_x(x[0], x[1]);
          v[1] = Function_<DataType_>::der_y(x[0], x[1]);
        }

        template<typename Value_, typename Point_>
        static void hess(Value_& v, const Point_& x)
        {
          v[0][0] = Function_<DataType_>::der_xx(x[0], x[1]);
          v[0][1] = Function_<DataType_>::der_xy(x[0], x[1]);
          v[1][0] = Function_<DataType_>::der_yx(x[0], x[1]);
          v[1][1] = Function_<DataType_>::der_yy(x[0], x[1]);
        }
      };

      // specialisation for 3D functions
      template<template<typename> class Function_, typename DataType_>
      struct StaticFunctionWrapper<Function_, DataType_, 3, true>
      {
        template<typename Value_, typename Point_>
        static void eval(Value_& v, const Point_& x)
        {
          v = Function_<DataType_>::eval(x[0], x[1],x[2]);
        }

        template<typename Value_, typename Point_>
        static void grad(Value_& v, const Point_& x)
        {
          v[0] = Function_<DataType_>::der_x(x[0], x[1], x[2]);
          v[1] = Function_<DataType_>::der_y(x[0], x[1], x[2]);
          v[2] = Function_<DataType_>::der_z(x[0], x[1], x[2]);
        }

        template<typename Value_, typename Point_>
        static void hess(Value_& v, const Point_& x)
        {
          v[0][0] = Function_<DataType_>::der_xx(x[0], x[1], x[2]);
          v[0][1] = Function_<DataType_>::der_xy(x[0], x[1], x[2]);
          v[0][2] = Function_<DataType_>::der_xz(x[0], x[1], x[2]);
          v[1][0] = Function_<DataType_>::der_yx(x[0], x[1], x[2]);
          v[1][1] = Function_<DataType_>::der_yy(x[0], x[1], x[2]);
          v[1][2] = Function_<DataType_>::der_yz(x[0], x[1], x[2]);
          v[2][0] = Function_<DataType_>::der_zx(x[0], x[1], x[2]);
          v[2][1] = Function_<DataType_>::der_zy(x[0], x[1], x[2]);
          v[2][2] = Function_<DataType_>::der_zz(x[0], x[1], x[2]);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief StaticFunction wrapper class template for Analytic::Function interface
     *
     * This class template implements the Analytic::Function interface by wrapping around a StaticFunction
     * class template.
     *
     * \tparam domain_dim_
     * The domain dimension of the function. Must be 1 <= domain_dim_ <= 3.
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
      int domain_dim_,
      template<typename> class Function_,
      bool can_value_ = true,
      bool can_grad_ = false,
      bool can_hess_ = false>
    class StaticWrapperFunction :
      public Analytic::Function
    {
    public:
      static_assert((1 <= domain_dim_) && (domain_dim_ <= 3), "invalid domain dimension");

      /// our domain dimension
      static constexpr int domain_dim = domain_dim_;

      /// this is always a scalar function
      typedef Analytic::Image::Scalar ImageType;

      static constexpr bool can_value = can_value_;
      static constexpr bool can_grad = can_grad_;
      static constexpr bool can_hess = can_hess_;

      StaticWrapperFunction() {}
      StaticWrapperFunction(const StaticWrapperFunction&) {}
      StaticWrapperFunction(StaticWrapperFunction&&) {}

      /** \copydoc Analytic::Function::Evaluator */
      template<typename Traits_>
      class Evaluator :
        public Analytic::Function::Evaluator<Traits_>
      {
      public:
        /// coefficient data type
        typedef typename Traits_::DataType DataType;
        /// domain point type
        typedef typename Traits_::PointType PointType;
        /// value type
        typedef typename Traits_::ValueType ValueType;
        /// gradient type
        typedef typename Traits_::GradientType GradientType;
        /// hessian type
        typedef typename Traits_::HessianType HessianType;

      public:
        /// Constructor
        explicit Evaluator(const StaticWrapperFunction&)
        {
        }

        ValueType value(const PointType& point) const
        {
          ValueType val;
          Intern::StaticFunctionWrapper<Function_, DataType, domain_dim, can_value_>::eval(val, point);
          return val;
        }

        GradientType gradient(const PointType& point) const
        {
          GradientType grad;
          Intern::StaticFunctionWrapper<Function_, DataType, domain_dim, can_grad_>::grad(grad, point);
          return grad;
        }

        HessianType hessian(const PointType& point) const
        {
          HessianType hess;
          Intern::StaticFunctionWrapper<Function_, DataType, domain_dim, can_hess_>::hess(hess, point);
          return hess;
        }
      }; // class StaticWrapperFunction::Evaluator<...>
    }; // class StaticWrapperFunction
  } // namespace Analytic
} // namespace FEAT

#endif // KERNEL_ANALYTIC_STATIC_WRAPPER_FUNCTION_HPP
