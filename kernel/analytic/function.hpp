// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ANALYTIC_FUNCTION_HPP
#define KERNEL_ANALYTIC_FUNCTION_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  /**
   * \brief Analytic namespace
   */
  namespace Analytic
  {
    /// Analytic Image namespace
    namespace Image
    {
      /**
       * \brief Scalar Function Image tag class
       *
       * This tag class represents a scalar function.
       */
      struct Scalar
      {
        /// this is a scalar function
        static constexpr bool is_scalar = true;
        /// this is not a vector field
        static constexpr bool is_vector = false;
        /// this function has 1 scalar component
        static constexpr int scalar_components = 1;
      };

      /**
       * \brief Vector Field Image tag class
       *
       * This tag class represents a vector field.
       *
       * \tparam m_
       * The image dimension of the vector field. Must be > 0.
       */
      template<int m_>
      struct Vector
      {
        static_assert(m_ > 0, "invalid vector field dimension");

        /// this is not a scalar function
        static constexpr bool is_scalar = false;
        /// this is a vector field
        static constexpr bool is_vector = true;
        /// this function has m_ scalar components
        static constexpr int scalar_components = m_;

        /// this is the image dimension of the vector field
        static constexpr int image_dim = m_;
      };
    } // namespace Image

    /// analytic evaluation traits base-class
    template<typename DataType_, int domain_dim_, typename ImageType_>
    struct EvalTraitsBase;

    /// analytic evaluation traits for scalar functions
    template<typename DataType_, int domain_dim_>
    struct EvalTraitsBase<DataType_, domain_dim_, Image::Scalar>
    {
      typedef DataType_ DataType;
      static constexpr int domain_dim = domain_dim_;

      typedef Tiny::Vector<DataType_, domain_dim_> PointType;

      typedef DataType_ ValueType;
      typedef Tiny::Vector<DataType_, domain_dim_> GradientType;
      typedef Tiny::Matrix<DataType_, domain_dim_, domain_dim_> HessianType;
    };

    /// analytic evaluation traits for vector-valued functions
    template<typename DataType_, int domain_dim_, int image_dim_>
    struct EvalTraitsBase<DataType_, domain_dim_, Image::Vector<image_dim_>>
    {
      typedef DataType_ DataType;
      static constexpr int domain_dim = domain_dim_;
      static constexpr int image_dim = image_dim_;

      typedef Tiny::Vector<DataType_, domain_dim_> PointType;

      typedef Tiny::Vector<DataType_, image_dim_> ValueType;
      typedef Tiny::Matrix<DataType_, image_dim_, domain_dim_> GradientType;
      typedef Tiny::Tensor3<DataType_, image_dim_, domain_dim_, domain_dim_> HessianType;
    };

    template<typename DataType_, typename Function_>
    struct EvalTraits :
      public EvalTraitsBase<DataType_, Function_::domain_dim, typename Function_::ImageType>
    {
    };

    /**
     * \brief Analytic Function interface
     *
     * This class acts as a base-class and interface description for analytic functions,
     * which can be used for various assembly purposes such as right-hand-side and boundary
     * condition assembly as well as in post-processing.
     *
     * \author Peter Zajac
     */
    class Function
    {
    public:
      // Note:
      // The following block is only visible for doxygen. Its contents must be specified
      // by the derived classes.
#ifdef DOXYGEN
      /**
       * \brief Specifies the domain dimension of the function.
       *
       * \attention This member must be specified by each derived class.
       */
      static constexpr int domain_dim = ...;

      /**
       * \brief Specifies the image type of the function.
       *
       * This is a typedef for either Image::Scalar or an instance of Image::Vector, which
       * specifies whether this is a scalar function or a vector field.
       *
       * \attention This member must be specified by each derived class.
       */
      typedef ... ImageType;

      /// Specifies whether the function's evaluator can compute function values.
      static constexpr bool can_value = true | false;
      /// Specifies whether the function's evaluator can compute function gradients.
      static constexpr bool can_grad = true | false;
      /// Specifies whether the function's evaluator can compute function hessians.
      static constexpr bool can_hess = true | false;
#else
      static constexpr bool can_value = false;
      static constexpr bool can_grad = false;
      static constexpr bool can_hess = false;
#endif // DOXYGEN

      /**
       * \brief Analytic Function Evaluator base-class template
       *
       * \tparam Traits_
       * A traits class that contains various typedefs which specify the data types
       * for the evaluation. See Analytic::EvalTraitsBase and its specializations for
       * all contained types.
       */
      template<typename Traits_>
      class Evaluator
      {
      public:
        /**
         * \brief The underlying floating point data type.
         *
         * This is the fundamental scalar datatype which is used by all other typedefs
         * in the traits class.
         */
        typedef typename Traits_::DataType DataType;

        /**
         * \brief The type of the domain evaluation point.
         *
         * This type specifies the type of the point in which the function is to be
         * evaluated. This is always a Tiny::Vector<DataType, domain_dim>.
         */
        typedef typename Traits_::PointType PointType;

        /**
         * \brief The type of the function value.
         *
         * This type specifies the type of the function's value:
         * - If the function is scalar, then ValueType is the same as DataType.
         * - If the function is a vector field, then ValueType is Tiny::Vector<DataType, image_dim>,
         *   where image_dim is the dimension of the vector field.
         */
        typedef typename Traits_::ValueType ValueType;

        /**
         * \brief The type of the function gradient.
         *
         * This type specifies the type of the function's gradient:
         * - If the function is scalar, then GradientType is Tiny::Vector<DataType, domain_dim>.
         * - If the function is a vector field, the GradientType is Tiny::Matrix<DataType, image_dim, domain_dim>.
         */
        typedef typename Traits_::GradientType GradientType;

        /**
         * \brief The type of the function hessian.
         *
         * This type specifies the type of the function's hessian:
         * - If the function is scalar, then HessianType is Tiny::Matrix<DataType, domain_dim, domain_dim>.
         * - If the function is a vector field, the HessianType is Tiny::Tensor3<DataType, image_dim, domain_dim, domain_dim>.
         */
        typedef typename Traits_::HessianType HessianType;

      public:
#ifdef DOXYGEN
        /**
         * \brief Mandatory constructor
         *
         * \param[in] function
         * A const-reference to the function object.
         */
        explicit Evaluator(const Function& function)
        {
        }

        /**
         * \brief Computes the function value in a given point.
         *
         * \param[in] point
         * The point in which the function is to be evaluated.
         *
         * \returns
         * The function value.
         */
        ValueType value(const PointType& point)
        {
        }

        /**
         * \brief Computes the function gradient in a given point.
         *
         * \param[in] point
         * The point in which the function is to be evaluated.
         *
         * \returns
         * The function gradient.
         */
        GradientType gradient(const PointType& point)
        {
        }

        /**
         * \brief Computes the function hessian in a given point.
         *
         * \param[in] point
         * The point in which the function is to be evaluated.
         *
         * \returns
         * The function hessian.
         */
        HessianType hessian(const PointType& point)
        {
        }
#endif // DOXYGEN
      }; // class Function::Evaluator<...>
    }; // class Function

    /**
     * \brief Helper function to quickly evaluate a function value in a given point
     *
     * \attention
     * This function is intended only for quick testing and debugging, because it recreates
     * the function evaluator in each function call, which is an absolute performance killer.\n
     * You have been warned.
     *
     * \param[in] function
     * The function whose value is to be evaluated.
     *
     * \param[in] point
     * The point in which the function is to be evaluated.
     *
     * \returns The function value of the given function in the given point.
     */
    template<typename Function_, typename DT_, int dim_, int s_>
    typename Analytic::EvalTraits<DT_, Function_>::ValueType
    eval_value(const Function_& function, const Tiny::Vector<DT_, dim_, s_>& point)
    {
      typedef Analytic::EvalTraits<DT_, Function_> Traits;
      static_assert(Function_::can_value, "function does not support evaluation of values");
      static_assert(dim_ == Traits::domain_dim, "invalid point dimension");
      typename Function_::template Evaluator<Traits> evaluator(function);
      return evaluator.value(point);
    }

    // \see eval_value()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::ValueType
    eval_value_x(const Function_& function, const DT_ x)
    {
      Tiny::Vector<DT_, 1> p;
      p.v[0] = x;
      return eval_value(function, p);
    }

    // \see eval_value()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::ValueType
    eval_value_x(const Function_& function, const DT_ x, const DT_ y)
    {
      Tiny::Vector<DT_, 2> p;
      p.v[0] = x;
      p.v[1] = y;
      return eval_value(function, p);
    }

    // \see eval_value()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::ValueType
    eval_value_x(const Function_& function, const DT_ x, const DT_ y, const DT_ z)
    {
      Tiny::Vector<DT_, 3> p;
      p.v[0] = x;
      p.v[1] = y;
      p.v[2] = z;
      return eval_value(function, p);
    }

    /**
     * \brief Helper function to quickly evaluate a function gradient in a given point
     *
     * \attention
     * This function is intended only for quick testing and debugging, because it recreates
     * the function evaluator in each function call, which is an absolute performance killer.\n
     * You have been warned.
     *
     * \param[in] function
     * The function whose gradient is to be evaluated.
     *
     * \param[in] point
     * The point in which the function is to be evaluated.
     *
     * \returns The function gradient of the given function in the given point.
     */
    template<typename Function_, typename DT_, int dim_, int s_>
    typename Analytic::EvalTraits<DT_, Function_>::GradientType
    eval_gradient(const Function_& function, const Tiny::Vector<DT_, dim_, s_>& point)
    {
      typedef Analytic::EvalTraits<DT_, Function_> Traits;
      static_assert(Function_::can_grad, "function does not support evaluation of gradients");
      static_assert(dim_ == Traits::domain_dim, "invalid point dimension");
      typename Function_::template Evaluator<Traits> evaluator(function);
      return evaluator.gradient(point);
    }

    // \see eval_gradient()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::GradientType
    eval_gradient_x(const Function_& function, const DT_ x)
    {
      Tiny::Vector<DT_, 1> p;
      p.v[0] = x;
      return eval_gradient(function, p);
    }

    // \see eval_gradient()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::GradientType
    eval_gradient_x(const Function_& function, const DT_ x, const DT_ y)
    {
      Tiny::Vector<DT_, 2> p;
      p.v[0] = x;
      p.v[1] = y;
      return eval_gradient(function, p);
    }

    // \see eval_gradient()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::GradientType
    eval_gradient_x(const Function_& function, const DT_ x, const DT_ y, const DT_ z)
    {
      Tiny::Vector<DT_, 3> p;
      p.v[0] = x;
      p.v[1] = y;
      p.v[2] = z;
      return eval_gradient(function, p);
    }

    /**
     * \brief Helper function to quickly evaluate a function hessian in a given point
     *
     * \attention
     * This function is intended only for quick testing and debugging, because it recreates
     * the function evaluator in each function call, which is an absolute performance killer.\n
     * You have been warned.
     *
     * \param[in] function
     * The function whose hessian is to be evaluated.
     *
     * \param[in] point
     * The point in which the function is to be evaluated.
     *
     * \returns The function hessian of the given function in the given point.
     */
    template<typename Function_, typename DT_, int dim_, int s_>
    typename Analytic::EvalTraits<DT_, Function_>::HessianType
    eval_hessian(const Function_& function, const Tiny::Vector<DT_, dim_, s_>& point)
    {
      typedef Analytic::EvalTraits<DT_, Function_> Traits;
      static_assert(Function_::can_hess, "function does not support evaluation of hessians");
      static_assert(dim_ == Traits::domain_dim, "invalid point dimension");
      typename Function_::template Evaluator<Traits> evaluator(function);
      return evaluator.hessian(point);
    }

    // \see eval_hessian()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::HessianType
    eval_hessian_x(const Function_& function, const DT_ x)
    {
      Tiny::Vector<DT_, 1> p;
      p.v[0] = x;
      return eval_hessian(function, p);
    }

    // \see eval_hessian()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::HessianType
    eval_hessian_x(const Function_& function, const DT_ x, const DT_ y)
    {
      Tiny::Vector<DT_, 2> p;
      p.v[0] = x;
      p.v[1] = y;
      return eval_hessian(function, p);
    }

    // \see eval_hessian()
    template<typename Function_, typename DT_>
    typename Analytic::EvalTraits<DT_, Function_>::HessianType
    eval_hessian_x(const Function_& function, const DT_ x, const DT_ y, const DT_ z)
    {
      Tiny::Vector<DT_, 3> p;
      p.v[0] = x;
      p.v[1] = y;
      p.v[2] = z;
      return eval_hessian(function, p);
    }

  } // namespace Analytic
} // namespace FEAT

#endif // KERNEL_ANALYTIC_FUNCTION_HPP
