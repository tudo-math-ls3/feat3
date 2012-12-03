#pragma once
#ifndef KERNEL_ANALYTIC_FUNCTOR_HPP
#define KERNEL_ANALYTIC_FUNCTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>

namespace FEAST
{
  /**
   * \brief Analytic namespace
   *
   * This namespace defines functor interfaces for (user-defined) analytic functions.
   */
  namespace Analytic
  {
    /**
     * \brief Analytic functor interface class
     *
     * This class acts as a base class and interface documentation for basic analytic functors, which are used
     * by various assemblies.
     *
     * \author Peter Zajac
     */
    class Functor
    {
    public:
      /**
       * \brief Functor capabilities enumeration
       */
      enum FunctorCapabilities
      {
        /**
         * \brief Function Value capability
         *
         * This entry specifies whether the functor is capable of computing function values.\n
         * If this value is non-zero, the functor offers a ValueEvaluator class template.
         */
        can_value = 0,

        /**
         * \brief Gradient capability
         *
         * This entry specifies whether the functor is capable of computing gradients.\n
         * If this value is non-zero, the functor offers a GradientEvaluator class template.
         */
        can_grad = 0
      };

      /**
       * \brief Analytic functor value evaluator class template
       *
       * \tparam EvalTraits_
       * Specifies the traits for the evaluation.
       */
      template<typename EvalTraits_>
      class ValueEvaluator
      {
      public:
        /// trafo evaluator data
        typedef typename EvalTraits_::TrafoEvaluator TrafoEvaluator;
        /// trafo data type
        typedef typename EvalTraits_::TrafoData TrafoData;
        /// coefficient type
        typedef typename EvalTraits_::CoeffType CoeffType;
        /// value type
        typedef typename EvalTraits_::ValueType ValueType;

        /**
         * \brief Prepares the evaluator for a given cell
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& trafo_eval)
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
         * \brief Evaluation operator
         *
         * This operator evaluates the analytic functor for a given point.
         *
         * \param[out] value
         * A reference to the value that is to be computed.
         *
         * \param[in] tau
         * The transformation data in the current evaluation point.
         */
        void operator()(ValueType& value, const TrafoData& tau) const;
#endif // DOXYGEN
      }; // class Functor::ValueEvaluator<...>

      /**
       * \brief Analytic functor gradient evaluator class template
       *
       * \tparam EvalTraits_
       * Specifies the traits for the evaluation.
       */
      template<typename EvalTraits_>
      class GradientEvaluator
      {
      public:
        /// trafo evaluator data
        typedef typename EvalTraits_::TrafoEvaluator TrafoEvaluator;
        /// trafo data type
        typedef typename EvalTraits_::TrafoData TrafoData;
        /// coefficient type
        typedef typename EvalTraits_::CoeffType CoeffType;
        /// value type
        typedef typename EvalTraits_::ValueType ValueType;

        /**
         * \brief Prepares the evaluator for a given cell
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& trafo_eval)
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
         * \brief Evaluation operator
         *
         * This operator evaluates the analytic functor for a given point.
         *
         * \param[out] value
         * A reference to the value that is to be computed.
         *
         * \param[in] tau
         * The transformation data in the current evaluation point.
         */
        void operator()(ValueType& value, const TrafoData& tau) const;
#endif // DOXYGEN
      }; // class Functor::GradientEvaluator<...>
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
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       */
      static DataType_ eval(DataType_ x);
      /** \copydoc eval(DataType_) */
      static DataType_ eval(DataType_ x, DataType_ y);
      /** \copydoc eval(DataType_) */
      static DataType_ eval(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the first order X-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       */
      static DataType_ der_x(DataType_ x);
      /** \copydoc der_x(DataType_) */
      static DataType_ der_x(DataType_ x, DataType_ y);
      static DataType_ der_x(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the first order Y-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       */
      static DataType_ der_y(DataType_ x, DataType_ y);
      /** \copydoc der_y(DataType_, DataType_) */
      static DataType_ der_y(DataType_ x, DataType_ y, DataType_ z);

      /**
       * \brief Evaluates the first order Z-derivative.
       *
       * \param[in] x,y,z
       * The coordinates of the evaluation point.
       */
      static DataType_ der_z(DataType_ x, DataType_ y, DataType_ z);
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

      template<template<typename> class Function_, typename Traits_>
      struct StaticFunctionWrapper<Function_, Traits_, true, 1>
      {
        static void eval(typename Traits_::ValueType& value, const typename Traits_::TrafoData& tau)
        {
          value = Function_<typename Traits_::CoeffType>::eval(tau.img_point[0]);
        }

        static void grad(typename Traits_::ValueType& value, const typename Traits_::TrafoData& tau)
        {
          value[0] = Function_<typename Traits_::CoeffType>::der_x(tau.img_point[0]);
        }
      };

      template<template<typename> class Function_, typename Traits_>
      struct StaticFunctionWrapper<Function_, Traits_, true, 2>
      {
        static void eval(typename Traits_::ValueType& value, const typename Traits_::TrafoData& tau)
        {
          value = Function_<typename Traits_::CoeffType>
            ::eval(tau.img_point[0], tau.img_point[1]);
        }

        static void grad(typename Traits_::ValueType& value, const typename Traits_::TrafoData& tau)
        {
          value[0] = Function_<typename Traits_::CoeffType>
            ::der_x(tau.img_point[0], tau.img_point[1]);
          value[1] = Function_<typename Traits_::CoeffType>
            ::der_y(tau.img_point[0], tau.img_point[1]);
        }
      };

      template<template<typename> class Function_, typename Traits_>
      struct StaticFunctionWrapper<Function_, Traits_, true, 3>
      {
        static void eval(typename Traits_::ValueType& value, const typename Traits_::TrafoData& tau)
        {
          value = Function_<typename Traits_::CoeffType>
            ::eval(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
        }

        static void grad(typename Traits_::ValueType& value, const typename Traits_::TrafoData& tau)
        {
          value[0] = Function_<typename Traits_::CoeffType>
            ::der_x(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value[1] = Function_<typename Traits_::CoeffType>
            ::der_y(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
          value[2] = Function_<typename Traits_::CoeffType>
            ::der_z(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
        }
      };
    } // namespace Intern
    /// \cond internal

    /**
     * \brief StaticFunction wrapper class template for Functor interface
     *
     * This class template implements the Analytic::Functor interface by wrapping around a Analytic::StaticFunction
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
     * \author Peter Zajac
     */
    template<
      template<typename> class Function_,
      bool can_value_ = true,
      bool can_grad_ = false>
    class StaticWrapperFunctor :
      public Functor
    {
    public:
      /** \copydoc Functor::FunctorCapabilities */
      enum
      {
        can_value = can_value_ ? 1 : 0,
        can_grad = can_grad_ ? 1 : 0
      };

      /** \copydoc Functor::ValueEvaluator */
      template<typename EvalTraits_>
      class ValueEvaluator :
        public Functor::ValueEvaluator<EvalTraits_>
      {
      public:
        explicit ValueEvaluator(const StaticWrapperFunctor&)
        {
        }

        void operator()(typename EvalTraits_::ValueType& value, const typename EvalTraits_::TrafoData& tau) const
        {
          return Intern::StaticFunctionWrapper<Function_, EvalTraits_, can_value_>::eval(value, tau);
        }
      }; // class StaticWrapperFunctor::ValueEvaluator<...>

      /** \copydoc Functor::GradientEvaluator */
      template<typename EvalTraits_>
      class GradientEvaluator :
        public Functor::GradientEvaluator<EvalTraits_>
      {
      public:
        explicit GradientEvaluator(const StaticWrapperFunctor&)
        {
        }

        void operator()(typename EvalTraits_::ValueType& value, const typename EvalTraits_::TrafoData& tau) const
        {
          return Intern::StaticFunctionWrapper<Function_, EvalTraits_, can_grad_>::grad(value, tau);
        }
      }; // class StaticWrapperFunctor::ValueEvaluator<...>
    }; // class StaticWrapperFunctor
  } // namespace Analytic
} // namespace FEAST

#endif // KERNEL_ANALYTIC_FUNCTOR_HPP
