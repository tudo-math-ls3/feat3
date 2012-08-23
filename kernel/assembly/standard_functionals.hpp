#pragma once
#ifndef KERNEL_ASSEMBLY_STANDARD_FUNCTIONALS_HPP
#define KERNEL_ASSEMBLY_STANDARD_FUNCTIONALS_HPP 1

// includes, FEAST
#include <kernel/assembly/linear_functor_base.hpp>
#include <kernel/assembly/linear_functional.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<
        template<typename> class Function_,
        typename Traits_,
        int dim_ = Traits_::image_dim>
      struct AnalyticFunctionWrapper;

      template<template<typename> class Function_, typename Traits_>
      struct AnalyticFunctionWrapper<Function_, Traits_, 1>
      {
        static typename Traits_::DataType eval(const typename Traits_::TrafoData& tau)
        {
          return Function_<typename Traits_::DataType>::eval(tau.img_point[0]);
        }
      };

      template<template<typename> class Function_, typename Traits_>
      struct AnalyticFunctionWrapper<Function_, Traits_, 2>
      {
        static typename Traits_::DataType eval(const typename Traits_::TrafoData& tau)
        {
          return Function_<typename Traits_::DataType>::eval(tau.img_point[0], tau.img_point[1]);
        }
      };

      template<template<typename> class Function_, typename Traits_>
      struct AnalyticFunctionWrapper<Function_, Traits_, 3>
      {
        static typename Traits_::DataType eval(const typename Traits_::TrafoData& tau)
        {
          return Function_<typename Traits_::DataType>::eval(tau.img_point[0], tau.img_point[1], tau.img_point[2]);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Linear scalar Integral functional functor implementation
     *
     * This functor implements the LinearFunctor interface with
     *   \f[ \ell(\varphi) := \int_\Omega f\cdot\varphi \f]
     * for a static analytic function \e f.
     *
     * \tparam Function_
     * A template class implementing the function \e f for the linearform.
     *
     * \author Peter Zajac
     */
    template<template<typename DataType_> class Function_>
    class LinearScalarIntegralFunctor :
      public LinearFunctorBase
    {
    public:
      /// trafo config tag
      struct TrafoConfig :
        public Trafo::ConfigBase
      {
        /// dummy enum
        enum
        {
          /// we need image points for the analytic function
          need_img_point = 1
        };
      };

      /// space config tag
      struct SpaceConfig :
        public Space::ConfigBase
      {
        /// dummy enum
        enum
        {
          /// we need basis function values
          need_value = 1
        };
      };

      /**
       * \brief Linear Functor Evaluator class template
       *
       * \tparam AsmTraits_
       * The assembly traits class.
       *
       * \author Peter Zajac
       */
      template<typename AsmTraits_>
      class Evaluator :
        public LinearFunctorBase::Evaluator<AsmTraits_>
      {
      public:
        /// data type
        typedef typename AsmTraits_::DataType DataType;
        /// trafo data type
        typedef typename AsmTraits_::TrafoData TrafoData;
        /// test-function data type
        typedef typename AsmTraits_::FuncData FuncData;

      public:
        /// constructor
        explicit Evaluator(const LinearScalarIntegralFunctor<Function_>&)
        {
        }

        /** \copydoc LinearFunctorBase::Evaluator::operator() */
        DataType operator()(const TrafoData& tau, const FuncData& psi) const
        {
          return Intern::AnalyticFunctionWrapper<Function_, AsmTraits_>::eval(tau) * psi.value;
        }
      }; // class LinearScalarIntegralFunctor::Evaluator<...>

      /**
       * \brief Assembles a vector using the analytic function.
       *
       * \param[in,out] vector
       * The vector that is to be assembled.
       *
       * \param[in] space
       * A reference to the finite-element space to be used.
       *
       * \param[in] cubature_name
       * A string containing the name of the cubature rule to be used for integration.
       */
      template<
        typename Vector_,
        typename Space_>
      static void assemble(Vector_& vector, const Space_& space, const String& cubature_name)
      {
        LinearScalarIntegralFunctor functor;
        LinearFunctional<Vector_, LinearScalarIntegralFunctor<Function_>, Space_>::
          assemble(vector, functor, space, cubature_name);
      }
    }; // class LinearScalarIntegralFunctor
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_STANDARD_FUNCTIONALS_HPP
