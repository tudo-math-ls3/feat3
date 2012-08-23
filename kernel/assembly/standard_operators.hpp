#pragma once
#ifndef KERNEL_ASSEMBLY_STANDARD_OPERATORS_HPP
#define KERNEL_ASSEMBLY_STANDARD_OPERATORS_HPP 1

// includes, FEAST
#include <kernel/assembly/bilinear_functor_base.hpp>
#include <kernel/assembly/bilinear_operator.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Bilinear scalar Laplace operator functor implementation
     *
     * This functor implements the weak formulation of the bilinear scalar Laplace operator, i.e.
     *   \f[ \nabla \varphi \cdot \nabla\psi \f]
     *
     * This functor can be used with the BilinearOperator assembly class template to assemble a
     * scalar Laplace/Stiffness matrix.
     *
     * \author Peter Zajac
     */
    class BilinearScalarLaplaceFunctor :
      public BilinearFunctorBase
    {
    public:
      /// test space configuration
      struct TestConfig :
        public Space::ConfigBase
      {
        /// test-space requirement enumeration
        enum
        {
          /// this functor requires test-function gradients
          need_grad = 1
        };
      };

      /// trial space configuration
      struct TrialConfig :
        public Space::ConfigBase
      {
        /// trial-space requirement enumeration
        enum
        {
          /// this functor requires trial-function gradients
          need_grad = 1
        };
      };

      /**
       * \brief Bilinear scalar Laplace evaluator class template
       *
       * \tparam AsmTraits_
       * The assembly traits class.
       *
       * \author Peter Zajac
       */
      template<typename AsmTraits_>
      class Evaluator :
          public BilinearFunctorBase::Evaluator<AsmTraits_>
      {
      public:
        /// the data type to be used
        typedef typename AsmTraits_::DataType DataType;
        /// the assembler's trafo data type
        typedef typename AsmTraits_::TrafoData TrafoData;
        /// the assembler's test-function data type
        typedef typename AsmTraits_::TestFuncData TestFuncData;
        /// the assembler's trial-function data type
        typedef typename AsmTraits_::TrialFuncData TrialFuncData;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] functor
         * A reference to the scalar Laplace functor.
         */
        explicit Evaluator(const BilinearScalarLaplaceFunctor& functor)
        {
        }

        /** \copydoc BilinearFunctorBase::Evaluator::operator() */
        DataType operator()(const TrafoData& /*tau*/, const TrialFuncData& phi, const TestFuncData& psi)
        {
          return dot(phi.grad, psi.grad);
        }
      }; // class BilinearScalarLaplaceFunctor::Evaluator<...>


    public:
      /**
       * \brief Assembles a stiffness matrix.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_name
       * A string containing the name of the cubature rule to be used for integration.
       */
      template<
        typename Matrix_,
        typename Space_>
      static void assemble(Matrix_& matrix, const Space_& space, const String& cubature_name)
      {
        BilinearScalarLaplaceFunctor functor;
        BilinearOperator<Matrix_, BilinearScalarLaplaceFunctor, Space_>::
          assemble(matrix, functor, space, cubature_name);
      }
    }; // class BilinearScalarLaplaceFunctor

    /**
     * \brief Bilinear scalar Identity operator functor implementation
     *
     * This functor implements the weak formulation of the bilinear scalar Identity operator, i.e.
     *   \f[ \varphi \cdot \psi \f]
     *
     * This functor can be used with the BilinearOperator assembler class template to assemble a
     * scalar mass matrix.
     *
     * \author Peter Zajac
     */
    class BilinearScalarIdentityFunctor :
      public BilinearFunctorBase
    {
    public:
      /// test space configuration
      struct TestConfig :
        public Space::ConfigBase
      {
        /// test-space requirement enumeration
        enum
        {
          /// this functor requires test-function values
          need_value = 1
        };
      };

      /// trial space configuration
      struct TrialConfig :
        public Space::ConfigBase
      {
        /// trial-space requirement enumeration
        enum
        {
          /// this functor requires trial-function values
          need_value = 1
        };
      };

      /**
       * \brief Bilinear scalar Identity evaluator class template
       *
       * \tparam AsmTraits_
       * The assembly traits class.
       *
       * \author Peter Zajac
       */
      template<typename AsmTraits_>
      class Evaluator :
          public BilinearFunctorBase::Evaluator<AsmTraits_>
      {
      public:
        /// the data type to be used
        typedef typename AsmTraits_::DataType DataType;
        /// the assembler's trafo data type
        typedef typename AsmTraits_::TrafoData TrafoData;
        /// the assembler's test-function data type
        typedef typename AsmTraits_::TestFuncData TestFuncData;
        /// the assembler's trial-function data type
        typedef typename AsmTraits_::TrialFuncData TrialFuncData;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] functor
         * A reference to the scalar Identity functor.
         */
        explicit Evaluator(const BilinearScalarIdentityFunctor& functor)
        {
        }

        /** \copydoc BilinearFunctorBase::Evaluator::operator() */
        DataType operator()(const TrafoData& /*tau*/, const TrialFuncData& phi, const TestFuncData& psi)
        {
          return phi.value * psi.value;
        }
      }; // class BilinearScalarIdentityFunctor::Evaluator<...>

    public:
      /**
       * \brief Assembles a mass matrix.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_name
       * A string containing the name of the cubature rule to be used for integration.
       */
      template<
        typename Matrix_,
        typename Space_>
      static void assemble(Matrix_& matrix, const Space_& space, const String& cubature_name)
      {
        BilinearScalarIdentityFunctor functor;
        BilinearOperator<Matrix_, BilinearScalarIdentityFunctor, Space_>::
          assemble(matrix, functor, space, cubature_name);
      }
    }; // class BilinearScalarIdentityFunctor
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_STANDARD_OPERATORS_HPP
