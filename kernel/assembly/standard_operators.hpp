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
        typedef typename AsmTraits_::TestBasisData TestBasisData;
        /// the assembler's trial-function data type
        typedef typename AsmTraits_::TrialBasisData TrialBasisData;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] functor
         * A reference to the scalar Laplace functor.
         */
        explicit Evaluator(const BilinearScalarLaplaceFunctor& DOXY(functor))
        {
        }

        /** \copydoc BilinearFunctorBase::Evaluator::operator() */
        DataType operator()(const TrafoData& DOXY(tau), const TrialBasisData& phi, const TestBasisData& psi)
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
       * \param[in] cubature_name
       * A string containing the name of the cubature rule to be used for integration.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       */
      template<
        typename Matrix_,
        typename Space_>
      static void assemble_matrix(
        Matrix_& matrix,
        const String& cubature_name,
        const Space_& space,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
      {
        BilinearScalarLaplaceFunctor functor;
        Cubature::DynamicFactory cubature_factory(cubature_name);
        BilinearOperator<BilinearScalarLaplaceFunctor>::
          assemble_matrix1(matrix, functor, cubature_factory, space, alpha);
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
        typedef typename AsmTraits_::TestBasisData TestBasisData;
        /// the assembler's trial-function data type
        typedef typename AsmTraits_::TrialBasisData TrialBasisData;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] functor
         * A reference to the scalar Identity functor.
         */
        explicit Evaluator(const BilinearScalarIdentityFunctor& DOXY(functor))
        {
        }

        /** \copydoc BilinearFunctorBase::Evaluator::operator() */
        DataType operator()(const TrafoData& DOXY(tau), const TrialBasisData& phi, const TestBasisData& psi)
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
      static void assemble_matrix1(
        Matrix_& matrix,
        const String& cubature_name,
        const Space_& space,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
      {
        BilinearScalarIdentityFunctor functor;
        Cubature::DynamicFactory cubature_factory(cubature_name);
        BilinearOperator<BilinearScalarIdentityFunctor>::
          assemble_matrix1(matrix, functor, cubature_factory, space, alpha);
      }

      /**
       * \brief Assembles a mass matrix.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] cubature_name
       * A string containing the name of the cubature rule to be used for integration.
       *
       * \param[in] test_space
       * A reference to the finite-element to be used as the test-space.
       *
       * \param[in] trial_space
       * A reference to the finite-element to be used as the trial-space.
       */
      template<
        typename Matrix_,
        typename TestSpace_,
        typename TrialSpace_>
      static void assemble_matrix2(
        Matrix_& matrix,
        const String& cubature_name,
        const TestSpace_& test_space,
        const TrialSpace_& trial_space,
        typename Matrix_::DataType alpha = typename Matrix_::DataType(1))
      {
        BilinearScalarIdentityFunctor functor;
        Cubature::DynamicFactory cubature_factory(cubature_name);
        BilinearOperator<BilinearScalarIdentityFunctor>::
          assemble_matrix2(matrix, functor, cubature_factory, test_space, trial_space, alpha);
      }
    }; // class BilinearScalarIdentityFunctor
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_STANDARD_OPERATORS_HPP
