#pragma once
#ifndef KERNEL_ASSEMBLY_COMMON_OPERATORS_HPP
#define KERNEL_ASSEMBLY_COMMON_OPERATORS_HPP 1

// includes, FEAT
#include <kernel/assembly/bilinear_operator.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Assembly Common namespace
     *
     * This namespace encapsulated commonly used operators and functionals for use with the
     * various assembly classes, which are often used in standard benchmark problems.
     */
    namespace Common
    {
      /**
       * \brief -Laplace operator implementation
       *
       * This functor implements the weak formulation of the bilinear scalar Laplace operator, i.e.
       *   \f[ \nabla \varphi \cdot \nabla\psi \f]
       *
       * This functor can be used with the BilinearOperator assembly class template to assemble a
       * scalar Laplace/Stiffness matrix.
       *
       * \author Peter Zajac
       */
      class LaplaceOperator :
        public BilinearOperator
      {
      public:
        static constexpr TrafoTags trafo_config = TrafoTags::none;
        static constexpr SpaceTags test_config = SpaceTags::grad;
        static constexpr SpaceTags trial_config = SpaceTags::grad;

        /**
         * \brief Laplace evaluator class template
         *
         * \tparam AsmTraits_
         * The assembly traits class.
         *
         * \author Peter Zajac
         */
        template<typename AsmTraits_>
        class Evaluator :
          public BilinearOperator::Evaluator<AsmTraits_>
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
           * \param[in] operat
           * A reference to the Laplace operator object.
           */
          explicit Evaluator(const LaplaceOperator& DOXY(operat))
          {
          }

          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
           * a single point.
           *
           * \param[in] phi
           * The trial function data in the current evaluation point. \see Space::EvalData
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the bilinear functor.
           **/
          DataType operator()(const TrialBasisData& phi, const TestBasisData& psi)
          {
            return dot(phi.grad, psi.grad);
          }
        }; // class LaplaceOperator::Evaluator<...>
      }; // class LaplaceOperator

      /**
       * \brief Identity operator implementation
       *
       * This functor implements the weak formulation of the bilinear scalar Identity operator, i.e.
       *   \f[ \varphi \cdot \psi \f]
       *
       * This functor can be used with the BilinearOperator assembler class template to assemble a
       * scalar mass matrix.
       *
       * \author Peter Zajac
       */
      class IdentityOperator :
        public BilinearOperator
      {
      public:
        static constexpr TrafoTags trafo_config = TrafoTags::none;
        static constexpr SpaceTags test_config = SpaceTags::value;
        static constexpr SpaceTags trial_config = SpaceTags::value;

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
          public BilinearOperator::Evaluator<AsmTraits_>
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
           * \param[in] operat
           * A reference to the Identity operator object.
           */
          explicit Evaluator(const IdentityOperator& DOXY(operat))
          {
          }

          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
           * a single point.
           *
           * \param[in] phi
           * The trial function data in the current evaluation point. \see Space::EvalData
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the bilinear functor.
           **/
          DataType operator()(const TrialBasisData& phi, const TestBasisData& psi)
          {
            return phi.value * psi.value;
          }
        }; // class IdentityOperator::Evaluator<...>
      }; // class IdentityOperator

      /**
       * \brief Vector-valued identity operator implementation
       *
       * This functor can be used with the BilinearOperator assembly class template to assemble one
       * scalar matrix corresponding to one block of the \f$ d \times d \f$ block matrix, where \f$ d \f$ is the
       * number of space dimensions.
       *
       * \author Jordi Paul
       */
      template<int dimension_>
      class IdentityOperatorBlocked :
        public BilinearOperator
      {
      public:
        /// Every block is a dimension x dimension matrix
        static constexpr int BlockHeight = dimension_;
        /// Every block is a dimension x dimension matrix
        static constexpr int BlockWidth = dimension_;

        static constexpr TrafoTags trafo_config = TrafoTags::none;
        static constexpr SpaceTags test_config = SpaceTags::value;
        static constexpr SpaceTags trial_config = SpaceTags::value;

        /**
         * \brief Vector-valued identity operator evaluator class template
         *
         * \tparam AsmTraits_
         * The assembly traits class.
         *
         * \author Jordi Paul
         */
        template<typename AsmTraits_>
        class Evaluator :
            public BilinearOperator::Evaluator<AsmTraits_>
        {
        public:
          /// the data type to be used
          typedef typename AsmTraits_::DataType DataType;
          /// the data type for the block system
          typedef typename AsmTraits_::OperatorValueType OperatorValueType;
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
           * \param[in] operat
           * A reference to the Du : Dv operator object.
           */
          explicit Evaluator(const IdentityOperatorBlocked& DOXY(operat))
          {
          }

          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
           * a single point.
           *
           * \param[in] phi
           * The trial function data in the current evaluation point. \see Space::EvalData
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the bilinear functor.
           **/
          OperatorValueType operator()(const TrialBasisData& phi, const TestBasisData& psi)
          {
            OperatorValueType r(DataType(0));
            for(int i(0); i < OperatorValueType::m; ++i)
            {
              r(i,i) = phi.value * psi.value;
            }
            return r;
          }
        }; // class IdentityOperatorBlocked::Evaluator<...>
      }; // class DuDVOperatorBlocked

      /**
       * \brief Trial-Derivative operator implementation
       *
       * This functor implements the weak formulation of the bilinear trial-function derivative operator, i.e.
       *   \f[ \partial_i \varphi \cdot \psi \f]
       *
       * This functor can be used with the BilinearOperator assembly class template to assemble a
       * scalar matrix for the pressure-gradient operator of the Stokes equation.
       *
       * \author Peter Zajac
       */
      class TrialDerivativeOperator :
        public BilinearOperator
      {
      public:
        /// the desired derivative
        int deriv;

        static constexpr TrafoTags trafo_config = TrafoTags::none;
        static constexpr SpaceTags test_config = SpaceTags::grad;
        static constexpr SpaceTags trial_config = SpaceTags::value;

        /**
         * \brief Trial-Derivative evaluator class template
         *
         * \tparam AsmTraits_
         * The assembly traits class.
         *
         * \author Peter Zajac
         */
        template<typename AsmTraits_>
        class Evaluator :
          public BilinearOperator::Evaluator<AsmTraits_>
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

        protected:
          /// the desired derivative
          int deriv;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] operat
           * A reference to the Laplace operator object.
           */
          explicit Evaluator(const TrialDerivativeOperator& operat) :
            deriv(operat.deriv)
          {
          }

          // copy pasted since Doxygen does not like the operator part in
          // \copydoc BilinearOperator::Evaluator::operator()
          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
           * a single point.
           *
           * \param[in] phi
           * The trial function data in the current evaluation point. \see Space::EvalData
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the bilinear functor.
           **/
          DataType operator()(const TrialBasisData& phi, const TestBasisData& psi)
          {
            return phi.value * psi.grad[deriv];
          }
        }; // class TrialDerivativeOperator::Evaluator<...>

        /*
         *
         * \param[in] derivative
         * The index of the derivative for this operator:
         *  - 0: X-derivative
         *  - 1: Y-derivative
         *  - 2: Z-derivative
         *  - ...
        */
        explicit TrialDerivativeOperator(int derivative) :
          deriv(derivative)
        {
        }
      }; // class TrialDerivativeOperator

      /**
       * \brief Test-Derivative operator implementation
       *
       * This functor implements the weak formulation of the bilinear test-function derivative operator, i.e.
       *   \f[ \varphi \cdot \partial_i \psi \f]
       *
       * This functor can be used with the BilinearOperator assembly class template to assemble a
       * scalar matrix for the pressure-gradient operator of the Stokes equation.
       *
       * \author Peter Zajac
       */
      class TestDerivativeOperator :
        public BilinearOperator
      {
      public:
        /// the desired derivative
        int deriv;

        static constexpr TrafoTags trafo_config = TrafoTags::none;
        static constexpr SpaceTags test_config = SpaceTags::grad;
        static constexpr SpaceTags trial_config = SpaceTags::value;

        /**
         * \brief Test-Derivative evaluator class template
         *
         * \tparam AsmTraits_
         * The assembly traits class.
         *
         * \author Peter Zajac
         */
        template<typename AsmTraits_>
        class Evaluator :
          public BilinearOperator::Evaluator<AsmTraits_>
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

        protected:
          /// the desired derivative
          int deriv;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] operat
           * A reference to the Laplace operator object.
           */
          explicit Evaluator(const TestDerivativeOperator& operat) :
            deriv(operat.deriv)
          {
          }

          // copy pasted since Doxygen does not like the operator part in
          // \copydoc BilinearOperator::Evaluator::operator()
          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
           * a single point.
           *
           * \param[in] phi
           * The trial function data in the current evaluation point. \see Space::EvalData
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the bilinear functor.
           **/
          DataType operator()(const TrialBasisData& phi, const TestBasisData& psi)
          {
            return phi.value * psi.grad[deriv];
          }
        }; // class TestDerivativeOperator::Evaluator<...>

        /*
         *
         * \param[in] derivative
         * The index of the derivative for this operator:
         *  - 0: X-derivative
         *  - 1: Y-derivative
         *  - 2: Z-derivative
         *  - ...
        */
        explicit TestDerivativeOperator(int derivative) :
          deriv(derivative)
        {
        }
      }; // class TestDerivativeOperator

      /**
       * \brief div(phi) * div(psi) operator implementation
       *
       * \author Peter Zajac
       */
      class DivDivOperator :
        public Assembly::BilinearOperator
      {
      public:
        /// Row index
        int ir;
        /// Column index
        int ic;

        static constexpr TrafoTags trafo_config = TrafoTags::none;
        static constexpr SpaceTags test_config = SpaceTags::grad;
        static constexpr SpaceTags trial_config = SpaceTags::grad;

        template<typename AsmTraits_>
        class Evaluator :
          public Assembly::BilinearOperator::Evaluator<AsmTraits_>
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

        protected:
          /// Row index
          int ir;
          /// Column index
          int ic;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] operat
           * A reference to the operator object.
           */
          explicit Evaluator(const DivDivOperator& operat) :
            ir(operat.ir), ic(operat.ic)
          {
          }

          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
           * a single point.
           *
           * \param[in] phi
           * The trial function data in the current evaluation point. \see Space::EvalData
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the bilinear functor.
           **/
          DataType operator()(const TrialBasisData& phi, const TestBasisData& psi)
          {
            return phi.grad[ic] * psi.grad[ir];
          }
        }; // class DivDivOperator::Evaluator<...>

        /**
         * \brief Constructor
         *
         * \param[in] ir_
         * The row index of the block
         *
         * \param[in] ic_
         * The column index of the block
         */
        explicit DivDivOperator(int ir_, int ic_) :
          ir(ir_), ic(ic_)
        {
        }
      }; // class DivDivOperator

      /**
       * \brief Du:Dv operator implementation
       *
       * This functor implements the weak formulation of the bilinear Du : Dv operator, i.e.
       * \f[
       *   \mathbf{D} \varphi : \mathbf{D} \psi
       * \f]
       * with the operator
       * \f[
       *   \mathbf{D} \varphi = \frac{1}{2} \left( \nabla \varphi + \left( \nabla \varphi \right)^T \right).
       * \f]
       *
       * Note that
       * \f[
       *   \mathbf{D} \varphi : \mathbf{D} \psi = 2 \left( \nabla \varphi : \nabla \psi + \nabla \varphi : \left( \nabla \psi \right)^T \right),
       * \f]
       * so the \f$ (k,l) \f$-block consists of entries corresponding to \f$ \partial_k \varphi \partial_l \psi \f$.
       *
       *
       * This functor can be used with the BilinearOperator assembly class template to assemble one
       * scalar matrix corresponding to one block of the \f$ d \times d \f$ block matrix, where \f$ d \f$ is the
       * number of space dimensions.
       *
       * \author Jordi Paul
       */
      class DuDvOperator :
        public BilinearOperator
      {
      public:
        /// Row index
        int ir;
        /// Column index
        int ic;

        static constexpr TrafoTags trafo_config = TrafoTags::none;
        static constexpr SpaceTags test_config = SpaceTags::grad;
        static constexpr SpaceTags trial_config = SpaceTags::grad;

        /**
         * \brief Du : Dv evaluator class template
         *
         * \tparam AsmTraits_
         * The assembly traits class.
         *
         * \author Jordi Paul
         */
        template<typename AsmTraits_>
        class Evaluator :
            public BilinearOperator::Evaluator<AsmTraits_>
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

        protected:
          /// Row index
          int ir;
          /// Column index
          int ic;

        public:
          /**
           * \brief Constructor
           *
           * \param[in] operat
           * A reference to the Du : Dv operator object.
           */
          explicit Evaluator(const DuDvOperator& operat) :
            ir(operat.ir), ic(operat.ic)
          {
          }

          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
           * a single point.
           *
           * \param[in] phi
           * The trial function data in the current evaluation point. \see Space::EvalData
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the bilinear functor.
           **/
          DataType operator()(const TrialBasisData& phi, const TestBasisData& psi)
          {
            return ((ir==ic) ? dot(phi.grad,psi.grad) : DataType(0)) + phi.grad[ir] * psi.grad[ic];
          }
        }; // class DuDvOperator::Evaluator<...>

        /**
         * \brief Constructor
         *
         * \param[in] ir_
         * The row index of the block
         *
         * \param[in] ic_
         * The column index of the block
         */
        explicit DuDvOperator(int ir_, int ic_) :
          ir(ir_), ic(ic_)
        {
        }
      }; // class DuDvOperator

      /**
       * \brief Du:Dv operator implementation
       *
       * This functor implements the weak formulation of the bilinear Du : Dv operator, i.e.
       * \f[
       *   \mathbf{D} \varphi : \mathbf{D} \psi
       * \f]
       * with the operator
       * \f[
       *   \mathbf{D} \varphi = \frac{1}{2} \left( \nabla \varphi + \left( \nabla \varphi \right)^T \right).
       * \f]
       *
       * This functor can be used with the BilinearOperator assembly class template to assemble one
       * scalar matrix corresponding to one block of the \f$ d \times d \f$ block matrix, where \f$ d \f$ is the
       * number of space dimensions.
       *
       * Note that
       * \f[
       *   \mathbf{D} \varphi : \mathbf{D} \psi = 2 \left( \nabla \varphi : \nabla \psi + \nabla \varphi : \left( \nabla \psi \right)^T \right),
       * \f]
       * so the \f$ (k,l) \f$-block consists of entries corresponding to \f$ \partial_k \varphi \partial_l \psi \f$.
       *
       * \author Jordi Paul
       */
      template<int dimension_>
      class DuDvOperatorBlocked :
        public BilinearOperator
      {
      public:
        /// Every block is a dimension x dimension matrix
        static constexpr int BlockHeight = dimension_;
        /// Every block is a dimension x dimension matrix
        static constexpr int BlockWidth = dimension_;

        static constexpr TrafoTags trafo_config = TrafoTags::none;
        static constexpr SpaceTags test_config = SpaceTags::grad;
        static constexpr SpaceTags trial_config = SpaceTags::grad;

        /**
         * \brief Du : Dv evaluator class template
         *
         * \tparam AsmTraits_
         * The assembly traits class.
         *
         * \author Jordi Paul
         */
        template<typename AsmTraits_>
        class Evaluator :
            public BilinearOperator::Evaluator<AsmTraits_>
        {
        public:
          /// the data type to be used
          typedef typename AsmTraits_::DataType DataType;
          /// the data type for the block system
          typedef typename AsmTraits_::OperatorValueType OperatorValueType;
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
           * \param[in] operat
           * A reference to the Du : Dv operator object.
           */
          explicit Evaluator(const DuDvOperatorBlocked& DOXY(operat))
          {
          }

          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the bilinear operator for a given combination of test- and trial-functions in
           * a single point.
           *
           * \param[in] phi
           * The trial function data in the current evaluation point. \see Space::EvalData
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the bilinear functor.
           **/
          OperatorValueType operator()(const TrialBasisData& phi, const TestBasisData& psi)
          {
            OperatorValueType r(DataType(0));
            for(int i(0); i < OperatorValueType::m; ++i)
            {
              r(i,i) = dot(phi.grad, psi.grad) + phi.grad[i]*psi.grad[i];
              for(int j(0); j < i; ++j)
              {
                r(i,j) = phi.grad[i] * psi.grad[j];
                r(j,i) = phi.grad[j] * psi.grad[i];
              }
            }
            return r;
          }
        }; // class DuDvOperatorBlocked::Evaluator<...>
      }; // class DuDVOperatorBlocked
    } // namespace Common
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_COMMON_OPERATORS_HPP
