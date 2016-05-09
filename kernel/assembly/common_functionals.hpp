#pragma once
#ifndef KERNEL_ASSEMBLY_COMMON_FUNCTIONALS_HPP
#define KERNEL_ASSEMBLY_COMMON_FUNCTIONALS_HPP 1

// includes, FEAST
#include <kernel/analytic/function.hpp>
#include <kernel/assembly/linear_functional.hpp>

namespace FEAST
{
  namespace Assembly
  {
    namespace Common
    {
      /**
       * \brief Force functional implementation
       *
       * This functor implements the LinearFunctional interface with
       *   \f[ \ell(\varphi) := \int_\Omega f\cdot\varphi \f]
       * for an analytic function \e f.
       *
       * \tparam Function_
       * A template class implementing the AnalyticFunction interface representing the function \e f.
       *
       * \author Peter Zajac
       */
      template<typename Function_>
      class ForceFunctional :
        public LinearFunctional
      {
      public:
        /// ensure that the function can compute values
        static_assert(Function_::can_value, "function can't compute values");

        /// trafo config tag
        static constexpr TrafoTags trafo_config = TrafoTags::img_point;

        /// test space config tag
        static constexpr SpaceTags test_config = SpaceTags::value;

        /**
         * \brief Force Functional Evaluator class template
         *
         * \tparam AsmTraits_
         * The assembly traits class.
         *
         * \author Peter Zajac
         */
        template<typename AsmTraits_>
        class Evaluator :
          public LinearFunctional::Evaluator<AsmTraits_>
        {
        public:
          /// trafo evaluator
          typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;
          /// data type
          typedef typename AsmTraits_::DataType DataType;
          /// trafo data type
          typedef typename AsmTraits_::TrafoData TrafoData;
          /// test-function data type
          typedef typename AsmTraits_::TestBasisData TestBasisData;

        protected:
          /// declare our analytic eval traits
          typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;
          /// the function evaluator
          typename Function_::template Evaluator<AnalyticEvalTraits> _func_eval;
          /// the function value in the current point
          typename AnalyticEvalTraits::ValueType _func_value;

        public:
          /// constructor
          explicit Evaluator(const ForceFunctional<Function_>& functional) :
            _func_eval(functional._function),
            _func_value(DataType(0))
          {
          }

          /** \copydoc LinearFunctional::Evaluator::set_point() */
          void set_point(const TrafoData& tau)
          {
            // evaluate function value
            _func_eval.value(_func_value, tau.img_point);
          }

          // copy pasted since Doxygen does not like the operator part in
          // \copydoc LinearFunctional::Evaluator::operator()
          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the linear functional for a given test function in a single point.
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the linear functional.
           **/
          DataType operator()(const TestBasisData& psi) const
          {
            return _func_value * psi.value;
          }
        }; // class ForceFunctional::Evaluator<...>

      protected:
        /// a reference to the analytic function
        const Function_& _function;

      public:
        /// constructor
        explicit ForceFunctional(const Function_& function) :
          _function(function)
        {
        }
      }; // class ForceFunctional

      /**
       * \brief Laplace functional implementation
       *
       * This functor implements the LinearFunctional interface with
       *   \f[ \ell(\varphi) := \int_\Omega -\Delta f\cdot\varphi \f]
       * for an analytic function \e f.
       *
       * \tparam Function_
       * A template class implementing the AnalyticFunction interface representing the function \e f.
       *
       * \author Peter Zajac
       */
      template<typename Function_>
      class LaplaceFunctional :
        public LinearFunctional
      {
      public:
        /// ensure that the function can compute hessian matrices
        static_assert(Function_::can_hess, "function can't compute hessians");

        /// trafo config tag
        static constexpr TrafoTags trafo_config = TrafoTags::img_point;

        /// test space config tag
        static constexpr SpaceTags test_config = SpaceTags::value;

        /**
         * \brief Force Functional Evaluator class template
         *
         * \tparam AsmTraits_
         * The assembly traits class.
         *
         * \author Peter Zajac
         */
        template<typename AsmTraits_>
        class Evaluator :
          public LinearFunctional::Evaluator<AsmTraits_>
        {
        public:
          /// trafo evaluator
          typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;
          /// data type
          typedef typename AsmTraits_::DataType DataType;
          /// trafo data type
          typedef typename AsmTraits_::TrafoData TrafoData;
          /// test-function data type
          typedef typename AsmTraits_::TestBasisData TestBasisData;

        protected:
          /// declare our analytic eval traits
          typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;
          /// the function evaluator
          typename Function_::template Evaluator<AnalyticEvalTraits> _func_eval;
          /// the function laplacian in the current point
          typename AnalyticEvalTraits::ValueType _func_laplace;

        public:
          /// constructor
          explicit Evaluator(const LaplaceFunctional<Function_>& functional) :
            _func_eval(functional._function),
            _func_laplace(DataType(0))
          {
          }

          /** \copydoc LinearFunctional::Evaluator::set_point() */
          void set_point(const TrafoData& tau)
          {
            // evaluate function hessian
            typename AnalyticEvalTraits::HessianType hess;
            _func_eval.hessian(hess, tau.img_point);

            // compute the hessian's trace (=laplacian)
            _func_laplace = - hess.trace();
          }

          // copy pasted since Doxygen does not like the operator part in
          // \copydoc LinearFunctional::Evaluator::operator()
          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the linear functional for a given test function in a single point.
           *
           * \param[in] psi
           * The test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the linear functional.
           **/
          DataType operator()(const TestBasisData& psi) const
          {
            return _func_laplace * psi.value;
          }
        }; // class LaplaceFunctional::Evaluator<...>

      protected:
        /// a reference to the analytic function
        const Function_& _function;

      public:
        /// constructor
        explicit LaplaceFunctional(const Function_& function) :
          _function(function)
        {
        }
      }; // class LaplaceFunctional
    } // namespace Common
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_COMMON_FUNCTIONALS_HPP
