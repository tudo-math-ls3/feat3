#pragma once
#ifndef KERNEL_ASSEMBLY_COMMON_FUNCTIONALS_HPP
#define KERNEL_ASSEMBLY_COMMON_FUNCTIONALS_HPP 1

// includes, FEAST
#include <kernel/assembly/linear_functional.hpp>
#include <kernel/assembly/analytic_function.hpp>

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
        static_assert(Function_::can_value != 0, "function can't compute values");

        /// function config tag
        struct FunctionConfig :
          public Trafo::AnalyticConfigBase
        {
          /// dummy enum
          enum
          {
            /// we need function values
            need_value = 1
          };
        };

        /// trafo config tag
        typedef typename Function_::template ConfigTraits<FunctionConfig>::TrafoConfig TrafoConfig;

        /// test space config tag
        struct TestConfig :
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
          /// the function evaluator
          typename Function_::template Evaluator<typename AsmTraits_::AnalyticEvalTraits> _func_eval;

        public:
          /// constructor
          explicit Evaluator(const ForceFunctional<Function_>& functional) :
            _func_eval(functional._function)
          {
          }

          /** \copydoc LinearFunctional::Evaluator::prepare() */
          void prepare(const TrafoEvaluator& trafo_eval)
          {
            // prepare function evaluator
            _func_eval.prepare(trafo_eval);
          }

          void finish()
          {
            // finish function evaluator
            _func_eval.finish();
          }

          /** \copydoc LinearFunctional::Evaluator::operator() */
          DataType operator()(const TrafoData& tau, const TestBasisData& psi) const
          {
            return _func_eval.value(tau) * psi.value;
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
    } // namespace Common
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_COMMON_FUNCTIONALS_HPP
