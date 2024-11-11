// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/analytic/function.hpp>
#include <kernel/assembly/linear_functional.hpp>

namespace FEAT
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
          /// Function return value, can be scalar or vector-valued
          typedef typename Analytic::EvalTraits<DataType, Function_>::ValueType ValueType;

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
            _func_value = _func_eval.value(tau.img_point);
          }

          // copy pasted since Doxygen does not like the operator part in
          // \copydoc LinearFunctional::Evaluator::operator()
          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the linear functional for a given test function in a single point.
           *
           * \param[in] psi
           * The \transient test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the linear functional.
           **/
          ValueType eval(const TestBasisData& psi) const
          {
            return _func_value * psi.value;
          }
        }; // class ForceFunctional::Evaluator<...>

      protected:
        /// a reference to the analytic function
        const Function_& _function;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] function
         * A \resident reference to the function that is to be wrapped into a force functional.
         */
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
          /// Function return value, can be scalar or vector-valued
          typedef typename Analytic::EvalTraits<DataType, Function_>::ValueType ValueType;

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
            // compute the function hessian's trace (=laplacian)
            //_func_laplace = - _func_eval.hessian(tau.img_point).trace();
            calc_trace(_func_eval.hessian(tau.img_point));
          }

          template<typename T_, int m_, int n_, int sm_, int sn_>
          void calc_trace(const Tiny::Matrix<T_, m_, n_, sm_, sn_>& hess)
          {
            _func_laplace = - hess.trace();
          }

          template<typename T_, int l_, int m_, int n_, int sl_, int sm_, int sn_>
          void calc_trace(const Tiny::Tensor3<T_, l_, m_, n_, sl_, sm_, sn_>& hess)
          {
            _func_laplace.format();
            for(int i(0); i < hess.l; ++i)
              for(int j(0); j < hess.m; ++j)
                _func_laplace[i] -= hess[i][j][j];
          }

          // copy pasted since Doxygen does not like the operator part in
          // \copydoc LinearFunctional::Evaluator::operator()
          /**
           * \brief Evaluation operator
           *
           * This operator evaluates the linear functional for a given test function in a single point.
           *
           * \param[in] psi
           * The \transient test function data in the current evaluation point. \see Space::EvalData
           *
           * \returns
           * The value of the linear functional.
           **/
          ValueType eval(const TestBasisData& psi) const
          {
            return _func_laplace * psi.value;
          }
        }; // class LaplaceFunctional::Evaluator<...>

      protected:
        /// a reference to the analytic function
        const Function_& _function;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] function
         * A \resident reference to the function that is to be wrapped into a Laplace functional.
         */
        explicit LaplaceFunctional(const Function_& function) :
          _function(function)
        {
        }
      }; // class LaplaceFunctional
    } // namespace Common
  } // namespace Assembly
} // namespace FEAT
