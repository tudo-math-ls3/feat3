#pragma once
#ifndef KERNEL_SPACE_EVAL_DATA_HPP
#define KERNEL_SPACE_EVAL_DATA_HPP 1

// includes, FEAST
#include <kernel/space/base.hpp>

namespace FEAST
{
  namespace Space
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Evaluator_,
        bool need_value_>
      struct EvalValueData
      {
        template<typename TrafoEvalData_>
        void eval(const Evaluator_&, const TrafoEvalData_&) {}
      };

      template<
        typename Evaluator_,
        bool need_grad_>
      struct EvalGradientData
      {
        template<typename TrafoEvalData_>
        void eval(const Evaluator_&, const TrafoEvalData_&) {}
      };

      template<typename Evaluator_>
      struct EvalValueData<Evaluator_, true>
      {
        static_assert(Evaluator_::can_value != 0, "space evaluator can't compute basis function values");

        /// value vector
        typename Evaluator_::BasisValueVectorType values;

        template<typename TrafoEvalData_>
        void eval(const Evaluator_& evaluator, const TrafoEvalData_& trafo_data)
        {
          evaluator.eval_values(values, trafo_data);
        }
      };

      template<typename Evaluator_>
      struct EvalGradientData<Evaluator_, true>
      {
        static_assert(Evaluator_::can_grad != 0, "space evaluator can't compute basis function gradients");

        /// gradient vector
        typename Evaluator_::BasisGradientVectorType grads;

        template<typename TrafoEvalData_>
        void eval(const Evaluator_& evaluator, const TrafoEvalData_& trafo_data)
        {
          evaluator.eval_gradients(grads, trafo_data);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Space evaluation data structure
     *
     * \tparam Evaluator_
     * The space evaluator that this evaluation data shall use.
     *
     * \tparam Cfg_
     * A space config class that specifies what data shall be supplied. See Space::ConfigBase for details.
     *
     * \author Peter Zajac
     */
    template<
      typename Evaluator_,
      typename Cfg_>
    class EvalData :
      public Intern::EvalValueData<Evaluator_, Cfg_::need_value != 0>,
      public Intern::EvalGradientData<Evaluator_, Cfg_::need_grad != 0>
    {
    public:
      /// support enumeration
      enum
      {
        /// specifies whether function values are given
        have_value = Cfg_::need_value,
        /// specifies whether gradients are given
        have_grad = Cfg_::need_grad
      };

      /// \cond internal
      typedef Intern::EvalValueData<Evaluator_, have_value != 0> EvalValueBase;
      typedef Intern::EvalGradientData<Evaluator_, have_grad != 0> EvalGradientBase;
      /// \endcond

      /**
       * \brief Evaluation operator
       *
       * \param[in] evaluator_
       * The space evaluator that is to be used for evaluation.
       *
       * \param[in] trafo_data
       * The trafo data structure that specifies the evaluation point.
       */
      template<typename TrafoEvalData_>
      void operator()(const Evaluator_& evaluator, const TrafoEvalData_& trafo_data)
      {
        EvalValueBase::eval(evaluator, trafo_data);
        EvalGradientBase::eval(evaluator, trafo_data);
      }
    }; // class EvalData<...>

    /* ***************************************************************************************** */

    /// \cond internal
    namespace Intern
    {
      template<typename Evaluator_, bool need_value_>
      struct FuncValueData
      {
        explicit FuncValueData(const EvalValueData<Evaluator_, false>&, Index) {}
      };

      template<typename Evaluator_, bool need_grad_>
      struct FuncGradientData
      {
        explicit FuncGradientData(const EvalGradientData<Evaluator_, false>&, Index) {}
      };

      template<typename Evaluator_>
      struct FuncValueData<Evaluator_, true>
      {
        /// value reference
        typename Evaluator_::BasisValueConstRef value;

        explicit FuncValueData(const EvalValueData<Evaluator_, true>& value_data, Index i) :
          value(value_data.values[i])
        {
        }
      };

      template<typename Evaluator_>
      struct FuncGradientData<Evaluator_, true>
      {
        /// gradient reference
        typename Evaluator_::BasisGradientConstRef grad;

        explicit FuncGradientData(const EvalGradientData<Evaluator_, true>& grad_data, Index i) :
          grad(grad_data.grads[i])
        {
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Basis function evaluation data structure
     *
     * \tparam Evaluator_
     * The space evaluator that this evaluation data shall use.
     *
     * \tparam Cfg_
     * A space config class that specifies what data shall be supplied. See Space::ConfigBase for details.
     *
     * \author Peter Zajac
     */
    template<
      typename Evaluator_,
      typename Cfg_>
    class FuncData :
      public Intern::FuncValueData<Evaluator_, Cfg_::need_value != 0>,
      public Intern::FuncGradientData<Evaluator_, Cfg_::need_grad != 0>
    {
    public:
      /// support enumeration
      enum
      {
        /// specifies whether function values are given
        have_value = Cfg_::need_value,
        /// specifies whether gradients are given
        have_grad = Cfg_::need_grad
      };

      /// \cond internal
      typedef Intern::FuncValueData<Evaluator_, have_value != 0> FuncValueBase;
      typedef Intern::FuncGradientData<Evaluator_, have_grad != 0> FuncGradientBase;
      /// \endcond

    public:
      /**
       * \brief Constructor
       *
       * \param[in] eval_data
       * The space evaluation data structure that contains the evaluated data.
       *
       * \param[in] i
       * The index of the basis function which values are to be referenced.
       */
      explicit FuncData(const EvalData<Evaluator_, Cfg_>& eval_data, Index i) :
        FuncValueBase(eval_data, i),
        FuncGradientBase(eval_data, i)
      {
      }
    }; // class FuncData<...>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_EVAL_DATA_HPP
