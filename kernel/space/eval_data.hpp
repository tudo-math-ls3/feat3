#pragma once
#ifndef KERNEL_SPACE_EVAL_DATA_HPP
#define KERNEL_SPACE_EVAL_DATA_HPP 1

// includes, FEAST
#include <kernel/space/base.hpp>
#include <kernel/trafo/eval_data.hpp>

namespace FEAST
{
  namespace Space
  {
    /// \cond internal
    namespace Intern
    {
      template<typename EvalTraits_, bool need_value_>
      struct BasisValueData
      {
        template<typename Data_, typename Evaluator_, typename TrafoEvalData_>
        static void eval(Data_&, const Evaluator_&, const TrafoEvalData_&) {}
      };

      template<typename EvalTraits_, bool need_grad_>
      struct BasisGradientData
      {
        template<typename Data_, typename Evaluator_, typename TrafoEvalData_>
        static void eval(Data_&, const Evaluator_&, const TrafoEvalData_&) {}
      };

      template<typename EvalTraits_, bool need_hess_>
      struct BasisHessianData
      {
        template<typename Data_, typename Evaluator_, typename TrafoEvalData_>
        static void eval(Data_&, const Evaluator_&, const TrafoEvalData_&) {}
      };

      template<typename EvalTraits_>
      struct BasisValueData<EvalTraits_, true>
      {
        /// basis function value object
        typename EvalTraits_::BasisValueType value;

        template<typename Data_, typename Evaluator_, typename TrafoEvalData_>
        static void eval(Data_& data, const Evaluator_& evaluator, const TrafoEvalData_& trafo_data)
        {
          static_assert(evaluator.can_value != 0, "space evaluator does not support basis function values");
          evaluator.eval_values(data, trafo_data);
        }
      };

      template<typename EvalTraits_>
      struct BasisGradientData<EvalTraits_, true>
      {
        /// basis gradient object
        typename EvalTraits_::BasisGradientType grad;

        template<typename Data_, typename Evaluator_, typename TrafoEvalData_>
        static void eval(Data_& data, const Evaluator_& evaluator, const TrafoEvalData_& trafo_data)
        {
          static_assert(evaluator.can_grad != 0, "space evaluator does not support basis function gradients");
          evaluator.eval_gradients(data, trafo_data);
        }
      };

      template<typename EvalTraits_>
      struct BasisHessianData<EvalTraits_, true>
      {
        /// basis hessian object
        typename EvalTraits_::BasisHessianType hess;

        template<typename Data_, typename Evaluator_, typename TrafoEvalData_>
        static void eval(Data_& data, const Evaluator_& evaluator, const TrafoEvalData_& trafo_data)
        {
          static_assert(evaluator.can_hess != 0, "space evaluator does not support basis function hessians");
          evaluator.eval_hessians(data, trafo_data);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Basis function evaluation data structure
     *
     * \tparam EvalTraits_
     * The space evaluator traits that this evaluation data shall use.
     *
     * \tparam Cfg_
     * A space config class that specifies what data shall be supplied. See Space::ConfigBase for details.
     *
     * \author Peter Zajac
     */
    template<
      typename EvalTraits_,
      typename Cfg_>
    class BasisData :
      public Intern::BasisHessianData<EvalTraits_, Cfg_::need_hess != 0>,
      public Intern::BasisGradientData<EvalTraits_, Cfg_::need_grad != 0>,
      public Intern::BasisValueData<EvalTraits_, Cfg_::need_value != 0>
    {
    public:
      /// support enumeration
      enum
      {
        /// specifies whether function values are given
        have_value = Cfg_::need_value,
        /// specifies whether gradients are given
        have_grad = Cfg_::need_grad,
        /// specifies whether hessians are given
        have_hess = Cfg_::need_hess
      };

      /// \cond internal
      typedef Intern::BasisValueData<EvalTraits_, have_value != 0> BasisValueBase;
      typedef Intern::BasisGradientData<EvalTraits_, have_grad != 0> BasisGradientBase;
      typedef Intern::BasisHessianData<EvalTraits_, have_hess != 0> BasisHessianBase;
      /// \endcond

      template<typename Data_, typename Evaluator_, typename TrafoEvalData_>
      static void eval(Data_& data, const Evaluator_& evaluator, const TrafoEvalData_& trafo_data)
      {
        BasisValueBase::eval(data, evaluator, trafo_data);
        BasisGradientBase::eval(data, evaluator, trafo_data);
        BasisHessianBase::eval(data, evaluator, trafo_data);
      }
    }; // class FuncData<...>

    /**
     * \brief Space evaluation data structure
     *
     * \tparam EvalTraits_
     * The space evaluator traits that this evaluation data shall use.
     *
     * \tparam Cfg_
     * A space config class that specifies what data shall be supplied. See Space::ConfigBase for details.
     *
     * \author Peter Zajac
     */
    template<
      typename EvalTraits_,
      typename Cfg_>
    class EvalData
    {
    public:
      /// support enumeration
      enum
      {
        /// specifies whether function values are given
        have_value = Cfg_::need_value,
        /// specifies whether gradients are given
        have_grad = Cfg_::need_grad,
        /// specifies whether hessians are given
        have_hess = Cfg_::need_hess,
        /// maximum number of local dofs
        max_local_dofs = EvalTraits_::max_local_dofs
      };

      /// the basis function data vector
      BasisData<EvalTraits_, Cfg_> phi[max_local_dofs];

      /**
       * \brief Evaluation operator
       *
       * \param[in] evaluator
       * The space evaluator that is to be used for evaluation.
       *
       * \param[in] trafo_data
       * The trafo data structure that specifies the evaluation point.
       */
      template<typename Evaluator_, typename TrafoEvalData_>
      void operator()(const Evaluator_& evaluator, const TrafoEvalData_& trafo_data)
      {
        BasisData<EvalTraits_, Cfg_>::eval(*this, evaluator, trafo_data);
      }
    }; // class EvalData<...>

  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_EVAL_DATA_HPP
