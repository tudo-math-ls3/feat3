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
      };

      template<typename EvalTraits_, bool need_grad_>
      struct BasisGradientData
      {
      };

      template<typename EvalTraits_, bool need_hess_>
      struct BasisHessianData
      {
      };

      template<typename EvalTraits_, bool need_ref_value_>
      struct BasisReferenceValueData
      {
      };

      template<typename EvalTraits_, bool need_ref_grad_>
      struct BasisReferenceGradientData
      {
      };

      template<typename EvalTraits_, bool need_ref_hess_>
      struct BasisReferenceHessianData
      {
      };

      template<typename EvalTraits_>
      struct BasisValueData<EvalTraits_, true>
      {
        /// basis function value object
        typename EvalTraits_::BasisValueType value;
      };

      template<typename EvalTraits_>
      struct BasisGradientData<EvalTraits_, true>
      {
        /// basis gradient object
        typename EvalTraits_::BasisGradientType grad;
      };

      template<typename EvalTraits_>
      struct BasisHessianData<EvalTraits_, true>
      {
        /// basis hessian object
        typename EvalTraits_::BasisHessianType hess;
      };

      template<typename EvalTraits_>
      struct BasisReferenceValueData<EvalTraits_, true>
      {
        /// basis reference value object
        typename EvalTraits_::BasisReferenceValueType ref_value;
      };

      template<typename EvalTraits_>
      struct BasisReferenceGradientData<EvalTraits_, true>
      {
        /// basis reference gradient object
        typename EvalTraits_::BasisReferenceGradientType ref_grad;
      };

      template<typename EvalTraits_>
      struct BasisReferenceHessianData<EvalTraits_, true>
      {
        /// basis reference hessian object
        typename EvalTraits_::BasisReferenceHessianType ref_hess;
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
      public Intern::BasisValueData<EvalTraits_, Cfg_::need_value != 0>,
      public Intern::BasisReferenceHessianData<EvalTraits_, Cfg_::need_ref_hess != 0>,
      public Intern::BasisReferenceGradientData<EvalTraits_, Cfg_::need_ref_grad != 0>,
      public Intern::BasisReferenceValueData<EvalTraits_, Cfg_::need_ref_value != 0>
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
        /// specifies whether reference values are given
        have_ref_value = Cfg_::need_ref_value,
        /// specifies whether reference gradients are given
        have_ref_grad = Cfg_::need_ref_grad,
        /// specifies whether reference hessians are given
        have_ref_hess = Cfg_::need_ref_hess
      };

      /// \cond internal
      typedef Intern::BasisValueData<EvalTraits_, have_value != 0> BasisValueBase;
      typedef Intern::BasisGradientData<EvalTraits_, have_grad != 0> BasisGradientBase;
      typedef Intern::BasisHessianData<EvalTraits_, have_hess != 0> BasisHessianBase;
      typedef Intern::BasisReferenceValueData<EvalTraits_, have_ref_value != 0> BasisReferenceValueBase;
      typedef Intern::BasisReferenceGradientData<EvalTraits_, have_ref_grad != 0> BasisReferenceGradientBase;
      typedef Intern::BasisReferenceHessianData<EvalTraits_, have_ref_hess != 0> BasisReferenceHessianBase;
      /// \endcond
    }; // class BasisData<...>

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
        /// maximum number of local dofs
        max_local_dofs = EvalTraits_::max_local_dofs
      };

      /// basis data type
      typedef BasisData<EvalTraits_, Cfg_> BasisDataType;

      /// the basis function data vector
      BasisDataType phi[max_local_dofs];

      /// support enumeration
      enum
      {
        /// specifies whether function values are given
        have_value = BasisDataType::have_value,
        /// specifies whether gradients are given
        have_grad = BasisDataType::have_grad,
        /// specifies whether hessians are given
        have_hess = BasisDataType::have_hess,
        /// specifies whether reference values are given
        have_ref_value = BasisDataType::have_ref_value,
        /// specifies whether reference gradients are given
        have_ref_grad = BasisDataType::have_ref_grad,
        /// specifies whether reference hessians are given
        have_ref_hess = BasisDataType::have_ref_hess
      };
    }; // class EvalData<...>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_EVAL_DATA_HPP
