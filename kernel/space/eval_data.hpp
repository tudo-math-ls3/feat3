#pragma once
#ifndef KERNEL_SPACE_EVAL_DATA_HPP
#define KERNEL_SPACE_EVAL_DATA_HPP 1

// includes, FEAST
#include <kernel/space/base.hpp>

namespace FEAST
{
  namespace Space
  {
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
      SpaceTags cfg_tags_>
    class BasisData
    {
    public:
      /// basis hessian object
      typename EvalTraits_::BasisHessianType hess;
      /// basis gradient object
      typename EvalTraits_::BasisGradientType grad;
      /// basis function value object
      typename EvalTraits_::BasisValueType value;

      /// basis reference hessian object
      typename EvalTraits_::BasisReferenceHessianType ref_hess;
      /// basis reference gradient object
      typename EvalTraits_::BasisReferenceGradientType ref_grad;
      /// basis reference value object
      typename EvalTraits_::BasisReferenceValueType ref_value;

      // our space config
      static constexpr SpaceTags config = cfg_tags_;
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
      SpaceTags cfg_tags_>
    class EvalData
    {
    public:
      /// maximum number of local dofs
      static constexpr int max_local_dofs = EvalTraits_::max_local_dofs;

      /// basis data type
      typedef BasisData<EvalTraits_, cfg_tags_> BasisDataType;

      /// the basis function data vector
      BasisDataType phi[max_local_dofs];

      // our space config
      static constexpr SpaceTags config = BasisDataType::config;
    }; // class EvalData<...>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_EVAL_DATA_HPP
