// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_EVAL_DATA_HPP
#define KERNEL_SPACE_EVAL_DATA_HPP 1

// includes, FEAT
#include <kernel/space/base.hpp>

namespace FEAT
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
     * \note
     * When compiling in debug mode, all values are initialised to NaN.
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

#ifdef DEBUG
      BasisData()
      {
        // format all values to NaN
        const auto qnan = Math::nan<typename EvalTraits_::DataType>();
        hess = qnan;
        grad = qnan;
        value = qnan;
        ref_hess = qnan;
        ref_grad = qnan;
        ref_value = qnan;
      }
#endif // DEBUG
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
     * \note
     * When compiling in debug mode, all values are initialised to NaN.
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
} // namespace FEAT

#endif // KERNEL_SPACE_EVAL_DATA_HPP
