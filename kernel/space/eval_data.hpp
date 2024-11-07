// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/space/base.hpp>
#ifndef __CUDA_ARCH__
#include <kernel/util/math.hpp>
#endif

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
     * When compiling in debug mode, all values are initialized to NaN.
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

#if defined(DEBUG) && !defined(__CUDA_ARCH__)
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
     * \brief Reduced Basis function evaluation data structure
     *
     * \tparam EvalTraits_
     * The space evaluator traits that this evaluation data shall use.
     *
     * \tparam Cfg_
     * A space config class that specifies what data shall be supplied. See Space::ConfigBase for details.
     *
     * \author Maximilian Esser
     */
    template<
      typename EvalTraits_,
      SpaceTags cfg_tags_>
    class BasisDataReduced
    {
    public:
      // union
      // {
        /// basis gradient object
        typename EvalTraits_::BasisGradientType grad;
        /// basis reference gradient object
        typename EvalTraits_::BasisReferenceGradientType ref_grad;
      // };
      union
      {
        /// basis function value object
        typename EvalTraits_::BasisValueType value;
        /// basis reference value object
        typename EvalTraits_::BasisReferenceValueType ref_value = value;
      };

      // our space config
      static constexpr SpaceTags config = cfg_tags_;

      BasisDataReduced() : grad(), value() {};

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
     * When compiling in debug mode, all values are initialized to NaN.
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

    /**
     * \brief Reduced Space evaluation data structure
     *
     * \tparam EvalTraits_
     * The space evaluator traits that this evaluation data shall use.
     *
     * \tparam Cfg_
     * A space config class that specifies what data shall be supplied. See Space::ConfigBase for details.
     *
     * \note
     * When compiling in debug mode, all values are initialized to NaN.
     *
     * \author Maximilian Esser
     */
    template<
      typename EvalTraits_,
      SpaceTags cfg_tags_>
    class EvalDataReduced
    {
    public:
      /// maximum number of local dofs
      static constexpr int max_local_dofs = EvalTraits_::max_local_dofs;

      /// basis data type
      typedef BasisDataReduced<EvalTraits_, cfg_tags_> BasisDataType;

      typedef EvalTraits_ EvalTraits;

      /// the basis function data vector
      BasisDataType phi[max_local_dofs];

      // our space config
      static constexpr SpaceTags config = BasisDataType::config;
    }; // class EvalData<...>
  } // namespace Space
} // namespace FEAT
