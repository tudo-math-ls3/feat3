// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_EVAL_DATA_HPP
#define KERNEL_TRAFO_EVAL_DATA_HPP 1

// includes, FEAT
#include <kernel/trafo/base.hpp>

namespace FEAT
{
  namespace Trafo
  {
    /**
     * \brief Trafo evaluation data structure.
     *
     * \tparam EvalTraits_
     * The trafo evaluator traits that this evaluation data shall use.
     *
     * \tparam cfg_tags_
     * A trafo data config class that specifies what data shall be supplied.
     *
     * \author Peter Zajac
     */
    template<
      typename EvalTraits_,
      TrafoTags cfg_tags_>
    class EvalData
    {
    public:
      /// trafo evaluation traits
      typedef EvalTraits_ EvalTraits;

      // Note:
      // The following members are ordered by non-ascending size
      // to avoid unnecessary data alignment padding.

      /// inverse hessian tensor
      typename EvalTraits::HessianInverseType hess_inv;
      /// hessian tensor
      typename EvalTraits::HessianTensorType hess_ten;
      /// inverse jacobian matrix
      typename EvalTraits::JacobianInverseType jac_inv;
      /// jacobian matrix
      typename EvalTraits::JacobianMatrixType jac_mat;
      /// image point
      typename EvalTraits::ImagePointType img_point;
      /// domain point
      typename EvalTraits::DomainPointType dom_point;
      /// jacobian determinant
      typename EvalTraits::JacobianDeterminantType jac_det;

      /// our trafo configuration
      static constexpr TrafoTags config = cfg_tags_;
    }; // class EvalData<...>
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_EVAL_DATA_HPP
