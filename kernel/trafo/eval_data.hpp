#pragma once
#ifndef KERNEL_TRAFO_EVAL_DATA_HPP
#define KERNEL_TRAFO_EVAL_DATA_HPP 1

// includes, FEAST
#include <kernel/trafo/base.hpp>

namespace FEAST
{
  namespace Trafo
  {
    /// \cond internal
    namespace Intern
    {
      template<typename EvalTraits_, bool need_dom_point_>
      struct DomPointData
      {
      };

      template<typename EvalTraits_, bool need_img_point_>
      struct ImgPointData
      {
      };

      template<typename EvalTraits_, bool need_jac_mat_>
      struct JacMatData
      {
      };

      template<typename EvalTraits_, bool need_jac_inv_>
      struct JacInvData
      {
      };

      template<typename EvalTraits_, bool need_jac_det_>
      struct JacDetData
      {
      };

      template<typename EvalTraits_, bool need_hess_ten_>
      struct HessTenData
      {
      };

      template<typename EvalTraits_, bool need_hess_inv_>
      struct HessInvData
      {
      };

      template<typename EvalTraits_>
      struct DomPointData<EvalTraits_, true>
      {
        /// domain point
        typename EvalTraits_::DomainPointType dom_point;
      };

      template<typename EvalTraits_>
      struct ImgPointData<EvalTraits_, true>
      {
        /// image point
        typename EvalTraits_::ImagePointType img_point;
      };

      template<typename EvalTraits_>
      struct JacMatData<EvalTraits_, true>
      {
        /// jacobian matrix
        typename EvalTraits_::JacobianMatrixType jac_mat;
      };

      template<typename EvalTraits_>
      struct JacInvData<EvalTraits_, true>
      {
        /// inverse jacobian matrix
        typename EvalTraits_::JacobianInverseType jac_inv;
      };

      template<typename EvalTraits_>
      struct JacDetData<EvalTraits_, true>
      {
        /// jacobian determinant
        typename EvalTraits_::JacobianDeterminantType jac_det;
      };

      template<typename EvalTraits_>
      struct HessTenData<EvalTraits_, true>
      {
        /// hessian tensor
        typename EvalTraits_::HessianTensorType hess_ten;
      };

      template<typename EvalTraits_>
      struct HessInvData<EvalTraits_, true>
      {
        /// inverse hessian tensor
        typename EvalTraits_::HessianInverseType hess_inv;
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Trafo evaluation data structure.
     *
     * \tparam EvalTraits_
     * The trafo evaluator traits that this evaluation data shall use.
     *
     * \tparam DataCfg_
     * A trafo data config class that specifies what data shall be supplied.
     *
     * \author Peter Zajac
     */
    template<
      typename EvalTraits_,
      typename DataCfg_>
    class EvalData :
      // Note: The following inheritance list is ordered by size of the objects;
      // this should optimise the alignment of the corresponding subobjects.
      public Intern::HessTenData<EvalTraits_, DataCfg_::need_hess_ten>,
      public Intern::HessInvData<EvalTraits_, DataCfg_::need_hess_inv>,
      public Intern::JacMatData<EvalTraits_, DataCfg_::need_jac_mat>,
      public Intern::JacInvData<EvalTraits_, DataCfg_::need_jac_inv>,
      public Intern::ImgPointData<EvalTraits_, DataCfg_::need_img_point>,
      public Intern::DomPointData<EvalTraits_, DataCfg_::need_dom_point>,
      public Intern::JacDetData<EvalTraits_, DataCfg_::need_jac_det>
    {
    public:
      /// trafo evaluation traits
      typedef EvalTraits_ EvalTraits;

      /// specifies whether domain point coordinates are given
      static constexpr bool have_dom_point = DataCfg_::need_dom_point;
      /// specifies whether image point coordinates are given
      static constexpr bool have_img_point = DataCfg_::need_img_point;
      /// specifies whether the jacobian matrix is given
      static constexpr bool have_jac_mat   = DataCfg_::need_jac_mat;
      /// specifies whether the jacobian inverse matrix is given
      static constexpr bool have_jac_inv   = DataCfg_::need_jac_inv;
      /// specifies whether the jacobian determinant is given
      static constexpr bool have_jac_det   = DataCfg_::need_jac_det;
      /// specifies whether the hessian tensor is given
      static constexpr bool have_hess_ten  = DataCfg_::need_hess_ten;
      /// specifies whether the inverse hessian tensor is given
      static constexpr bool have_hess_inv  = DataCfg_::need_hess_inv;

      /// \cond internal
      typedef Intern::DomPointData<EvalTraits_, have_dom_point> DomPointBase;
      typedef Intern::ImgPointData<EvalTraits_, have_img_point> ImgPointBase;
      typedef Intern::JacMatData<EvalTraits_, have_jac_mat> JacMatBase;
      typedef Intern::JacInvData<EvalTraits_, have_jac_inv> JacInvBase;
      typedef Intern::JacDetData<EvalTraits_, have_jac_det> JacDetBase;
      typedef Intern::HessTenData<EvalTraits_, have_hess_ten> HessTenBase;
      typedef Intern::HessInvData<EvalTraits_, have_hess_inv> HessInvBase;
      /// \endcond

    }; // class EvalData<...>
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_EVAL_DATA_HPP
