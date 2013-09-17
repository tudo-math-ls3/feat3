#pragma once
#ifndef KERNEL_TRAFO_BASE_HPP
#define KERNEL_TRAFO_BASE_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAST
{
  /**
   * \brief Transformation namespace
   *
   * This namespace encapsulates classes and class templates related to the transformation between reference
   * cells and the cells of a Geometry mesh object, which is a basic requirement for the definition of
   * finite element spaces.
   */
  namespace Trafo
  {
    /**
     * \brief Standard evaluation policy class template
     *
     * \tparam Shape_
     * The shape for which the evaluator shall be used.
     *
     * \tparam DataType_
     * The data-type that is to be used.
     *
     * \tparam image_dim_
     * The image dimension.
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename DataType_ = Real,
      int image_dim_ = Shape_::dimension>
    struct StandardEvalPolicy
    {
      /// shape type
      typedef Shape_ ShapeType;

      /// dummy enumeration
      enum
      {
        /// domain dimension
        domain_dim = Shape_::dimension > 0 ? Shape_::dimension : 1,
        /// image dimension
        image_dim = image_dim_
      };

      /* *************************************************************************************** */

      /// data type
      typedef DataType_ DataType;

      /// domain point type
      typedef Tiny::Vector<DataType, domain_dim> DomainPointType;
      /// image point type
      typedef Tiny::Vector<DataType, image_dim> ImagePointType;
      /// jacobian matrix type
      typedef Tiny::Matrix<DataType, image_dim, domain_dim> JacobianMatrixType;
      /// inverse jacobian matrix type
      typedef Tiny::Matrix<DataType, domain_dim, image_dim> JacobianInverseType;
      /// jacobian determinant type
      typedef DataType JacobianDeterminantType;
      /// hessian tensor type
      typedef Tiny::Tensor3<DataType, image_dim, domain_dim, domain_dim> HessianTensorType;
      /// inverse hessian tensor type
      typedef Tiny::Tensor3<DataType, image_dim, domain_dim, domain_dim> HessianInverseType;

    }; // class StandardEvalPolicy<...>

    /**
     * \brief Base class for trafo config tags
     *
     * \author Peter Zajac
     */
    struct ConfigBase
    {
      /**
       * \brief Trafo requirements enumeration
       */
      enum TrafoRequirements
      {
        /// specifies whether the trafo should supply domain point coordinates
        need_dom_point = 0,
        /// specifies whether the trafo should supply image point coordinates
        need_img_point = 0,
        /// specifies whether the trafo should supply jacobian matrices
        need_jac_mat = 0,
        /// specifies whether the trafo should supply inverse jacobian matrices
        need_jac_inv = 0,
        /// specifies whether the trafo should supply jacobian determinants
        need_jac_det = 0,
        /// specifies whether the trafo should supply hessian tensors
        need_hess_ten = 0,
        /// specifies whether the trafo should supply inverse hessian tensors
        need_hess_inv = 0
      };
    }; // struct ConfigBase

    /**
     * \brief Boolean OR operator class template for trafo config tags
     *
     * This templates applies a boolean OR operator on two trafo config tags.
     *
     * \tparam Cfg1_, Cfg2_
     * The two trafo config tags that shall be OR'ed.
     *
     * \author Peter Zajac
     */
    template<typename Cfg1_, typename Cfg2_>
    struct ConfigOr
    {
      /** \copydoc Trafo::ConfigBase::TrafoRequirements */
      enum TrafoRequirements
      {
        /// specifies whether the trafo should supply domain point coordinates
        need_dom_point = (Cfg1_::need_dom_point != 0) || (Cfg2_::need_dom_point != 0) ? 1 : 0,
        /// specifies whether the trafo should supply image point coordinates
        need_img_point = (Cfg1_::need_img_point != 0) | (Cfg2_::need_img_point != 0) ? 1 : 0,
        /// specifies whether the trafo should supply jacobian matrices
        need_jac_mat = (Cfg1_::need_jac_mat != 0) | (Cfg2_::need_jac_mat != 0) ? 1 : 0,
        /// specifies whether the trafo should supply inverse jacobian matrices
        need_jac_inv = (Cfg1_::need_jac_inv != 0) | (Cfg2_::need_jac_inv != 0) ? 1 : 0,
        /// specifies whether the trafo should supply jacobian determinants
        need_jac_det = (Cfg1_::need_jac_det != 0) | (Cfg2_::need_jac_det != 0) ? 1 : 0,
        /// specifies whether the trafo should supply hessian tensors
        need_hess_ten = (Cfg1_::need_hess_ten != 0) | (Cfg2_::need_hess_ten != 0) ? 1 : 0,
        /// specifies whether the trafo should supply inverse hessian tensors
        need_hess_inv = (Cfg1_::need_hess_inv != 0) | (Cfg2_::need_hess_inv != 0) ? 1 : 0
      };
    }; // struct ConfigOr<...>
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_BASE_HPP
