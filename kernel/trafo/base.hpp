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

      /// domain dimension
      static constexpr int domain_dim = Shape_::dimension > 0 ? Shape_::dimension : 1;
      /// image dimension
      static constexpr int image_dim = image_dim_;

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
      /// specifies whether the trafo should supply domain point coordinates
      static constexpr bool need_dom_point = false;
      /// specifies whether the trafo should supply image point coordinates
      static constexpr bool need_img_point = false;
      /// specifies whether the trafo should supply jacobian matrices
      static constexpr bool need_jac_mat   = false;
      /// specifies whether the trafo should supply inverse jacobian matrices
      static constexpr bool need_jac_inv   = false;
      /// specifies whether the trafo should supply jacobian determinants
      static constexpr bool need_jac_det   = false;
      /// specifies whether the trafo should supply hessian tensors
      static constexpr bool need_hess_ten  = false;
      /// specifies whether the trafo should supply inverse hessian tensors
      static constexpr bool need_hess_inv  = false;
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
      /// specifies whether the trafo should supply domain point coordinates
      static constexpr bool need_dom_point = Cfg1_::need_dom_point || Cfg2_::need_dom_point;
      /// specifies whether the trafo should supply image point coordinates
      static constexpr bool need_img_point = Cfg1_::need_img_point || Cfg2_::need_img_point;
      /// specifies whether the trafo should supply jacobian matrices
      static constexpr bool need_jac_mat   = Cfg1_::need_jac_mat   || Cfg2_::need_jac_mat;
      /// specifies whether the trafo should supply inverse jacobian matrices
      static constexpr bool need_jac_inv   = Cfg1_::need_jac_inv   || Cfg2_::need_jac_inv;
      /// specifies whether the trafo should supply jacobian determinants
      static constexpr bool need_jac_det   = Cfg1_::need_jac_det   || Cfg2_::need_jac_det;
      /// specifies whether the trafo should supply hessian tensors
      static constexpr bool need_hess_ten  = Cfg1_::need_hess_ten  || Cfg2_::need_hess_ten;
      /// specifies whether the trafo should supply inverse hessian tensors
      static constexpr bool need_hess_inv  = Cfg1_::need_hess_inv  || Cfg2_::need_hess_inv;
    }; // struct ConfigOr<...>
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_BASE_HPP
