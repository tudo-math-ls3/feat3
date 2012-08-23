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
        domain_dim = Shape_::dimension,
        /// image dimension
        image_dim = image_dim_
      };

      /* *************************************************************************************** */

      /// data type
      typedef DataType_ DataType;

      /// trafo coefficient type
      typedef DataType_ TrafoCoeffType;

      /// domain coordinate type
      typedef DataType_ DomainCoordType;
      /// domain point type
      typedef Tiny::Vector<DomainCoordType, domain_dim> DomainPointType;
      /// domain point reference
      typedef DomainPointType& DomainPointRef;
      /// domain point const reference
      typedef const DomainPointType& DomainPointConstRef;

      /// image coordinate type
      typedef DataType_ ImageCoordType;
      /// image point type
      typedef Tiny::Vector<ImageCoordType, image_dim> ImagePointType;
      /// image point reference
      typedef ImagePointType& ImagePointRef;
      /// image point const reference
      typedef const ImagePointType& ImagePointConstRef;

      /// jacobian matrix coefficient type
      typedef DataType_ JacMatCoeff;
      /// jacobian matrix type
      typedef Tiny::Matrix<JacMatCoeff, image_dim, domain_dim> JacMatType;
      /// jacobian matrix reference
      typedef JacMatType& JacMatRef;
      /// jacobian matrix const reference
      typedef const JacMatType& JacMatConstRef;

      /// jacobian inverse matrix coefficient type
      typedef DataType_ JacInvCoeff;
      /// jacobian inverse matrix type
      typedef Tiny::Matrix<JacInvCoeff, domain_dim, image_dim> JacInvType;
      /// jacobian inverse matrix reference
      typedef JacInvType& JacInvRef;
      /// jacobian inverse matrix const reference
      typedef const JacInvType& JacInvConstRef;

      /// jacobian determinant type
      typedef DataType_ JacDetType;

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
        need_jac_det = 0
      };
    }; // struct ConfigBase

    /**
     * \brief Boolean OR operator class template for trafo config tags
     *
     * This templates applies a boolean OR operator on two trafo config tags.
     *
     * \author Peter Zajac
     */
    template<typename Cfg1_, typename Cfg2_>
    struct ConfigOr
    {
      enum
      {
        need_dom_point = int(Cfg1_::need_dom_point) | int(Cfg2_::need_dom_point),
        need_img_point = int(Cfg1_::need_img_point) | int(Cfg2_::need_img_point),
        need_jac_mat = int(Cfg1_::need_jac_mat) | int(Cfg2_::need_jac_mat),
        need_jac_inv = int(Cfg1_::need_jac_inv) | int(Cfg2_::need_jac_inv),
        need_jac_det = int(Cfg1_::need_jac_det) | int(Cfg2_::need_jac_det),
      };
    }; // struct ConfigOr<...>
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_BASE_HPP
