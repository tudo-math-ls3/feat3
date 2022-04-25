// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_BASE_HPP
#define KERNEL_TRAFO_BASE_HPP 1

// includes, FEAT
#include <kernel/shape.hpp>
#include <kernel/eval_tags.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
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
      /// normal vector type (only for facet trafo)
      typedef Tiny::Vector<DataType, image_dim> NormalVectorType;
    }; // class StandardEvalPolicy<...>
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_BASE_HPP
