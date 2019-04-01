// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_BASE_HPP
#define KERNEL_SPACE_BASE_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/eval_tags.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  /**
   * \brief Finite Element Space namespace
   *
   * This namespaces encapsulates classes and class templates which implement Finite Element spaces.
   */
  namespace Space
  {
    /**
     * \brief Standard scalar evaluator traits class template.
     *
     * This class implements a simple evaluator traits which should suffice for most scalar elements.
     *
     * \tparam EvalPolicy_
     * The evaluation policy that is to be used. See Trafo::StandardEvalPolicy.
     *
     * \tparam max_local_dofs_
     * The maximum number of local degrees of freedom.
     *
     * \tparam DataType_
     * The data type that is to be used.
     *
     * \author Peter Zajac
     */
    template<
      typename EvalPolicy_,
      int max_local_dofs_,
      typename DataType_ = typename EvalPolicy_::DataType>
    class StandardScalarEvalTraits :
      public EvalPolicy_
    {
    public:
      /// evaluation policy typedef
      typedef EvalPolicy_ EvalPolicy;

      /// evaluation data type
      typedef DataType_ DataType;

      /// domain dimension
      static constexpr int domain_dim = EvalPolicy::domain_dim;
      /// image dimension
      static constexpr int image_dim = EvalPolicy::image_dim;
      /// maximum number of local dofs
      static constexpr int max_local_dofs = max_local_dofs_;

      /// basis function value type
      typedef DataType BasisValueType;
      /// basis function value type on reference element
      typedef DataType BasisReferenceValueType;

      /// basis gradient type
      typedef Tiny::Vector<DataType, image_dim> BasisGradientType;
      /// basis gradient type on reference element
      typedef Tiny::Vector<DataType, domain_dim> BasisReferenceGradientType;

      /// basis hessian matrix type
      typedef Tiny::Matrix<DataType, image_dim, image_dim> BasisHessianType;
      /// basis hessian type on reference element
      typedef Tiny::Matrix<DataType, domain_dim, domain_dim> BasisReferenceHessianType;
    }; // class StandardScalarEvalTraits<...>

    /**
     * \brief Standard vector evaluator traits class template.
     *
     * This class implements a simple evaluator traits which should suffice for most vector-valued elements.
     *
     * \tparam EvalPolicy_
     * The evaluation policy that is to be used. See Trafo::StandardEvalPolicy.
     *
     * \tparam max_local_dofs_
     * The maximum number of local degrees of freedom.
     *
     * \tparam DataType_
     * The data type that is to be used.
     *
     * \author Peter Zajac
     */
    template<
      typename EvalPolicy_,
      int max_local_dofs_,
      typename DataType_ = typename EvalPolicy_::DataType>
    class StandardVectorEvalTraits :
      public EvalPolicy_
    {
    public:
      /// evaluation policy typedef
      typedef EvalPolicy_ EvalPolicy;

      /// evaluation data type
      typedef DataType_ DataType;

      /// domain dimension
      static constexpr int domain_dim = EvalPolicy::domain_dim;
      /// image dimension
      static constexpr int image_dim = EvalPolicy::image_dim;
      /// maximum number of local dofs
      static constexpr int max_local_dofs = max_local_dofs_;

      /// basis function value type
      typedef Tiny::Vector<DataType, image_dim> BasisValueType;

      /// basis gradient type
      typedef Tiny::Matrix<DataType, image_dim, image_dim> BasisGradientType;

      /// basis hessian type
      typedef Tiny::Tensor3<DataType, image_dim, image_dim, image_dim> BasisHessianType;
    }; // class StandardScalarEvalTraits<...>
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_BASE_HPP
