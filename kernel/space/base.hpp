#pragma once
#ifndef KERNEL_SPACE_BASE_HPP
#define KERNEL_SPACE_BASE_HPP 1

// includes, FEAST
#include <kernel/trafo/base.hpp>

namespace FEAST
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

      /// dummy enumeration
      enum
      {
        /// domain dimension
        domain_dim = EvalPolicy::domain_dim,
        /// image dimension
        image_dim = EvalPolicy::image_dim,
        /// maximum number of local dofs
        max_local_dofs = max_local_dofs_
      };

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

      /// dummy enumeration
      enum
      {
        /// domain dimension
        domain_dim = EvalPolicy::domain_dim,
        /// image dimension
        image_dim = EvalPolicy::image_dim,
        /// maximum number of local dofs
        max_local_dofs = max_local_dofs_
      };

      /// basis function value type
      typedef Tiny::Vector<DataType, image_dim> BasisValueType;

      /// basis gradient type
      typedef Tiny::Matrix<DataType, image_dim, image_dim> BasisGradientType;

      /// basis hessian type
      typedef Tiny::Tensor3<DataType, image_dim, image_dim, image_dim> BasisHessianType;
    }; // class StandardScalarEvalTraits<...>

    /**
     * \brief Base class for space config tags
     */
    struct ConfigBase
    {
      /**
       * \brief Space requirements enumeration
       */
      enum SpaceRequirements
      {
        /// specifies whether the space should supply basis function values
        need_value = 0,
        /// specifies whether the space should supply basis function gradients
        need_grad = 0,
        /// specifies whether the space should supply basis function hessians
        need_hess = 0,
        /// specifies whether the space should supply reference basis function values
        need_ref_value = 0,
        /// specifies whether the space should supply reference basis function gradients
        need_ref_grad = 0,
        /// specifies whether the space should supply reference basis function hessians
        need_ref_hess = 0
      };
    }; // struct ConfigBase

    /**
     * \brief Boolean OR operator class template for space config tags
     *
     * This template applies a boolean OR operator on two space config tags.
     *
     * \tparam Cfg1_, Cfg2_
     * The two space config tags that shall be OR'ed.
     *
     * \author Peter Zajac
     */
    template<typename Cfg1_, typename Cfg2_>
    struct ConfigOr
    {
      /** \copydoc Space::ConfigBasae::SpaceRequirements */
      enum SpaceRequirements
      {
        /// specifies whether the space should supply basis function values
        need_value = (Cfg1_::need_value != 0) || (Cfg2_::need_value != 0) ? 1 : 0,
        /// specifies whether the space should supply basis function gradients
        need_grad = (Cfg1_::need_grad != 0) || (Cfg2_::need_grad != 0) ? 1 : 0,
        /// specifies whether the space should supply basis function hessians
        need_hess = (Cfg1_::need_hess != 0) || (Cfg2_::need_hess != 0) ? 1 : 0,
        /// specifies whether the space should supply reference basis function values
        need_ref_value = (Cfg1_::need_ref_value != 0) || (Cfg2_::need_ref_value != 0) ? 1 : 0,
        /// specifies whether the space should supply reference basis function gradients
        need_ref_grad = (Cfg1_::need_ref_grad != 0) || (Cfg2_::need_ref_grad != 0) ? 1 : 0,
        /// specifies whether the space should supply reference basis function hessians
        need_ref_hess = (Cfg1_::need_ref_hess != 0) || (Cfg2_::need_ref_hess != 0) ? 1 : 0
      };
    }; // struct ConfigOr<...>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_BASE_HPP
