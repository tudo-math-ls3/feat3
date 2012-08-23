#pragma once
#ifndef KERNEL_SPACE_BASE_HPP
#define KERNEL_SPACE_BASE_HPP 1

// includes, FEAST
#include <kernel/trafo/base.hpp>

namespace FEAST
{
  /**
   * \brief Element namespace
   */
  namespace Space
  {
    /**
     * \brief Standard scalar evaluator traits class template.
     *
     * This class implements a simple evaluator traits which should suffice for most scalar elements.
     *
     * \tparam EvalPolicy_
     * The evaluation policy that is to be used. See #Trafo::StandardEvalPolicy.
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
    class StandardScalarEvalTraits
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

      /// basis function value coefficient type
      typedef DataType BasisValueCoeff;
      /// basis function value type
      typedef BasisValueCoeff BasisValueType;
      /// basis function value reference
      typedef BasisValueType& BasisValueRef;
      /// basis function value const-reference
      typedef const BasisValueType& BasisValueConstRef;
      /// basis function vector type
      typedef BasisValueType BasisValueVectorType[max_local_dofs];
      /// basis function vector reference
      typedef BasisValueVectorType& BasisValueVectorRef;
      /// basis function value vector const-reference
      typedef const BasisValueVectorType& BasisValueVectorConstRef;

      /// basis gradient coefficient type
      typedef DataType BasisGradientCoeff;
      /// basis gradient type
      typedef Tiny::Vector<BasisGradientCoeff, domain_dim> BasisGradientType;
      /// basis gradient reference
      typedef BasisGradientType& BasisGradientRef;
      /// basis gradient const-reference
      typedef const BasisGradientType& BasisGradientConstRef;
      /// basis gradient vector type
      typedef BasisGradientType BasisGradientVectorType[max_local_dofs];
      /// basis gradient vector reference
      typedef BasisGradientVectorType& BasisGradientVectorRef;
      /// basis function gradient vector cosnt-reference
      typedef const BasisGradientVectorType& BasisGradientVectorConstRef;
    }; // class StandardScalarEvalTraits<...>

    /**
     * \brief Base class for space config tags
     */
    struct ConfigBase
    {
      /// dummy enumeration
      enum
      {
        /// specifies whether the space should supply basis function values
        need_value = 0,
        /// specifies whether the space should supply basis function gradients
        need_grad = 0
      };
    }; // struct ConfigBase

    /**
     * \brief Boolean OR operator class template for space config tags
     *
     * This template applies a boolean OR operator on two space config tags.
     */
    template<typename Cfg1_, typename Cfg2_>
    struct ConfigOr
    {
      enum
      {
        need_value = int(Cfg1_::need_value) | int(Cfg2_::need_value),
        need_grad = int(Cfg1_::need_grad) | int(Cfg2_::need_grad)
      };
    }; // struct ConfigOr<...>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_BASE_HPP
