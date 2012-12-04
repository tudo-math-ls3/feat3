#pragma once
#ifndef KERNEL_SPACE_EVALUATOR_BASE_HPP
#define KERNEL_SPACE_EVALUATOR_BASE_HPP 1

// includes, FEAST
#include <kernel/space/eval_data.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Basic Space Evaluator CRTP base-class template
     *
     * This class is a CRTP base-class used by various finite-element space evaluators to outsource
     * common wrapper code which is independent of the actual finite-element in use.
     *
     * \tparam SpaceEvaluator_
     * The space evaluator class that derives from this class template.
     *
     * \tparam TrafoEvaluator_
     * The trafo evaluator class to be used by this space evaluator.
     *
     * \tparam SpaceEvalTraits_
     * The evaluator traits of the space evaluator.
     *
     * \author Peter Zajac
     */
    template<
      typename SpaceEvaluator_,
      typename TrafoEvaluator_,
      typename SpaceEvalTraits_>
    class EvaluatorBase
    {
    public:
      /// space evaluator type
      typedef SpaceEvaluator_ SpaceEvaluator;

      /// trafo evaluator type
      typedef TrafoEvaluator_ TrafoEvaluator;

      /// space evaluator traits
      typedef SpaceEvalTraits_ SpaceEvalTraits;

      /// evaluation policy
      typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

      /// dummy enum
      enum
      {
        /// trafo domain dimension
        domain_dim = EvalPolicy::domain_dim,
        /// trafo image dimension
        image_dim = EvalPolicy::image_dim,

        /// maximum number of local DOFs
        max_local_dofs = SpaceEvalTraits::max_local_dofs
      };

      /// basis function value coefficient
      typedef typename SpaceEvalTraits::BasisValueCoeff BasisValueCoeff;
      /// basis function value type
      typedef typename SpaceEvalTraits::BasisValueType BasisValueType;
      /// basis function value reference
      typedef typename SpaceEvalTraits::BasisValueRef BasisValueRef;
      /// basis function value const-reference
      typedef typename SpaceEvalTraits::BasisValueConstRef BasisValueConstRef;
      /// basis function value vector type
      typedef typename SpaceEvalTraits::BasisValueVectorType BasisValueVectorType;
      /// basis function value vector reference
      typedef typename SpaceEvalTraits::BasisValueVectorRef BasisValueVectorRef;
      /// basis function value vector const-reference
      typedef typename SpaceEvalTraits::BasisValueVectorConstRef BasisValueVectorConstRef;

      /// basis function gradient coefficient type
      typedef typename SpaceEvalTraits::BasisGradientCoeff BasisGradientCoeff;
      /// basis function gradient type
      typedef typename SpaceEvalTraits::BasisGradientType BasisGradientType;
      /// basis function gradient reference
      typedef typename SpaceEvalTraits::BasisGradientRef BasisGradientRef;
      /// basis function gradient const-reference
      typedef typename SpaceEvalTraits::BasisGradientConstRef BasisGradientConstRef;
      /// basis function gradient vector type
      typedef typename SpaceEvalTraits::BasisGradientVectorType BasisGradientVectorType;
      /// basis function gradient vector reference
      typedef typename SpaceEvalTraits::BasisGradientVectorRef BasisGradientVectorRef;
      /// basis function gradient vector const-reference
      typedef typename SpaceEvalTraits::BasisGradientVectorConstRef BasisGradientVectorConstRef;

      /**
       * \brief Capability enumeration
       *
       * This enumeration specifies the capabilites of the evaluator, i.e. its values specifiy whether
       * the evaluator is capable of computing basis function values, gradients, etc.
       */
      enum EvaluatorCapabilities
      {
        /**
         * \brief Basis Function-Values capability
         *
         * This entry specifies whether the evaluator is capable of computing basis function values.\n
         * If this value is non-zero, the evaluator implements the #eval_values member function.\n
         * See #eval_values for details.
         */
        can_value = 0,

        /**
         * \brief Basis Gradients capability
         *
         * This entry specifies whether the evaluator is capable of computing basis function gradients.\n
         * If this value is non-zero, the evaluator implements the #eval_grads member function.\n
         * See #eval_grads for details.
         */
        can_grad = 0,
      };

    protected:
      /// \cond internal
      SpaceEvaluator& me()
      {
        return static_cast<SpaceEvaluator&>(*this);
      }

      const SpaceEvaluator& me() const
      {
        return static_cast<const SpaceEvaluator&>(*this);
      }
      /// \endcond

    public:
      /**
       * \brief Prepares the evaluator for a given cell.
       *
       * \param[in] trafo_eval
       * A reference to the trafo evaluator containing the cell information.
       */
      void prepare(const TrafoEvaluator& /*trafo_eval*/)
      {
      }

      /**
       * \brief Releases the evaluator from the current cell.
       */
      void finish()
      {
      }

      // Note:
      // The following block serves as an element interface documentation and is therefore only
      // visible to doxygen. The actual functionality has to be supplied by the implementation.
#ifdef DOXYGEN
      /**
       * \brief Evaluates the basis function values on the real cell.
       *
       * \tparam TrafoEvalData_
       * The trafo evaluation data structure type. Is determined automatically by the compiler.\n
       * See TrafoEvalData for details.
       *
       * \param[out] values
       * A reference to a basis value vector that shall receive the basis function values.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<typename TrafoEvalData_>
      void eval_values(
        BasisValueVectorRef values,
        const TrafoEvalData_& trafo_data) const;

      /**
       * \brief Evaluates the basis function gradients on the real cell.
       *
       * \tparam TrafoEvalData_
       * The trafo evaluation data structure type. Is determined automatically by the compiler.\n
       * See TrafoEvalData for details.
       *
       * \param[out] grads
       * A reference to a basis gradient vector that shall receive the basis function gradients.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<typename TrafoEvalData_>
      void eval_gradients(
        BasisGradientVectorRef grads,
        const TrafoEvalData_& trafo_data) const;
#endif // DOXYGEN
    }; // class EvaluatorBase<...>

    /**
     * \brief Finite-Element Parametric Evaluator CRTP base-class template
     *
     * This class is a CRTP base-class for parametric Finite-Element space evaluators.
     *
     * \tparam SpaceEvaluator_
     * The space evaluator class which derives from this class template.
     *
     * \tparam TrafoEvaluator_
     * The trafo evaluator class to be used by the space evaluator.
     *
     * \tparam SpaceEvalTraits_
     * The space evaluator traits of the space evaluator.
     *
     * \tparam ReferenceCapabilities_
     * A tag class containing the capabilities of the space evaluator for reference value computation.
     *
     * \author Peter Zajac
     */
    template<
      typename SpaceEvaluator_,
      typename TrafoEvaluator_,
      typename SpaceEvalTraits_,
      typename ReferenceCapabilities_>
    class EvaluatorParametric :
      public EvaluatorBase<SpaceEvaluator_, TrafoEvaluator_, SpaceEvalTraits_>
    {
    public:
      /// base class typedef
      typedef EvaluatorBase<SpaceEvaluator_, TrafoEvaluator_, SpaceEvalTraits_> BaseClass;

      /** \copydoc EvaluatorBase::EvaluatorCapabilities */
      enum EvaluatorCapabilities
      {
        /// can compute function values if the trafo can compute domain points
        can_value =
          int(ReferenceCapabilities_::can_ref_value) |
          int(TrafoEvaluator_::can_dom_point),

        /// can compute gradients if the trafo can compute jacobian matrices and domain points
        can_grad =
          int(ReferenceCapabilities_::can_ref_grad) |
          int(TrafoEvaluator_::can_dom_point) |
          int(TrafoEvaluator_::can_jac_inv)
      };

      /// basis function value vector reference
      typedef typename SpaceEvalTraits_::BasisValueVectorRef BasisValueVectorRef;
      /// basis function gradient vector reference
      typedef typename SpaceEvalTraits_::BasisGradientVectorRef BasisGradientVectorRef;

    protected:
      /// \cond internal
      SpaceEvaluator_& me()
      {
        return static_cast<SpaceEvaluator_&>(*this);
      }

      const SpaceEvaluator_& me() const
      {
        return static_cast<const SpaceEvaluator_&>(*this);
      }
      /// \endcond

    private:
      /// auxiliary basis gradient vector
      mutable typename SpaceEvalTraits_::BasisGradientVectorType _aux_grads;

    public:
      /** \copydoc EvaluatorBase::eval_values() */
      template<typename TrafoEvalData_>
      void eval_values(
        BasisValueVectorRef values,
        const TrafoEvalData_& trafo_data) const
      {
        static_assert(can_value != 0, "space evaluator can't compute function values");

        // compute reference values in the domain point;
        // these coincide with the real values in the image point
        me().eval_ref_values(values, trafo_data.dom_point);
      }

      /** \copydoc EvaluatorBase::eval_grads() */
      template<typename TrafoEvalData_>
      void eval_gradients(
        BasisGradientVectorRef grads,
        const TrafoEvalData_& trafo_data) const
      {
        static_assert(can_grad != 0, "space evaluator can't compute gradients");

        // compute reference gradients in the domain point and store them in the
        // auxiliary gradient vector
        me().eval_ref_gradients(_aux_grads, trafo_data.dom_point);

        // loop over all basis functions
        for(Index j(0); j < me().get_num_local_dofs(); ++j)
        {
          // apply the chain rule on the reference gradient, i.e.
          // mutliply it by the transposed jacobian matrix inverse
          grads[j].set_vec_mat_mult(_aux_grads[j], trafo_data.jac_inv);
        }
      }
    }; // class EvaluatorParametric<...>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_EVALUATOR_BASE_HPP
