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

      /// trafo evaluator traits
      typedef typename TrafoEvaluator::EvalTraits TrafoEvalTraits;

      /// dummy enum
      enum
      {
        /// trafo domain dimension
        domain_dim = SpaceEvalTraits::domain_dim,
        /// trafo image dimension
        image_dim = SpaceEvalTraits::image_dim,

        /// maximum number of local DOFs
        max_local_dofs = SpaceEvalTraits::max_local_dofs
      };

      /// basis function value coefficient
      typedef typename SpaceEvalTraits::BasisValueCoeff BasisValueCoeff;
      /// basis function value type
      typedef typename SpaceEvalTraits::BasisValueType BasisValueType;

      /// basis function gradient coefficient type
      typedef typename SpaceEvalTraits::BasisGradientCoeff BasisGradientCoeff;
      /// basis function gradient type
      typedef typename SpaceEvalTraits::BasisGradientType BasisGradientType;

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
         * If this value is non-zero, the evaluator implements the #eval_gradients member function.\n
         * See #eval_gradients for details.
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
      void prepare(const TrafoEvaluator& DOXY(trafo_eval))
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
       * \param[out] data
       * A reference to an evaluation data structure that shall receive the basis function values.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<typename SpaceCfg_, typename TrafoCfg_>
      void eval_values(
        EvalData<SpaceEvalTraits, SpaceCfg_>& data,
        const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const

      /**
       * \brief Evaluates the basis function gradients on the real cell.
       *
       * \tparam TrafoEvalData_
       * The trafo evaluation data structure type. Is determined automatically by the compiler.\n
       * See TrafoEvalData for details.
       *
       * \param[out] data
       * A reference to an evaluation data structure that shall receive the basis function gradients.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<typename SpaceCfg_, typename TrafoCfg_>
      void eval_gradients(
        EvalData<SpaceEvalTraits, SpaceCfg_>& data,
        const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const;

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

      /// space evaluator type
      typedef SpaceEvaluator_ SpaceEvaluator;

      /// trafo evaluator type
      typedef TrafoEvaluator_ TrafoEvaluator;

      /// space evaluator traits
      typedef SpaceEvalTraits_ SpaceEvalTraits;

      /// trafo evaluator traits
      typedef typename TrafoEvaluator::EvalTraits TrafoEvalTraits;

      /** \copydoc EvaluatorBase::EvaluatorCapabilities */
      enum EvaluatorCapabilities
      {
        /// can compute function values if the trafo can compute domain points
        can_value =
          (ReferenceCapabilities_::can_ref_value != 0) &&
          (TrafoEvaluator_::can_dom_point != 0) ? 1 : 0,

        /// can compute gradients if the trafo can compute jacobian matrices and domain points
        can_grad =
          (ReferenceCapabilities_::can_ref_grad != 0) &&
          (TrafoEvaluator_::can_dom_point != 0) &&
          (TrafoEvaluator_::can_jac_inv != 0) ? 1 : 0
      };

      /// dummy enumeration
      enum
      {
        /// domain dimension
        domain_dim = SpaceEvalTraits::domain_dim,
        /// image dimension
        image_dim = SpaceEvalTraits::image_dim,
        /// maximum number of local dofs
        max_local_dofs = SpaceEvalTraits::max_local_dofs
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

      /// \cond internal
      struct ReferenceGradientData
      {
        struct PhiData
        {
          Tiny::Vector<typename SpaceEvalTraits::BasisGradientCoeff, domain_dim> grad;
        } phi[max_local_dofs];
      };
      /// \endcond

    public:
      /** \copydoc EvaluatorBase::eval_values() */
      template<typename SpaceCfg_, typename TrafoCfg_>
      void eval_values(
        EvalData<SpaceEvalTraits, SpaceCfg_>& data,
        const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const
      {
        static_assert(can_value != 0, "space evaluator can't compute function values");
        static_assert(trafo_data.have_dom_point != 0, "trafo data does not offer domain point coordinates");

        // compute reference values in the domain point;
        // these coincide with the real values in the image point
        me().eval_ref_values(data, trafo_data.dom_point);
      }

      /** \copydoc EvaluatorBase::eval_gradients() */
      template<typename SpaceCfg_, typename TrafoCfg_>
      void eval_gradients(
        EvalData<SpaceEvalTraits, SpaceCfg_>& data,
        const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const
      {
        static_assert(can_grad != 0, "space evaluator can't compute function gradients");
        static_assert(trafo_data.have_dom_point != 0, "trafo data does not offer domain point coordinates");
        static_assert(trafo_data.have_jac_inv != 0, "trafo data does not offer inverse jacobian matrix");

        // declare auxiliary reference gradient data
        ReferenceGradientData aux_grads;

        // compute reference gradients in the domain point and store them in the
        // auxiliary gradient vector
        me().eval_ref_gradients(aux_grads, trafo_data.dom_point);

        // loop over all basis functions
        for(Index j(0); j < me().get_num_local_dofs(); ++j)
        {
          // apply the chain rule on the reference gradient, i.e.
          // mutliply it by the transposed jacobian matrix inverse
          data.phi[j].grad.set_vec_mat_mult(aux_grads.phi[j].grad, trafo_data.jac_inv);
        }
      }
    }; // class EvaluatorParametric<...>
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_EVALUATOR_BASE_HPP
