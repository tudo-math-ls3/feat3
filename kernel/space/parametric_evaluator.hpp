#pragma once
#ifndef KERNEL_SPACE_PARAMETRIC_EVALUATOR_HPP
#define KERNEL_SPACE_PARAMETRIC_EVALUATOR_HPP 1

// includes, FEAST
#include <kernel/space/evaluator_base.hpp>

namespace FEAST
{
  namespace Space
  {
    /// \cond internal
    namespace Intern
    {
      template<bool _enable>
      struct ParamBasisEvalHelper;
    }; // namespace Intern
    /// \endcond

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
    class ParametricEvaluator :
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
          (TrafoEvaluator_::can_jac_inv != 0) ? 1 : 0,

        /// can compute hessians if the trafo can compute inverse hessian tensors,
        /// inverse jacobian matrices and domain points
        can_hess =
          (ReferenceCapabilities_::can_ref_grad != 0) &&
          (ReferenceCapabilities_::can_ref_hess != 0) &&
          (TrafoEvaluator_::can_dom_point != 0) &&
          (TrafoEvaluator_::can_hess_inv != 0) &&
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

      /**
       * \brief Space configuration traits class template.
       */
      template<typename Cfg_>
      struct ConfigTraits
      {
        /// evaluation data configuration
        struct EvalDataConfig :
          public Space::ConfigBase
        {
          /** \copydoc Space::ConfigBase::SpaceRequirements */
          enum
          {
            /// need hessians if they are desired
            need_hess = (Cfg_::need_hess != 0) ? 1 : 0,
            /// need gradients if they are desired
            need_grad = (Cfg_::need_grad != 0) ? 1 : 0,
            /// need values if they are desired
            need_value = (Cfg_::need_value != 0) ? 1 : 0,
            /// need reference hessians for hessians
            need_ref_hess = (Cfg_::need_ref_hess != 0) || (need_hess != 0) ? 1 : 0,
            /// need reference gradients for gradients and hessians
            need_ref_grad = (Cfg_::need_ref_grad != 0) || (need_grad != 0) || (need_hess != 0) ? 1 : 0,
            /// need reference values for values
            need_ref_value = (Cfg_::need_ref_value != 0) || (need_value != 0) ? 1 : 0
          };
        };

        /// trafo configuration
        struct TrafoConfig :
          public Trafo::ConfigBase
        {
          /** \copydoc Trafo::ConfigBase::TrafoRequirements */
          enum
          {
            /// we always need domain point coordinates
            need_dom_point = 1,
            /// we need inverse jacobians for basis gradients and hessians
            need_jac_inv = (EvalDataConfig::need_grad != 0) ? 1 : 0,
            /// we need inverse hessians for basis hessians
            need_hess_inv = (EvalDataConfig::need_hess != 0) ? 1 : 0
          };
        };

        /// evaluation data typedef
        typedef Space::EvalData<SpaceEvalTraits, EvalDataConfig> EvalDataType;
      };

    protected:
      /// \cond internal
      /// CRTP downcast
      SpaceEvaluator& cast()
      {
        return static_cast<SpaceEvaluator&>(*this);
      }

      /// CRTP downcast
      const SpaceEvaluator& cast() const
      {
        return static_cast<const SpaceEvaluator&>(*this);
      }
      /// \endcond

    public:
      /**
       * \brief Space reference evaluation function
       *
       * This function evaluates the basis functions in a point on the reference cell.
       *
       * \param[out] space_ref_data
       * A reference to the space reference data that is to be computed.
       *
       * \param[in] dom_point
       * A reference to the domain point in which the evaluation shall take place.
       */
      template<typename SpaceCfg_>
      void reference_eval(
        EvalData<SpaceEvalTraits, SpaceCfg_>& space_ref_data,
        const typename SpaceEvalTraits::DomainPointType& dom_point) const
      {
        // compute reference basis values
        Intern::ParamBasisEvalHelper<SpaceCfg_::need_ref_value != 0>::eval_ref_values(space_ref_data, cast(), dom_point);
        // compute reference basis values
        Intern::ParamBasisEvalHelper<SpaceCfg_::need_ref_grad != 0>::eval_ref_gradients(space_ref_data, cast(), dom_point);
        // compute reference basis values
        Intern::ParamBasisEvalHelper<SpaceCfg_::need_ref_hess != 0>::eval_ref_hessians(space_ref_data, cast(), dom_point);
      }

      /** \copydoc EvaluatorBase::operator()() */
      template<typename SpaceCfg_, typename TrafoCfg_>
      void operator()(
        EvalData<SpaceEvalTraits, SpaceCfg_>& space_data,
        const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const
      {
        // evaluate reference data
        cast().reference_eval(space_data, trafo_data.dom_point);

        // transform basis values
        Intern::ParamBasisEvalHelper<SpaceCfg_::need_value != 0>::trans_values(space_data, trafo_data);
        // transform basis gradients
        Intern::ParamBasisEvalHelper<SpaceCfg_::need_grad != 0>::trans_gradients(space_data, trafo_data);
        // transform basis hessians
        Intern::ParamBasisEvalHelper<SpaceCfg_::need_hess != 0>::trans_hessians(space_data, trafo_data);
      }
    }; // class EvaluatorParametric<...>

    /// \cond internal
    namespace Intern
    {
      template<bool _enable>
      struct ParamBasisEvalHelper
      {
        template<typename SpaceData_, typename Evaluator_, typename DomPoint_>
        static void eval_ref_values(SpaceData_&, const Evaluator_&, const DomPoint_&) {}

        template<typename SpaceData_, typename Evaluator_, typename DomPoint_>
        static void eval_ref_gradients(SpaceData_&, const Evaluator_&, const DomPoint_&) {}

        template<typename SpaceData_, typename Evaluator_, typename DomPoint_>
        static void eval_ref_hessians(SpaceData_&, const Evaluator_&, const DomPoint_&) {}

        template<typename SpaceData_, typename TrafoData_>
        static void trans_values(SpaceData_&, const TrafoData_&) {}

        template<typename SpaceData_, typename TrafoData_>
        static void trans_gradients(SpaceData_&, const TrafoData_&) {}

        template<typename SpaceData_, typename TrafoData_>
        static void trans_hessians(SpaceData_&, const TrafoData_&) {}
      };

      template<>
      struct ParamBasisEvalHelper<true>
      {
        template<typename SpaceData_, typename Evaluator_, typename DomPoint_>
        static void eval_ref_values(SpaceData_& space_data, const Evaluator_& evaluator, const DomPoint_& dom_point)
        {
          evaluator.eval_ref_values(space_data, dom_point);
        }

        template<typename SpaceData_, typename Evaluator_, typename DomPoint_>
        static void eval_ref_gradients(SpaceData_& space_data, const Evaluator_& evaluator, const DomPoint_& dom_point)
        {
          evaluator.eval_ref_gradients(space_data, dom_point);
        }

        template<typename SpaceData_, typename Evaluator_, typename DomPoint_>
        static void eval_ref_hessians(SpaceData_& space_data, const Evaluator_& evaluator, const DomPoint_& dom_point)
        {
          evaluator.eval_ref_hessians(space_data, dom_point);
        }

        template<typename SpaceData_, typename TrafoData_>
        static void trans_values(SpaceData_& space_data, const TrafoData_&)
        {
          // loop over all basis functions
          for(int i(0); i < SpaceData_::max_local_dofs; ++i)
          {
            // and copy the basis function value
            space_data.phi[i].value = space_data.phi[i].ref_value;
          }
        }

        template<typename SpaceData_, typename TrafoData_>
        static void trans_gradients(SpaceData_& space_data, const TrafoData_& trafo_data)
        {
          // loop over all basis functions
          for(int i(0); i < SpaceData_::max_local_dofs; ++i)
          {
            // and apply the first-order chain rule
            space_data.phi[i].grad.set_vec_mat_mult(space_data.phi[i].ref_grad, trafo_data.jac_inv);
          }
        }

        template<typename SpaceData_, typename TrafoData_>
        static void trans_hessians(SpaceData_& space_data, const TrafoData_& trafo_data)
        {
          // loop over all basis functions
          for(int i(0); i < SpaceData_::max_local_dofs; ++i)
          {
            // and apply the second-order chain rule
            space_data.phi[i].hess.set_double_mat_mult(space_data.phi[i].ref_hess, trafo_data.jac_inv, trafo_data.jac_inv);
            space_data.phi[i].hess.add_vec_tensor_mult(space_data.phi[i].ref_grad, trafo_data.hess_inv);
          }
        }
      };
    }; // namespace Intern
    /// \endcond
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_PARAMETRIC_EVALUATOR_HPP
