// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_PARAMETRIC_EVALUATOR_HPP
#define KERNEL_SPACE_PARAMETRIC_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/space/evaluator_base.hpp>

namespace FEAT
{
  namespace Space
  {
    /// \cond internal
    namespace Intern
    {
      template<bool _enable>
      struct ParamBasisEvalHelper;
    } // namespace Intern
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
     * \tparam ref_caps_
     * A tag class containing the capabilities of the space evaluator for reference value computation.
     *
     * \author Peter Zajac
     */
    template<
      typename SpaceEvaluator_,
      typename TrafoEvaluator_,
      typename SpaceEvalTraits_,
      SpaceTags ref_caps_>
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

      /// the capabilities of the trafo evaluator
      static constexpr TrafoTags trafo_caps = TrafoEvaluator::eval_caps;

      /// our evaluation capabilities
      static constexpr SpaceTags eval_caps =
        // reference values, gradients and hessians are supplied by our derived CRTP class
        // (if available), so we mask that out from the reference capabilities
        (ref_caps_ & (SpaceTags::ref_value | SpaceTags::ref_grad | SpaceTags::ref_hess)) |
        // we supply values if the derived class supplies reference values and the
        // trafo supplies domain point coordinates:
        ((*(ref_caps_ & SpaceTags::ref_value) && *(trafo_caps & TrafoTags::dom_point)) ?
          SpaceTags::value : SpaceTags::none) |
        // we supply gradients if the derived class supplies reference gradients and
        // the trafo supplies domain point coordinates and jacobian inverse matrices
        // (which we require for the first-order chain rule)
        ((*(ref_caps_ & SpaceTags::ref_grad) && *(trafo_caps & (TrafoTags::dom_point | TrafoTags::jac_inv))) ?
          SpaceTags::grad : SpaceTags::none) |
        // we supply hessians if the derived class supplies reference hessians and
        // the trafo supplies domain point coordinates, jacobian inverse matrices and
        // hessian inverse tensors (which we require for the second-order chain rule)
        ((*(ref_caps_ & SpaceTags::ref_hess) && *(trafo_caps & (TrafoTags::dom_point | TrafoTags::jac_inv | TrafoTags::hess_inv))) ?
          SpaceTags::hess : SpaceTags::none);

      /// domain dimension
      static constexpr int domain_dim = SpaceEvalTraits::domain_dim;
      /// image dimension
      static constexpr int image_dim = SpaceEvalTraits::image_dim;
      /// maximum number of local dofs
      static constexpr int max_local_dofs = SpaceEvalTraits::max_local_dofs;

      /**
       * \brief Space configuration traits class template.
       */
      template<SpaceTags cfg_>
      struct ConfigTraits
      {
        /// space evaluation config
        static constexpr SpaceTags config =
          (  cfg_ & (SpaceTags::value | SpaceTags::grad | SpaceTags::hess)) |
          // we need reference hessians if hessians are required
          (*(cfg_ & (SpaceTags::ref_hess | SpaceTags::hess)) ? SpaceTags::ref_hess : SpaceTags::none) |
          // we need reference gradients if gradients or hessians are required
          (*(cfg_ & (SpaceTags::ref_grad | SpaceTags::grad | SpaceTags::hess)) ? SpaceTags::ref_grad : SpaceTags::none) |
          // we need reference values if values are required
          (*(cfg_ & (SpaceTags::ref_value | SpaceTags::value)) ? SpaceTags::ref_value : SpaceTags::none);

        /// trafo evaluation configuration
        static constexpr TrafoTags trafo_config = TrafoTags::dom_point |
          // we need jacobian inverse matrices for gradient or hessian computation by the chain rule
          (*(config & (SpaceTags::grad | SpaceTags::hess)) ? TrafoTags::jac_inv : TrafoTags::none) |
          // we need hessan inverse tensors for hessian computation by the chain rule
          (*(config & SpaceTags::hess) ? TrafoTags::hess_inv : TrafoTags::none);

        /// evaluation data typedef
        typedef Space::EvalData<SpaceEvalTraits, config> EvalDataType;
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
      template<SpaceTags space_cfg_>
      void reference_eval(
        EvalData<SpaceEvalTraits, space_cfg_>& space_ref_data,
        const typename SpaceEvalTraits::DomainPointType& dom_point) const
      {
        // compute reference basis values
        Intern::ParamBasisEvalHelper<*(space_cfg_ & SpaceTags::ref_value)>::eval_ref_values(space_ref_data, cast(), dom_point);
        // compute reference basis values
        Intern::ParamBasisEvalHelper<*(space_cfg_ & SpaceTags::ref_grad)>::eval_ref_gradients(space_ref_data, cast(), dom_point);
        // compute reference basis values
        Intern::ParamBasisEvalHelper<*(space_cfg_ & SpaceTags::ref_hess)>::eval_ref_hessians(space_ref_data, cast(), dom_point);
      }

      /** \copydoc EvaluatorBase::operator()() */
      template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
      void operator()(
        EvalData<SpaceEvalTraits, space_cfg_>& space_data,
        const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& trafo_data) const
      {
        // evaluate reference data
        cast().reference_eval(space_data, trafo_data.dom_point);

        // transform basis values
        Intern::ParamBasisEvalHelper<*(space_cfg_ & SpaceTags::value)>::trans_values(space_data, trafo_data);
        // transform basis gradients
        Intern::ParamBasisEvalHelper<*(space_cfg_ & SpaceTags::grad)>::trans_gradients(space_data, trafo_data);
        // transform basis hessians
        Intern::ParamBasisEvalHelper<*(space_cfg_ & SpaceTags::hess)>::trans_hessians(space_data, trafo_data);
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
    } // namespace Intern
    /// \endcond
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_PARAMETRIC_EVALUATOR_HPP
