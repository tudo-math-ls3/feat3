#pragma once
#ifndef KERNEL_SPACE_EVALUATOR_BASE_HPP
#define KERNEL_SPACE_EVALUATOR_BASE_HPP 1

// includes, FEAT
#include <kernel/space/eval_data.hpp>
#include <kernel/trafo/eval_data.hpp>

namespace FEAT
{
  namespace Space
  {
    /// \cond internal
    namespace Intern
    {
      template<bool _enable>
      struct BasisEvalHelper;
    } // namespace Intern
    /// \endcond

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

      /// trafo domain dimension
      static constexpr int domain_dim = SpaceEvalTraits::domain_dim;
      /// trafo image dimension
      static constexpr int image_dim = SpaceEvalTraits::image_dim;
      /// maximum number of local DOFs
      static constexpr int max_local_dofs = SpaceEvalTraits::max_local_dofs;

      template<SpaceTags cfg_>
      struct ConfigTraits
      {
        /// evaluation data configuration
        static constexpr SpaceTags config = cfg_;

        /// trafo configuration
        static constexpr TrafoTags trafo_config = TrafoTags::none;

        /// evaluation data typedef
        typedef Space::EvalData<SpaceEvalTraits, config> EvalDataType;
      };

#ifdef DOXYGEN
      static constexpr SpaceTags eval_caps = ...;
#endif

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

      /**
       * \brief Space evaluation operator
       *
       * This operator evaluates the basis functions in a point on a mesh cell.
       *
       * \param[out] space_data
       * A reference to the space data that is to be computed.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
      void operator()(
        EvalData<SpaceEvalTraits, space_cfg_>& space_data,
        const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& trafo_data) const
      {
        // typedef mumbo-jumbo
        typedef EvalData<SpaceEvalTraits, space_cfg_> EvalDataType;


        // \compilerhack GCC 6.1 template compiler bug
        // The following static constexpr bool variables are required to
        // circumvent a compiler bug in GCC 6.1
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 70000)
        // compute basis values
        static constexpr bool want_value = *(EvalDataType::config & SpaceTags::value);
        Intern::BasisEvalHelper<want_value>::eval_values(space_data, cast(), trafo_data);
        // compute basis gradients
        static constexpr bool want_grad = *(EvalDataType::config & SpaceTags::grad);
        Intern::BasisEvalHelper<want_grad>::eval_gradients(space_data, cast(), trafo_data);
        // compute basis hessians
        static constexpr bool want_hess = *(EvalDataType::config & SpaceTags::hess);
        Intern::BasisEvalHelper<want_hess>::eval_hessians(space_data, cast(), trafo_data);
#else
        // compute basis values
        Intern::BasisEvalHelper<*(EvalDataType::config & SpaceTags::value)>::eval_values(space_data, cast(), trafo_data);
        // compute basis gradients
        Intern::BasisEvalHelper<*(EvalDataType::config & SpaceTags::grad)>::eval_gradients(space_data, cast(), trafo_data);
        // compute basis hessians
        Intern::BasisEvalHelper<*(EvalDataType::config & SpaceTags::hess)>::eval_hessians(space_data, cast(), trafo_data);
#endif
      }

      /**
       * \brief Evaluates the basis function values on the real cell.
       *
       * \tparam space_cfg_, trafo_cfg_
       * The space and trafo configuration tags. These are determined automatically by the compiler.
       *
       * \param[out] space_data
       * A reference to an evaluation data structure that shall receive the basis function values.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
      void eval_values(
        EvalData<SpaceEvalTraits, space_cfg_>& DOXY(space_data),
        const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& DOXY(trafo_data)) const
      {
        throw InternalError(__func__, __FILE__, __LINE__, "space evaluator does not support basis function values");
      }

      /**
       * \brief Evaluates the basis function gradients on the real cell.
       *
       * \tparam space_cfg_, trafo_cfg_
       * The space and trafo configuration tags. These are determined automatically by the compiler.
       *
       * \param[out] space_data
       * A reference to an evaluation data structure that shall receive the basis function gradients.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
      void eval_gradients(
        EvalData<SpaceEvalTraits, space_cfg_>& DOXY(space_data),
        const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& DOXY(trafo_data)) const
      {
        throw InternalError(__func__, __FILE__, __LINE__, "space evaluator does not support basis function gradients");
      }

      /**
       * \brief Evaluates the basis function hessians on the real cell.
       *
       * \tparam space_cfg_, trafo_cfg_
       * The space and trafo configuration tags. These are determined automatically by the compiler.
       *
       * \param[out] space_data
       * A reference to an evaluation data structure that shall receive the basis function hessian matrices.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
      void eval_hessians(
        EvalData<SpaceEvalTraits, space_cfg_>& DOXY(space_data),
        const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& DOXY(trafo_data)) const
      {
        throw InternalError(__func__, __FILE__, __LINE__, "space evaluator does not support basis function hessians");
      }
    }; // class EvaluatorBase<...>

    /// \cond internal
    namespace Intern
    {
      template<bool _enable>
      struct BasisEvalHelper
      {
        template<typename SpaceData_, typename Evaluator_, typename TrafoData_>
        static void eval_values(SpaceData_&, const Evaluator_&, const TrafoData_&) {}

        template<typename SpaceData_, typename Evaluator_, typename TrafoData_>
        static void eval_gradients(SpaceData_&, const Evaluator_&, const TrafoData_&) {}

        template<typename SpaceData_, typename Evaluator_, typename TrafoData_>
        static void eval_hessians(SpaceData_&, const Evaluator_&, const TrafoData_&) {}
      };

      template<>
      struct BasisEvalHelper<true>
      {
        template<typename SpaceData_, typename Evaluator_, typename TrafoData_>
        static void eval_values(SpaceData_& space_data, const Evaluator_& evaluator, const TrafoData_& trafo_data)
        {
          if(!*(Evaluator_::eval_caps & SpaceTags::value))
            throw InternalError(__func__, __FILE__, __LINE__, "space evaluator does not support basis function values");
          evaluator.eval_values(space_data, trafo_data);
        }

        template<typename SpaceData_, typename Evaluator_, typename TrafoData_>
        static void eval_gradients(SpaceData_& space_data, const Evaluator_& evaluator, const TrafoData_& trafo_data)
        {
          if(!*(Evaluator_::eval_caps & SpaceTags::grad))
            throw InternalError(__func__, __FILE__, __LINE__, "space evaluator does not support basis function gradients");
          evaluator.eval_gradients(space_data, trafo_data);
        }

        template<typename SpaceData_, typename Evaluator_, typename TrafoData_>
        static void eval_hessians(SpaceData_& space_data, const Evaluator_& evaluator, const TrafoData_& trafo_data)
        {
          if(!*(Evaluator_::eval_caps & SpaceTags::hess))
            throw InternalError(__func__, __FILE__, __LINE__, "space evaluator does not support basis function hessians");
          evaluator.eval_hessians(space_data, trafo_data);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_EVALUATOR_BASE_HPP
