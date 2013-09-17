#pragma once
#ifndef KERNEL_SPACE_EVALUATOR_BASE_HPP
#define KERNEL_SPACE_EVALUATOR_BASE_HPP 1

// includes, FEAST
#include <kernel/space/eval_data.hpp>

namespace FEAST
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

      template<typename Cfg_>
      struct ConfigTraits
      {
        /// evaluation data configuration
        typedef Cfg_ EvalDataConfig;

        /// trafo configuration
        typedef Trafo::ConfigBase TrafoConfig;

        /// evaluation data typedef
        typedef Space::EvalData<SpaceEvalTraits, EvalDataConfig> EvalDataType;
      };

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

        /**
         * \brief Basis Hessians capability
         *
         * This entry specifies whether the evaluator is capable of computing basis function hessian matrices.\n
         * If this value is non-zero, the evaluator implements the #eval_hessians member function.\n
         * See #eval_hessians for details.
         */
        can_hess = 0
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
      template<typename SpaceCfg_, typename TrafoCfg_>
      void operator()(
        EvalData<SpaceEvalTraits, SpaceCfg_>& space_data,
        const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const
      {
        // typedef mumbo-jumbo
        typedef EvalData<SpaceEvalTraits, SpaceCfg_> EvalDataType;

        // compute basis values
        Intern::BasisEvalHelper<EvalDataType::have_value != 0>::eval_values(space_data, cast(), trafo_data);
        // compute basis gradients
        Intern::BasisEvalHelper<EvalDataType::have_grad != 0>::eval_gradients(space_data, cast(), trafo_data);
        // compute basis hessians
        Intern::BasisEvalHelper<EvalDataType::have_hess != 0>::eval_hessians(space_data, cast(), trafo_data);
      }

      // Note:
      // The following block serves as an element interface documentation and is therefore only
      // visible to doxygen. The actual functionality has to be supplied by the implementation.
#ifdef DOXYGEN
      /**
       * \brief Evaluates the basis function values on the real cell.
       *
       * \tparam SpaceCfg_, TrafoCfg_
       * The space and trafo configuration classes. These are determined automatically by the compiler.
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
       * \tparam SpaceCfg_, TrafoCfg_
       * The space and trafo configuration classes. These are determined automatically by the compiler.
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

      /**
       * \brief Evaluates the basis function hessians on the real cell.
       *
       * \tparam SpaceCfg_, TrafoCfg_
       * The space and trafo configuration classes. These are determined automatically by the compiler.
       *
       * \param[out] data
       * A reference to an evaluation data structure that shall receive the basis function hessian matrices.
       *
       * \param[in] trafo_data
       * The trafo evaluation data containing information about the evaluation point.
       */
      template<typename SpaceCfg_, typename TrafoCfg_>
      void eval_hessians(
        EvalData<SpaceEvalTraits, SpaceCfg_>& data,
        const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const;
#endif // DOXYGEN
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
          static_assert(evaluator.can_value != 0, "space evaluator does not support basis function values");
          evaluator.eval_values(space_data, trafo_data);
        }

        template<typename SpaceData_, typename Evaluator_, typename TrafoData_>
        static void eval_gradients(SpaceData_& space_data, const Evaluator_& evaluator, const TrafoData_& trafo_data)
        {
          static_assert(evaluator.can_grad != 0, "space evaluator does not support basis function gradients");
          evaluator.eval_gradients(space_data, trafo_data);
        }

        template<typename SpaceData_, typename Evaluator_, typename TrafoData_>
        static void eval_hessians(SpaceData_& space_data, const Evaluator_& evaluator, const TrafoData_& trafo_data)
        {
          static_assert(evaluator.can_hess != 0, "space evaluator does not support basis function hessians");
          evaluator.eval_hessians(space_data, trafo_data);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_EVALUATOR_BASE_HPP
