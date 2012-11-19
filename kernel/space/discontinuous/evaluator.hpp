#pragma once
#ifndef KERNEL_SPACE_DISCONTINUOUS_EVALUATOR_HPP
#define KERNEL_SPACE_DISCONTINUOUS_EVALUATOR_HPP 1

// includes, FEAST
#include <kernel/space/evaluator_base.hpp>
#include <kernel/space/discontinuous/variant.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Discontinuous
    {
      template<
        typename Space_,
        typename SpaceEvalTraits_,
        typename TrafoEvaluator_,
        typename VariantTag_,
        typename Shape_ = typename Space_::ShapeType>
      class Evaluator DOXY({});

      template<
        typename Space_,
        typename SpaceEvalTraits_,
        typename TrafoEvaluator_,
        typename Shape_>
      class Evaluator<Space_, SpaceEvalTraits_, TrafoEvaluator_, Variant::StdPolyP<0>, Shape_> :
        public EvaluatorBase<
          Evaluator<Space_, SpaceEvalTraits_, TrafoEvaluator_, Variant::StdPolyP<0>, Shape_>,
          TrafoEvaluator_,
          SpaceEvalTraits_>
      {
      public:
        /// basis value coefficient type
        typedef typename SpaceEvalTraits_::BasisValueCoeff BasisValueCoeff;
        /// basis value vector reference
        typedef typename SpaceEvalTraits_::BasisValueVectorRef BasisValueVectorRef;

        /** \copydoc EvaluatorBase::EvaluatorCapabilities */
        enum EvaluatorCapabilities
        {
          /// can compute function values
          can_value = 1
        };

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const Space_& /*space*/)
        {
        }

        /// virtual destructor
        virtual ~Evaluator()
        {
        }

        /**
         * \brief Returns the number of local DOFs.
         *
         * \returns
         * The number of local dofs.
         */
        Index get_num_local_dofs() const
        {
          return 1;
        }

        /** \copydoc Space::EvaluatorBase::eval_values() */
        template<typename TrafoEvalData_>
        void eval_values(
          BasisValueVectorRef values,
          const TrafoEvalData_& /*trafo_data*/) const
        {
          values[0] = BasisValueCoeff(1);
        }
      };
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DISCONTINUOUS_EVALUATOR_HPP
