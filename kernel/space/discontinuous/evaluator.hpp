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
        /// space evaluation traits
        typedef SpaceEvalTraits_ SpaceEvalTraits;

        /// trafo evaluation traits
        typedef typename TrafoEvaluator_::EvalTraits TrafoEvalTraits;

        /// basis value coefficient type
        typedef typename SpaceEvalTraits::DataType DataType;

        /// can compute function values
        static constexpr bool can_value = true;

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const Space_& DOXY(space))
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
        template<typename SpaceCfg_, typename TrafoCfg_>
        void eval_values(
          EvalData<SpaceEvalTraits, SpaceCfg_>& data,
          const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& DOXY(trafo_data)) const
        {
          data.phi[0].value = DataType(1);
        }
      };
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_DISCONTINUOUS_EVALUATOR_HPP
