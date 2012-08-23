#pragma once
#ifndef KERNEL_SPACE_LAGRANGE1_EVALUATOR_HPP
#define KERNEL_SPACE_LAGRANGE1_EVALUATOR_HPP 1

// includes, FEAST
#include <kernel/space/evaluator_base.hpp>
#include <kernel/space/dof_mapping_common.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Lagrange1
    {
      /**
       * \brief Lagrange-1 Element Evaluator reference capabilities
       *
       * \author Peter Zajac
       */
      struct ReferenceCapabilities
      {
        /// dummy enum
        enum
        {
          /// can compute reference function values
          can_ref_value = 1,
          /// can compute reference gradients
          can_ref_grad = 1
        };
      };

      /**
       * \brief Lagrange-1 Element Evaluator class template declaration.
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_,
        typename Shape_ = typename Space_::ShapeType>
      class Evaluator DOXY({});

      /**
       * \brief Lagrange-1 Element evaluator implementation for Quadrilateral shape
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Hypercube<2> > :
        public EvaluatorParametric<
          Evaluator<
            Space_,
            TrafoEvaluator_,
            SpaceEvalTraits_,
            Shape::Hypercube<2> >,
          TrafoEvaluator_,
          SpaceEvalTraits_,
          ReferenceCapabilities>
      {
      public:
        /// base-class typedef
        typedef EvaluatorParametric<Evaluator, TrafoEvaluator_, SpaceEvalTraits_, ReferenceCapabilities> BaseClass;

        /// space type
        typedef Space_ SpaceType;

        /// space evaluation traits
        typedef SpaceEvalTraits_ SpaceEvalTraits;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// const domain point reference
        typedef typename EvalPolicy::DomainPointConstRef DomainPointConstRef;

        /// basis value coefficient type
        typedef typename SpaceEvalTraits::BasisValueCoeff BasisValueCoeff;
        /// basis value vector reference
        typedef typename SpaceEvalTraits::BasisValueVectorRef BasisValueVectorRef;

        /// basis gradient coefficient type
        typedef typename SpaceEvalTraits::BasisGradientCoeff BasisGradientCoeff;
        /// basis gradent vector reference
        typedef typename SpaceEvalTraits::BasisGradientVectorRef BasisGradientVectorRef;

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const SpaceType& /*space*/)
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
          return 4;
        }

        /**
         * \brief Evaluates the basis function values on the reference cell.
         *
         * \param[out] values
         * A reference to a basis value vector receiving the result.
         *
         * \param[in] point
         * A reference to the point on the reference cell where to evaluate.
         */
        void eval_ref_values(
          BasisValueVectorRef values,
          DomainPointConstRef point) const
        {
          values[0] = BasisValueCoeff(0.25 * (1.0 - point[0]) * (1.0 - point[1]));
          values[1] = BasisValueCoeff(0.25 * (1.0 + point[0]) * (1.0 - point[1]));
          values[2] = BasisValueCoeff(0.25 * (1.0 - point[0]) * (1.0 + point[1]));
          values[3] = BasisValueCoeff(0.25 * (1.0 + point[0]) * (1.0 + point[1]));
        }

        /**
         * \brief Evaluates the basis function gradients on the reference cell.
         *
         * \param[out] grads
         * A reference to a basis gradient vector receiveing the result.
         *
         * \param[in] point
         * A reference to the point on the reference cell where to evaluate.
         */
        void eval_ref_gradients(
          BasisGradientVectorRef grads,
          DomainPointConstRef point) const
        {
          grads[0][0] = BasisGradientCoeff(-0.25 * (1.0 - point[1]));
          grads[0][1] = BasisGradientCoeff(-0.25 * (1.0 - point[0]));
          grads[1][0] = BasisGradientCoeff( 0.25 * (1.0 - point[1]));
          grads[1][1] = BasisGradientCoeff(-0.25 * (1.0 + point[0]));
          grads[2][0] = BasisGradientCoeff(-0.25 * (1.0 + point[1]));
          grads[2][1] = BasisGradientCoeff( 0.25 * (1.0 - point[0]));
          grads[3][0] = BasisGradientCoeff( 0.25 * (1.0 + point[1]));
          grads[3][1] = BasisGradientCoeff( 0.25 * (1.0 + point[0]));
        }
      }; // class Evaluator<...,Hypercube<2>>
    } // namespace Lagrange1
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_LAGRANGE1_EVALUATOR_HPP
