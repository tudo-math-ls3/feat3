// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_CAI_DOU_SAN_SHE_YE_EVALUATOR_HPP
#define KERNEL_SPACE_CAI_DOU_SAN_SHE_YE_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/space/parametric_evaluator.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace CaiDouSanSheYe
    {
      /**
       * \brief Cai-Douglas-Santos-Sheen-Ye Element Evaluator reference capabilities
       *
       * \author Peter Zajac
       */
      static constexpr SpaceTags ref_caps = SpaceTags::ref_value | SpaceTags::ref_grad;

      /**
       * \brief Cai-Douglas-Santos-Sheen-Ye Element Evaluator class template declaration.
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
       * \brief Implementation of standard parametric quadrilateral Cai-Douglas-Santos-Sheen-Ye evaluator
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Quadrilateral> :
        public ParametricEvaluator<
          Evaluator< Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Quadrilateral>,
          TrafoEvaluator_,
          SpaceEvalTraits_,
          ref_caps>
      {
      public:
        /// base-class typedef
        typedef ParametricEvaluator<Evaluator, TrafoEvaluator_, SpaceEvalTraits_, ref_caps> BaseClass;

        /// space type
        typedef Space_ SpaceType;

        /// space evaluation traits
        typedef SpaceEvalTraits_ SpaceEvalTraits;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const SpaceType& DOXY(space))
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
        int get_num_local_dofs() const
        {
          return 5;
        }

        /**
        * \brief Evaluates the basis function values on the reference cell.
        *
        * \param[out] data
        * A reference to a basis value vector receiving the result.
        *
        * \param[in] point
        * A reference to the point on the reference cell where to evaluate.
        */
        template<typename EvalData_>
        void eval_ref_values(
          EvalData_& data,
          const DomainPointType& point) const
        {
          const DataType x = point[0];
          const DataType y = point[1];
          const DataType x2 = x*x;
          const DataType y2 = y*y;
          // edge midpoints dofs
          data.phi[0].ref_value = DataType(0.25) - DataType(0.5)*y + ( DataType(0.375) - DataType(0.625)*x2)*x2 + (-DataType(0.375) + DataType(0.625)*y2)*y2;
          data.phi[1].ref_value = DataType(0.25) + DataType(0.5)*y + ( DataType(0.375) - DataType(0.625)*x2)*x2 + (-DataType(0.375) + DataType(0.625)*y2)*y2;
          data.phi[2].ref_value = DataType(0.25) - DataType(0.5)*x + (-DataType(0.375) + DataType(0.625)*x2)*x2 + ( DataType(0.375) - DataType(0.625)*y2)*y2;
          data.phi[3].ref_value = DataType(0.25) + DataType(0.5)*x + (-DataType(0.375) + DataType(0.625)*x2)*x2 + ( DataType(0.375) - DataType(0.625)*y2)*y2;
          // bubble dof
          data.phi[4].ref_value = DataType(2.25)*x*y;
        }

        /**
        * \brief Evaluates the basis function gradients on the reference cell.
        *
        * \param[out] data
        * A reference to a basis gradient vector receiving the result.
        *
        * \param[in] point
        * A reference to the point on the reference cell where to evaluate.
        */
        template<typename EvalData_>
        void eval_ref_gradients(
          EvalData_& data,
          const DomainPointType& point) const
        {
          const DataType x = point[0];
          const DataType y = point[1];
          const DataType x2 = x*x;
          const DataType y2 = y*y;
          // edge midpoint dofs
          data.phi[0].ref_grad[0] = ( DataType(0.75) - DataType(2.5)*x2)*x;
          data.phi[0].ref_grad[1] = (-DataType(0.75) + DataType(2.5)*y2)*y - DataType(0.5);
          data.phi[1].ref_grad[0] = ( DataType(0.75) - DataType(2.5)*x2)*x;
          data.phi[1].ref_grad[1] = (-DataType(0.75) + DataType(2.5)*y2)*y + DataType(0.5);
          data.phi[2].ref_grad[0] = (-DataType(0.75) + DataType(2.5)*x2)*x - DataType(0.5);
          data.phi[2].ref_grad[1] = ( DataType(0.75) - DataType(2.5)*y2)*y;
          data.phi[3].ref_grad[0] = (-DataType(0.75) + DataType(2.5)*x2)*x + DataType(0.5);
          data.phi[3].ref_grad[1] = ( DataType(0.75) - DataType(2.5)*y2)*y;
          // bubble dof
          data.phi[4].ref_grad[0] = DataType(2.25)*y;
          data.phi[4].ref_grad[1] = DataType(2.25)*x;
        }
      }; // Evaluator<..., Shape::Quadrilateral>
    } // namespace CaiDouSanSheYe
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_CAI_DOU_SAN_SHE_YE_EVALUATOR_HPP
