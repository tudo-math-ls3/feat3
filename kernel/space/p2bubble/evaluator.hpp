// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_P2BUBBLE_EVALUATOR_HPP
#define KERNEL_SPACE_P2BUBBLE_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/space/parametric_evaluator.hpp>
#include <kernel/space/dof_mapping_common.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace P2Bubble
    {
      /**
       * \brief P2-Bubble Element Evaluator reference capabilities
       *
       * \author Peter Zajac
       */
      static constexpr SpaceTags ref_caps = SpaceTags::ref_value | SpaceTags::ref_grad | SpaceTags::ref_hess;

      /**
       * \brief P2-Bubble Element Evaluator class template declaration.
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
       * \brief P2-Bubble Element evaluator implementation for Triangle shape
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Simplex<2> > :
        public ParametricEvaluator<
          Evaluator<
            Space_,
            TrafoEvaluator_,
            SpaceEvalTraits_,
            Shape::Simplex<2> >,
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

        /**
         * \brief Returns the number of local DOFs.
         *
         * \returns
         * The number of local dofs.
         */
        int get_num_local_dofs() const
        {
          return 7;
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
          // vertex dofs
          data.phi[0].ref_value = (DataType(1) + DataType(3)*x*y - DataType(2)*(x + y))*(DataType(1) - x - y);
          data.phi[1].ref_value = x*(DataType(2)*x + DataType(3)*y*(DataType(1) - x - y) - DataType(1));
          data.phi[2].ref_value = y*(DataType(2)*y + DataType(3)*x*(DataType(1) - x - y) - DataType(1));
          // edge dofs
          data.phi[3].ref_value = DataType(4)*x*y*(DataType(3)*(x+y)-DataType(2));
          data.phi[4].ref_value = DataType(4)*y*(DataType(3)*x-DataType(1))*(x+y-DataType(1));
          data.phi[5].ref_value = DataType(4)*x*(DataType(3)*y-DataType(1))*(x+y-DataType(1));
          // centre dof
          data.phi[6].ref_value = DataType(27)*x*y*(DataType(1)-x-y);
        }

        /**
         * \brief Evaluates the basis function gradients on the reference cell.
         *
         * \param[out] data
         * A reference to a basis gradient vector receiveing the result.
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
          // vertex dofs
          data.phi[0].ref_grad[0] = y*(DataType(7)-DataType(3)*y-DataType(6)*x)+DataType(4)*x-DataType(3);
          data.phi[0].ref_grad[1] = x*(DataType(7)-DataType(3)*x-DataType(6)*y)+DataType(4)*y-DataType(3);
          data.phi[1].ref_grad[0] = DataType(3)*y*(DataType(1)-y-DataType(2)*x)+DataType(4)*x-DataType(1);
          data.phi[1].ref_grad[1] = DataType(3)*x*(DataType(1)-x-DataType(2)*y);
          data.phi[2].ref_grad[0] = DataType(3)*y*(DataType(1)-y-DataType(2)*x);
          data.phi[2].ref_grad[1] = DataType(3)*x*(DataType(1)-x-DataType(2)*y)+DataType(4)*y-DataType(1);
          // edge dofs
          data.phi[3].ref_grad[0] = DataType(4)*y*(DataType(6)*x+DataType(3)*y-DataType(2));
          data.phi[3].ref_grad[1] = DataType(4)*x*(DataType(3)*x+DataType(6)*y-DataType(2));
          data.phi[4].ref_grad[0] = DataType(4)*y*(DataType(6)*x+DataType(3)*y-DataType(4));
          data.phi[4].ref_grad[1] = DataType(4)*(DataType(3)*x-DataType(1))*(x+DataType(2)*y-DataType(1));
          data.phi[5].ref_grad[0] = DataType(4)*(DataType(3)*y-DataType(1))*(y+DataType(2)*x-DataType(1));
          data.phi[5].ref_grad[1] = DataType(4)*x*(DataType(3)*x+DataType(6)*y-DataType(4));
          // centre dof
          data.phi[6].ref_grad[0] = DataType(27)*y*(DataType(1)-y-DataType(2)*x);
          data.phi[6].ref_grad[1] = DataType(27)*x*(DataType(1)-x-DataType(2)*y);
        }

        /**
         * \brief Evaluates the basis function hessians on the reference cell.
         *
         * \param[out] data
         * A reference to a basis hessian vector receiveing the result.
         *
         * \param[in] point
         * A reference to the point on the reference cell where to evaluate.
         */
        template<typename EvalData_>
        void eval_ref_hessians(
          EvalData_& data,
          const DomainPointType& point) const
        {
          const DataType x = point[0];
          const DataType y = point[1];
          // vertex dofs
          data.phi[0].ref_hess[0][0] = -DataType(6)*y+DataType(4);
          data.phi[0].ref_hess[1][1] =-DataType(6)*x+DataType(4);
          data.phi[0].ref_hess[0][1] = data.phi[0].ref_hess[1][0] = DataType(7)-DataType(6)*(x+y);
          data.phi[1].ref_hess[0][0] = -DataType(6)*y+DataType(4);
          data.phi[1].ref_hess[1][1] = -DataType(6)*x;
          data.phi[1].ref_hess[0][1] = data.phi[1].ref_hess[1][0] = DataType(3)-DataType(6)*(x+y);
          data.phi[2].ref_hess[0][0] = -DataType(6)*y;
          data.phi[2].ref_hess[1][1] = -DataType(6)*x+DataType(4);
          data.phi[2].ref_hess[0][1] = data.phi[2].ref_hess[1][0] = DataType(3)-DataType(6)*(x+y);
          // edge dofs
          data.phi[3].ref_hess[0][0] = DataType(24)*y;
          data.phi[3].ref_hess[1][1] = DataType(24)*x;
          data.phi[3].ref_hess[0][1] = data.phi[3].ref_hess[1][0] = DataType(8)*(DataType(3)*(x+y)-DataType(1));
          data.phi[4].ref_hess[0][0] = DataType(24)*y;
          data.phi[4].ref_hess[1][1] = DataType(24)*x-DataType(8);
          data.phi[4].ref_hess[0][1] = data.phi[4].ref_hess[1][0] = DataType(8)*(DataType(3)*(x+y)-DataType(2));
          data.phi[5].ref_hess[0][0] = DataType(24)*y-DataType(8);
          data.phi[5].ref_hess[1][1] = DataType(24)*x;
          data.phi[5].ref_hess[0][1] = data.phi[5].ref_hess[1][0] = DataType(8)*(DataType(3)*(x+y)-DataType(2));
          // centre dof
          data.phi[6].ref_hess[0][0] = -DataType(54)*y;
          data.phi[6].ref_hess[1][1] = -DataType(54)*x;
          data.phi[6].ref_hess[0][1] = data.phi[6].ref_hess[1][0] = DataType(27)*(DataType(1)-DataType(2)*(x+y));
        }
      }; // class Evaluator<...,Simplex<2>>
    } // namespace P2Bubble
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_P2BUBBLE_EVALUATOR_HPP
