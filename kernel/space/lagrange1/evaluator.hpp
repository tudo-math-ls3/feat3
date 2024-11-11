// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/space/parametric_evaluator.hpp>
#include <kernel/space/dof_mapping_common.hpp>

namespace FEAT
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
      static constexpr SpaceTags ref_caps = SpaceTags::ref_value | SpaceTags::ref_grad;

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
       * \brief Lagrange-1 Element evaluator implementation for Triangle shape
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
          return 3;
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
          data.phi[0].ref_value = DataType(1) - point[0] - point[1];
          data.phi[1].ref_value = point[0];
          data.phi[2].ref_value = point[1];
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
          const DomainPointType& DOXY(point)) const
        {
          data.phi[0].ref_grad[0] = -DataType(1);
          data.phi[0].ref_grad[1] = -DataType(1);
          data.phi[1].ref_grad[0] = DataType(1);
          data.phi[1].ref_grad[1] = DataType(0);
          data.phi[2].ref_grad[0] = DataType(0);
          data.phi[2].ref_grad[1] = DataType(1);
        }
      }; // class Evaluator<...,Simplex<2>>

      /**
       * \brief Lagrange-1 Element evaluator implementation for Tetrahedron shape
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Simplex<3> > :
        public ParametricEvaluator<
          Evaluator<
            Space_,
            TrafoEvaluator_,
            SpaceEvalTraits_,
            Shape::Simplex<3> >,
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
          return 4;
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
          data.phi[0].ref_value = DataType(1) - point[0] - point[1] - point[2];
          data.phi[1].ref_value = point[0];
          data.phi[2].ref_value = point[1];
          data.phi[3].ref_value = point[2];
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
          const DomainPointType& DOXY(point)) const
        {
          data.phi[0].ref_grad[0] = -DataType(1);
          data.phi[0].ref_grad[1] = -DataType(1);
          data.phi[0].ref_grad[2] = -DataType(1);
          data.phi[1].ref_grad[0] = DataType(1);
          data.phi[1].ref_grad[1] = DataType(0);
          data.phi[1].ref_grad[2] = DataType(0);
          data.phi[2].ref_grad[0] = DataType(0);
          data.phi[2].ref_grad[1] = DataType(1);
          data.phi[2].ref_grad[2] = DataType(0);
          data.phi[3].ref_grad[0] = DataType(0);
          data.phi[3].ref_grad[1] = DataType(0);
          data.phi[3].ref_grad[2] = DataType(1);
        }
      }; // class Evaluator<...,Simplex<3>>

      /**
       * \brief Lagrange-1 Element evaluator implementation for 1D Hypercube shape
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Hypercube<1> > :
        public ParametricEvaluator<
          Evaluator<
            Space_,
            TrafoEvaluator_,
            SpaceEvalTraits_,
            Shape::Hypercube<1> >,
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
          return 2;
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
          data.phi[0].ref_value = DataType(0.5) * (DataType(1) - point[0]);
          data.phi[1].ref_value = DataType(0.5) * (DataType(1) + point[0]);
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
          const DomainPointType& DOXY(point)) const
        {
          data.phi[0].ref_grad[0] = DataType(-0.5);
          data.phi[1].ref_grad[0] = DataType( 0.5);
        }
      }; // class Evaluator<...,Hypercube<1>>

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
        public ParametricEvaluator<
          Evaluator<
            Space_,
            TrafoEvaluator_,
            SpaceEvalTraits_,
            Shape::Hypercube<2> >,
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
          return 4;
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
          data.phi[0].ref_value = DataType(0.25) * (DataType(1) - point[0]) * (DataType(1) - point[1]);
          data.phi[1].ref_value = DataType(0.25) * (DataType(1) + point[0]) * (DataType(1) - point[1]);
          data.phi[2].ref_value = DataType(0.25) * (DataType(1) - point[0]) * (DataType(1) + point[1]);
          data.phi[3].ref_value = DataType(0.25) * (DataType(1) + point[0]) * (DataType(1) + point[1]);
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
          data.phi[0].ref_grad[0] = DataType(-0.25) * (DataType(1) - point[1]);
          data.phi[0].ref_grad[1] = DataType(-0.25) * (DataType(1) - point[0]);
          data.phi[1].ref_grad[0] = DataType( 0.25) * (DataType(1) - point[1]);
          data.phi[1].ref_grad[1] = DataType(-0.25) * (DataType(1) + point[0]);
          data.phi[2].ref_grad[0] = DataType(-0.25) * (DataType(1) + point[1]);
          data.phi[2].ref_grad[1] = DataType( 0.25) * (DataType(1) - point[0]);
          data.phi[3].ref_grad[0] = DataType( 0.25) * (DataType(1) + point[1]);
          data.phi[3].ref_grad[1] = DataType( 0.25) * (DataType(1) + point[0]);
        }
      }; // class Evaluator<...,Hypercube<2>>

      /**
       * \brief Lagrange-1 Element evaluator implementation for Hexahedron shape
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Hypercube<3> > :
        public ParametricEvaluator<
          Evaluator<
            Space_,
            TrafoEvaluator_,
            SpaceEvalTraits_,
            Shape::Hypercube<3> >,
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
          return 8;
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
          data.phi[0].ref_value = DataType(0.125) * (DataType(1) - point[0]) * (DataType(1) - point[1]) * (DataType(1) - point[2]);
          data.phi[1].ref_value = DataType(0.125) * (DataType(1) + point[0]) * (DataType(1) - point[1]) * (DataType(1) - point[2]);
          data.phi[2].ref_value = DataType(0.125) * (DataType(1) - point[0]) * (DataType(1) + point[1]) * (DataType(1) - point[2]);
          data.phi[3].ref_value = DataType(0.125) * (DataType(1) + point[0]) * (DataType(1) + point[1]) * (DataType(1) - point[2]);
          data.phi[4].ref_value = DataType(0.125) * (DataType(1) - point[0]) * (DataType(1) - point[1]) * (DataType(1) + point[2]);
          data.phi[5].ref_value = DataType(0.125) * (DataType(1) + point[0]) * (DataType(1) - point[1]) * (DataType(1) + point[2]);
          data.phi[6].ref_value = DataType(0.125) * (DataType(1) - point[0]) * (DataType(1) + point[1]) * (DataType(1) + point[2]);
          data.phi[7].ref_value = DataType(0.125) * (DataType(1) + point[0]) * (DataType(1) + point[1]) * (DataType(1) + point[2]);
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
          data.phi[0].ref_grad[0] = DataType(-0.125) * (DataType(1) - point[1]) * (DataType(1) - point[2]);
          data.phi[0].ref_grad[1] = DataType(-0.125) * (DataType(1) - point[2]) * (DataType(1) - point[0]);
          data.phi[0].ref_grad[2] = DataType(-0.125) * (DataType(1) - point[0]) * (DataType(1) - point[1]);

          data.phi[1].ref_grad[0] = DataType( 0.125) * (DataType(1) - point[1]) * (DataType(1) - point[2]);
          data.phi[1].ref_grad[1] = DataType(-0.125) * (DataType(1) - point[2]) * (DataType(1) + point[0]);
          data.phi[1].ref_grad[2] = DataType(-0.125) * (DataType(1) + point[0]) * (DataType(1) - point[1]);

          data.phi[2].ref_grad[0] = DataType(-0.125) * (DataType(1) + point[1]) * (DataType(1) - point[2]);
          data.phi[2].ref_grad[1] = DataType( 0.125) * (DataType(1) - point[2]) * (DataType(1) - point[0]);
          data.phi[2].ref_grad[2] = DataType(-0.125) * (DataType(1) - point[0]) * (DataType(1) + point[1]);

          data.phi[3].ref_grad[0] = DataType( 0.125) * (DataType(1) + point[1]) * (DataType(1) - point[2]);
          data.phi[3].ref_grad[1] = DataType( 0.125) * (DataType(1) - point[2]) * (DataType(1) + point[0]);
          data.phi[3].ref_grad[2] = DataType(-0.125) * (DataType(1) + point[0]) * (DataType(1) + point[1]);

          data.phi[4].ref_grad[0] = DataType(-0.125) * (DataType(1) - point[1]) * (DataType(1) + point[2]);
          data.phi[4].ref_grad[1] = DataType(-0.125) * (DataType(1) + point[2]) * (DataType(1) - point[0]);
          data.phi[4].ref_grad[2] = DataType( 0.125) * (DataType(1) - point[0]) * (DataType(1) - point[1]);

          data.phi[5].ref_grad[0] = DataType( 0.125) * (DataType(1) - point[1]) * (DataType(1) + point[2]);
          data.phi[5].ref_grad[1] = DataType(-0.125) * (DataType(1) + point[2]) * (DataType(1) + point[0]);
          data.phi[5].ref_grad[2] = DataType( 0.125) * (DataType(1) + point[0]) * (DataType(1) - point[1]);

          data.phi[6].ref_grad[0] = DataType(-0.125) * (DataType(1) + point[1]) * (DataType(1) + point[2]);
          data.phi[6].ref_grad[1] = DataType( 0.125) * (DataType(1) + point[2]) * (DataType(1) - point[0]);
          data.phi[6].ref_grad[2] = DataType( 0.125) * (DataType(1) - point[0]) * (DataType(1) + point[1]);

          data.phi[7].ref_grad[0] = DataType( 0.125) * (DataType(1) + point[1]) * (DataType(1) + point[2]);
          data.phi[7].ref_grad[1] = DataType( 0.125) * (DataType(1) + point[2]) * (DataType(1) + point[0]);
          data.phi[7].ref_grad[2] = DataType( 0.125) * (DataType(1) + point[0]) * (DataType(1) + point[1]);
        }
      }; // class Evaluator<...,Hypercube<3>>
    } // namespace Lagrange1
  } // namespace Space
} // namespace FEAT
