#pragma once
#ifndef KERNEL_SPACE_LAGRANGE2_EVALUATOR_HPP
#define KERNEL_SPACE_LAGRANGE2_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/space/parametric_evaluator.hpp>
#include <kernel/space/dof_mapping_common.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Lagrange2
    {
      /**
       * \brief Lagrange-2 Element Evaluator reference capabilities
       *
       * \author Peter Zajac
       */
      static constexpr SpaceTags ref_caps = SpaceTags::ref_value | SpaceTags::ref_grad | SpaceTags::ref_hess;

      /// \cond internal
      namespace Intern
      {
        // p0, p1 and p2 are the 1D basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the Hypercube evaluators.

        // basis function for left vertex
        template<typename T_>
        inline T_ p0(T_ x)
        {
          return T_(0.5) * x * (x - T_(1));
        }

        // basis function for right vertex
        template<typename T_>
        inline T_ p1(T_ x)
        {
          return T_(0.5) * x * (x + T_(1));
        }

        // basis function for edge midpoint
        template<typename T_>
        inline T_ p2(T_ x)
        {
          return (T_(1) - x) * (T_(1) + x);
        }

        // first order derivatives

        template<typename T_>
        inline T_ d1p0(T_ x)
        {
          return x - T_(0.5);
        }

        template<typename T_>
        inline T_ d1p1(T_ x)
        {
          return x + T_(0.5);
        }

        template<typename T_>
        inline T_ d1p2(T_ x)
        {
          return -T_(2) * x;
        }

        // second order derivatives

        template<typename T_>
        inline T_ d2p0(T_)
        {
          return T_(1);
        }

        template<typename T_>
        inline T_ d2p1(T_)
        {
          return T_(1);
        }

        template<typename T_>
        inline T_ d2p2(T_)
        {
          return -T_(2);
        }
      } // namespace Intern
      /// \endcond

      /**
       * \brief Lagrange-2 Element Evaluator class template declaration.
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
       * \brief Lagrange-2 Element evaluator implementation for Triangle shape
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
          return 6;
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
          // vertex dofs
          data.phi[0].ref_value = DataType(2) * (point[0] + point[1] - DataType(0.5)) * (point[0] + point[1] - DataType(1));
          data.phi[1].ref_value = DataType(2) * point[0] * (point[0] - DataType(0.5));
          data.phi[2].ref_value = DataType(2) * point[1] * (point[1] - DataType(0.5));
          // edge dofs
          data.phi[3].ref_value = DataType(4) * point[0] * point[1];
          data.phi[4].ref_value = DataType(4) * point[1] * (DataType(1) - point[0] - point[1]);
          data.phi[5].ref_value = DataType(4) * point[0] * (DataType(1) - point[0] - point[1]);
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
          // vertex dofs
          data.phi[0].ref_grad[0] = DataType(4) * (point[0] + point[1]) - DataType(3);
          data.phi[0].ref_grad[1] = DataType(4) * (point[0] + point[1]) - DataType(3);
          data.phi[1].ref_grad[0] = DataType(4) * point[0] - DataType(1);
          data.phi[1].ref_grad[1] = DataType(0);
          data.phi[2].ref_grad[0] = DataType(0);
          data.phi[2].ref_grad[1] = DataType(4) * point[1] - DataType(1);
          // edge dofs
          data.phi[3].ref_grad[0] = DataType(4) * point[1];
          data.phi[3].ref_grad[1] = DataType(4) * point[0];
          data.phi[4].ref_grad[0] = -DataType(4) * point[1];
          data.phi[4].ref_grad[1] = -DataType(4) * (DataType(2)*point[1] + point[0] - DataType(1));
          data.phi[5].ref_grad[0] = -DataType(4) * (DataType(2)*point[0] + point[1] - DataType(1));
          data.phi[5].ref_grad[1] = -DataType(4) * point[0];
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
          const DomainPointType& DOXY(point)) const
        {
          // vertex dofs
          data.phi[0].ref_hess[0][0] = data.phi[0].ref_hess[1][1] =
          data.phi[0].ref_hess[0][1] = data.phi[0].ref_hess[1][0] = DataType(4);
          data.phi[1].ref_hess[0][0] = DataType(4);
          data.phi[1].ref_hess[1][1] = data.phi[1].ref_hess[0][1] = data.phi[1].ref_hess[1][0] = DataType(0);
          data.phi[2].ref_hess[1][1] = DataType(4);
          data.phi[2].ref_hess[0][0] = data.phi[2].ref_hess[0][1] = data.phi[2].ref_hess[1][0] = DataType(0);
          // edge dofs
          data.phi[3].ref_hess[0][0] = data.phi[3].ref_hess[1][1] = DataType(0);
          data.phi[3].ref_hess[0][1] = data.phi[3].ref_hess[1][0] = DataType(4);
          data.phi[4].ref_hess[0][0] = DataType(0);
          data.phi[4].ref_hess[1][1] = -DataType(8);
          data.phi[4].ref_hess[0][1] = data.phi[4].ref_hess[1][0] = -DataType(4);
          data.phi[5].ref_hess[0][0] = -DataType(8);
          data.phi[5].ref_hess[1][1] = DataType(0);
          data.phi[5].ref_hess[0][1] = data.phi[5].ref_hess[1][0] = -DataType(4);
        }
      }; // class Evaluator<...,Simplex<2>>

      /**
       * \brief Lagrange-2 Element evaluator implementation for Tetrahedron shape
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
          return 10;
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
          // vertex dofs
          data.phi[0].ref_value = DataType(2) * (point[0] + point[1] + point[2] - DataType(0.5)) * (point[0] + point[1] + point[2] - DataType(1));
          data.phi[1].ref_value = DataType(2) * point[0] * (point[0] - DataType(0.5));
          data.phi[2].ref_value = DataType(2) * point[1] * (point[1] - DataType(0.5));
          data.phi[3].ref_value = DataType(2) * point[2] * (point[2] - DataType(0.5));
          // edge dofs
          data.phi[4].ref_value = DataType(4) * point[0] * (DataType(1) - point[0] - point[1] - point[2]);
          data.phi[5].ref_value = DataType(4) * point[1] * (DataType(1) - point[0] - point[1] - point[2]);
          data.phi[6].ref_value = DataType(4) * point[2] * (DataType(1) - point[0] - point[1] - point[2]);
          data.phi[7].ref_value = DataType(4) * point[0] * point[1];
          data.phi[8].ref_value = DataType(4) * point[0] * point[2];
          data.phi[9].ref_value = DataType(4) * point[1] * point[2];
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
          // vertex dofs
          data.phi[0].ref_grad[0] = DataType(4) * (point[0] + point[1] + point[2]) - DataType(3);
          data.phi[0].ref_grad[1] = DataType(4) * (point[0] + point[1] + point[2]) - DataType(3);
          data.phi[0].ref_grad[2] = DataType(4) * (point[0] + point[1] + point[2]) - DataType(3);
          data.phi[1].ref_grad[0] = DataType(4) * point[0] - DataType(1);
          data.phi[1].ref_grad[1] = DataType(0);
          data.phi[1].ref_grad[2] = DataType(0);
          data.phi[2].ref_grad[0] = DataType(0);
          data.phi[2].ref_grad[1] = DataType(4) * point[1] - DataType(1);
          data.phi[2].ref_grad[2] = DataType(0);
          data.phi[3].ref_grad[0] = DataType(0);
          data.phi[3].ref_grad[1] = DataType(0);
          data.phi[3].ref_grad[2] = DataType(4) * point[2] - DataType(1);
          // edge dofs
          data.phi[4].ref_grad[0] = -DataType(4) * (DataType(2)*point[0] + point[1] + point[2] - DataType(1));
          data.phi[4].ref_grad[1] = -DataType(4) * point[0];
          data.phi[4].ref_grad[2] = -DataType(4) * point[0];
          data.phi[5].ref_grad[0] = -DataType(4) * point[1];
          data.phi[5].ref_grad[1] = -DataType(4) * (DataType(2)*point[1] + point[0] + point[2] - DataType(1));
          data.phi[5].ref_grad[2] = -DataType(4) * point[1];
          data.phi[6].ref_grad[0] = -DataType(4) * point[2];
          data.phi[6].ref_grad[1] = -DataType(4) * point[2];
          data.phi[6].ref_grad[2] = -DataType(4) * (DataType(2)*point[2] + point[0] + point[1] - DataType(1));
          data.phi[7].ref_grad[0] = DataType(4) * point[1];
          data.phi[7].ref_grad[1] = DataType(4) * point[0];
          data.phi[7].ref_grad[2] = DataType(0);
          data.phi[8].ref_grad[0] = DataType(4) * point[2];
          data.phi[8].ref_grad[1] = DataType(0);
          data.phi[8].ref_grad[2] = DataType(4) * point[0];
          data.phi[9].ref_grad[0] = DataType(0);
          data.phi[9].ref_grad[1] = DataType(4) * point[2];
          data.phi[9].ref_grad[2] = DataType(4) * point[1];
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
          const DomainPointType& DOXY(point)) const
        {
          // vertex dofs
          data.phi[0].ref_hess[0][0] = DataType(4);
          data.phi[0].ref_hess[1][1] = DataType(4);
          data.phi[0].ref_hess[2][2] = DataType(4);
          data.phi[0].ref_hess[0][1] = data.phi[0].ref_hess[1][0] = DataType(4);
          data.phi[0].ref_hess[0][2] = data.phi[0].ref_hess[2][0] = DataType(4);
          data.phi[0].ref_hess[1][2] = data.phi[0].ref_hess[2][1] = DataType(4);

          data.phi[1].ref_hess[0][0] = DataType(4);
          data.phi[1].ref_hess[1][1] = DataType(0);
          data.phi[1].ref_hess[2][2] = DataType(0);
          data.phi[1].ref_hess[0][1] = data.phi[1].ref_hess[1][0] = DataType(0);
          data.phi[1].ref_hess[0][2] = data.phi[1].ref_hess[2][0] = DataType(0);
          data.phi[1].ref_hess[1][2] = data.phi[1].ref_hess[2][1] = DataType(0);

          data.phi[2].ref_hess[0][0] = DataType(0);
          data.phi[2].ref_hess[1][1] = DataType(4);
          data.phi[2].ref_hess[2][2] = DataType(0);
          data.phi[2].ref_hess[0][1] = data.phi[2].ref_hess[1][0] = DataType(0);
          data.phi[2].ref_hess[0][2] = data.phi[2].ref_hess[2][0] = DataType(0);
          data.phi[2].ref_hess[1][2] = data.phi[2].ref_hess[2][1] = DataType(0);

          data.phi[3].ref_hess[0][0] = DataType(0);
          data.phi[3].ref_hess[1][1] = DataType(0);
          data.phi[3].ref_hess[2][2] = DataType(4);
          data.phi[3].ref_hess[0][1] = data.phi[3].ref_hess[1][0] = DataType(0);
          data.phi[3].ref_hess[0][2] = data.phi[3].ref_hess[2][0] = DataType(0);
          data.phi[3].ref_hess[1][2] = data.phi[3].ref_hess[2][1] = DataType(0);

          // edge dofs
          data.phi[4].ref_hess[0][0] = -DataType(8);
          data.phi[4].ref_hess[1][1] = DataType(0);
          data.phi[4].ref_hess[2][2] = DataType(0);
          data.phi[4].ref_hess[0][1] = data.phi[4].ref_hess[1][0] = -DataType(4);
          data.phi[4].ref_hess[0][2] = data.phi[4].ref_hess[2][0] = -DataType(4);
          data.phi[4].ref_hess[1][2] = data.phi[4].ref_hess[2][1] = DataType(0);

          data.phi[5].ref_hess[0][0] = DataType(0);
          data.phi[5].ref_hess[1][1] = -DataType(8);
          data.phi[5].ref_hess[2][2] = DataType(0);
          data.phi[5].ref_hess[0][1] = data.phi[5].ref_hess[1][0] = -DataType(4);
          data.phi[5].ref_hess[0][2] = data.phi[5].ref_hess[2][0] = DataType(0);
          data.phi[5].ref_hess[1][2] = data.phi[5].ref_hess[2][1] = -DataType(4);

          data.phi[6].ref_hess[0][0] = DataType(0);
          data.phi[6].ref_hess[1][1] = DataType(0);
          data.phi[6].ref_hess[2][2] = -DataType(8);
          data.phi[6].ref_hess[0][1] = data.phi[6].ref_hess[1][0] = DataType(0);
          data.phi[6].ref_hess[0][2] = data.phi[6].ref_hess[2][0] = -DataType(4);
          data.phi[6].ref_hess[1][2] = data.phi[6].ref_hess[2][1] = -DataType(4);

          data.phi[7].ref_hess[0][0] = DataType(0);
          data.phi[7].ref_hess[1][1] = DataType(0);
          data.phi[7].ref_hess[2][2] = DataType(0);
          data.phi[7].ref_hess[0][1] = data.phi[7].ref_hess[1][0] = DataType(4);
          data.phi[7].ref_hess[0][2] = data.phi[7].ref_hess[2][0] = DataType(0);
          data.phi[7].ref_hess[1][2] = data.phi[7].ref_hess[2][1] = DataType(0);

          data.phi[8].ref_hess[0][0] = DataType(0);
          data.phi[8].ref_hess[1][1] = DataType(0);
          data.phi[8].ref_hess[2][2] = DataType(0);
          data.phi[8].ref_hess[0][1] = data.phi[8].ref_hess[1][0] = DataType(0);
          data.phi[8].ref_hess[0][2] = data.phi[8].ref_hess[2][0] = DataType(4);
          data.phi[8].ref_hess[1][2] = data.phi[8].ref_hess[2][1] = DataType(0);

          data.phi[9].ref_hess[0][0] = DataType(0);
          data.phi[9].ref_hess[1][1] = DataType(0);
          data.phi[9].ref_hess[2][2] = DataType(0);
          data.phi[9].ref_hess[0][1] = data.phi[9].ref_hess[1][0] = DataType(0);
          data.phi[9].ref_hess[0][2] = data.phi[9].ref_hess[2][0] = DataType(0);
          data.phi[9].ref_hess[1][2] = data.phi[9].ref_hess[2][1] = DataType(4);

        }
      }; // class Evaluator<...,Simplex<3>>

      /**
       * \brief Lagrange-2 Element evaluator implementation for 1D Hypercube shape
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
          data.phi[0].ref_value = Intern::p0(point[0]);
          data.phi[1].ref_value = Intern::p1(point[0]);
          data.phi[2].ref_value = Intern::p2(point[0]);
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
          data.phi[0].ref_grad[0] = Intern::d1p0(point[0]);
          data.phi[1].ref_grad[0] = Intern::d1p1(point[0]);
          data.phi[2].ref_grad[0] = Intern::d1p2(point[0]);
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
          data.phi[0].ref_hess[0][0] = Intern::d2p0(point[0]);
          data.phi[1].ref_hess[0][0] = Intern::d2p1(point[0]);
          data.phi[2].ref_hess[0][0] = Intern::d2p2(point[0]);
        }
      }; // class Evaluator<...,Hypercube<1>>

      /**
       * \brief Lagrange-2 Element evaluator implementation for Quadrilateral shape
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
          return 9;
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
          using namespace Lagrange2::Intern;

          // vertex dofs
          data.phi[0].ref_value = p0(point[0]) * p0(point[1]);
          data.phi[1].ref_value = p1(point[0]) * p0(point[1]);
          data.phi[2].ref_value = p0(point[0]) * p1(point[1]);
          data.phi[3].ref_value = p1(point[0]) * p1(point[1]);
          // edge dofs
          data.phi[4].ref_value = p2(point[0]) * p0(point[1]);
          data.phi[5].ref_value = p2(point[0]) * p1(point[1]);
          data.phi[6].ref_value = p0(point[0]) * p2(point[1]);
          data.phi[7].ref_value = p1(point[0]) * p2(point[1]);
          // center dof
          data.phi[8].ref_value = p2(point[0]) * p2(point[1]);
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
          using namespace Lagrange2::Intern;

          // vertex dofs
          data.phi[0].ref_grad[0] = d1p0(point[0]) * p0(point[1]);
          data.phi[0].ref_grad[1] = p0(point[0]) * d1p0(point[1]);
          data.phi[1].ref_grad[0] = d1p1(point[0]) * p0(point[1]);
          data.phi[1].ref_grad[1] = p1(point[0]) * d1p0(point[1]);
          data.phi[2].ref_grad[0] = d1p0(point[0]) * p1(point[1]);
          data.phi[2].ref_grad[1] = p0(point[0]) * d1p1(point[1]);
          data.phi[3].ref_grad[0] = d1p1(point[0]) * p1(point[1]);
          data.phi[3].ref_grad[1] = p1(point[0]) * d1p1(point[1]);
          // edge dofs
          data.phi[4].ref_grad[0] = d1p2(point[0]) * p0(point[1]);
          data.phi[4].ref_grad[1] = p2(point[0]) * d1p0(point[1]);
          data.phi[5].ref_grad[0] = d1p2(point[0]) * p1(point[1]);
          data.phi[5].ref_grad[1] = p2(point[0]) * d1p1(point[1]);
          data.phi[6].ref_grad[0] = d1p0(point[0]) * p2(point[1]);
          data.phi[6].ref_grad[1] = p0(point[0]) * d1p2(point[1]);
          data.phi[7].ref_grad[0] = d1p1(point[0]) * p2(point[1]);
          data.phi[7].ref_grad[1] = p1(point[0]) * d1p2(point[1]);
          // center dof
          data.phi[8].ref_grad[0] = d1p2(point[0]) * p2(point[1]);
          data.phi[8].ref_grad[1] = p2(point[0]) * d1p2(point[1]);
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
          using namespace Lagrange2::Intern;

          // vertex dofs
          data.phi[0].ref_hess[0][0] = d2p0(point[0]) * p0(point[1]);
          data.phi[0].ref_hess[1][1] = p0(point[0]) * d2p0(point[1]);
          data.phi[0].ref_hess[1][0] =
          data.phi[0].ref_hess[0][1] = d1p0(point[0]) * d1p0(point[1]);
          data.phi[1].ref_hess[0][0] = d2p1(point[0]) * p0(point[1]);
          data.phi[1].ref_hess[1][1] = p1(point[0]) * d2p0(point[1]);
          data.phi[1].ref_hess[1][0] =
          data.phi[1].ref_hess[0][1] = d1p1(point[0]) * d1p0(point[1]);
          data.phi[2].ref_hess[0][0] = d2p0(point[0]) * p1(point[1]);
          data.phi[2].ref_hess[1][1] = p0(point[0]) * d2p1(point[1]);
          data.phi[2].ref_hess[1][0] =
          data.phi[2].ref_hess[0][1] = d1p0(point[0]) * d1p1(point[1]);
          data.phi[3].ref_hess[0][0] = d2p1(point[0]) * p1(point[1]);
          data.phi[3].ref_hess[1][1] = p1(point[0]) * d2p1(point[1]);
          data.phi[3].ref_hess[1][0] =
          data.phi[3].ref_hess[0][1] = d1p1(point[0]) * d1p1(point[1]);
          // edge dofs
          data.phi[4].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]);
          data.phi[4].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]);
          data.phi[4].ref_hess[1][0] =
          data.phi[4].ref_hess[0][1] = d1p2(point[0]) * d1p0(point[1]);
          data.phi[5].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]);
          data.phi[5].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]);
          data.phi[5].ref_hess[1][0] =
          data.phi[5].ref_hess[0][1] = d1p2(point[0]) * d1p1(point[1]);
          data.phi[6].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]);
          data.phi[6].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]);
          data.phi[6].ref_hess[1][0] =
          data.phi[6].ref_hess[0][1] = d1p0(point[0]) * d1p2(point[1]);
          data.phi[7].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]);
          data.phi[7].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]);
          data.phi[7].ref_hess[1][0] =
          data.phi[7].ref_hess[0][1] = d1p1(point[0]) * d1p2(point[1]);
          // center dof
          data.phi[8].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]);
          data.phi[8].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]);
          data.phi[8].ref_hess[1][0] =
          data.phi[8].ref_hess[0][1] = d1p2(point[0]) * d1p2(point[1]);
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
          return 27;
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
          using namespace Lagrange2::Intern;

          // vertex dofs
          data.phi[0].ref_value = p0(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[1].ref_value = p1(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[2].ref_value = p0(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[3].ref_value = p1(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[4].ref_value = p0(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[5].ref_value = p1(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[6].ref_value = p0(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[7].ref_value = p1(point[0]) * p1(point[1]) * p1(point[2]);
          // edge dofs
          data.phi[ 8].ref_value = p2(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 9].ref_value = p2(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[10].ref_value = p2(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[11].ref_value = p2(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[12].ref_value = p0(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[13].ref_value = p1(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[14].ref_value = p0(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[15].ref_value = p1(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[16].ref_value = p0(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[17].ref_value = p1(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[18].ref_value = p0(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[19].ref_value = p1(point[0]) * p1(point[1]) * p2(point[2]);
          // face dofs
          data.phi[20].ref_value = p2(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[21].ref_value = p2(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[22].ref_value = p2(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[23].ref_value = p2(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[24].ref_value = p0(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[25].ref_value = p1(point[0]) * p2(point[1]) * p2(point[2]);
          // center dof
          data.phi[26].ref_value = p2(point[0]) * p2(point[1]) * p2(point[2]);
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
          using namespace Lagrange2::Intern;
          // vertex dofs
          data.phi[0].ref_grad[0] = d1p0(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[0].ref_grad[1] = p0(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[0].ref_grad[2] = p0(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[1].ref_grad[0] = d1p1(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[1].ref_grad[1] = p1(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[1].ref_grad[2] = p1(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[2].ref_grad[0] = d1p0(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[2].ref_grad[1] = p0(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[2].ref_grad[2] = p0(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[3].ref_grad[0] = d1p1(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[3].ref_grad[1] = p1(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[3].ref_grad[2] = p1(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[4].ref_grad[0] = d1p0(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[4].ref_grad[1] = p0(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[4].ref_grad[2] = p0(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[5].ref_grad[0] = d1p1(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[5].ref_grad[1] = p1(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[5].ref_grad[2] = p1(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[6].ref_grad[0] = d1p0(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[6].ref_grad[1] = p0(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[6].ref_grad[2] = p0(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[7].ref_grad[0] = d1p1(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[7].ref_grad[1] = p1(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[7].ref_grad[2] = p1(point[0]) * p1(point[1]) * d1p1(point[2]);
          // edge dofs
          data.phi[ 8].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ 9].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[10].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[10].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[10].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[11].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[11].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[11].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[12].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[12].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[12].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[13].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[13].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[13].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[14].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[14].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[14].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[15].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[15].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[15].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[16].ref_grad[0] = d1p0(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[16].ref_grad[1] = p0(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[16].ref_grad[2] = p0(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[17].ref_grad[0] = d1p1(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[17].ref_grad[1] = p1(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[17].ref_grad[2] = p1(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[18].ref_grad[0] = d1p0(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[18].ref_grad[1] = p0(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[18].ref_grad[2] = p0(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[19].ref_grad[0] = d1p1(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[19].ref_grad[1] = p1(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[19].ref_grad[2] = p1(point[0]) * p1(point[1]) * d1p2(point[2]);
          // face dofs
          data.phi[20].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[20].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[20].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[21].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[21].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[21].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[22].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[22].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[22].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[23].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[23].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[23].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[24].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[24].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[24].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[25].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[25].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[25].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p2(point[2]);
          // center dof
          data.phi[26].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[26].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[26].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p2(point[2]);
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
          using namespace Lagrange2::Intern;
          // vertex dofs
          data.phi[ 0].ref_hess[0][0] = d2p0(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 0].ref_hess[1][1] = p0(point[0]) * d2p0(point[1]) * p0(point[2]);
          data.phi[ 0].ref_hess[2][2] = p0(point[0]) * p0(point[1]) * d2p0(point[2]);
          data.phi[ 0].ref_hess[1][0] =
          data.phi[ 0].ref_hess[0][1] = d1p0(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ 0].ref_hess[2][0] =
          data.phi[ 0].ref_hess[0][2] = d1p0(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ 0].ref_hess[2][1] =
          data.phi[ 0].ref_hess[1][2] = p0(point[0]) * d1p0(point[1]) * d1p0(point[2]);
          data.phi[ 1].ref_hess[0][0] = d2p1(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 1].ref_hess[1][1] = p1(point[0]) * d2p0(point[1]) * p0(point[2]);
          data.phi[ 1].ref_hess[2][2] = p1(point[0]) * p0(point[1]) * d2p0(point[2]);
          data.phi[ 1].ref_hess[1][0] =
          data.phi[ 1].ref_hess[0][1] = d1p1(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ 1].ref_hess[2][0] =
          data.phi[ 1].ref_hess[0][2] = d1p1(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ 1].ref_hess[2][1] =
          data.phi[ 1].ref_hess[1][2] = p1(point[0]) * d1p0(point[1]) * d1p0(point[2]);
          data.phi[ 2].ref_hess[0][0] = d2p0(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ 2].ref_hess[1][1] = p0(point[0]) * d2p1(point[1]) * p0(point[2]);
          data.phi[ 2].ref_hess[2][2] = p0(point[0]) * p1(point[1]) * d2p0(point[2]);
          data.phi[ 2].ref_hess[1][0] =
          data.phi[ 2].ref_hess[0][1] = d1p0(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ 2].ref_hess[2][0] =
          data.phi[ 2].ref_hess[0][2] = d1p0(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[ 2].ref_hess[2][1] =
          data.phi[ 2].ref_hess[1][2] = p0(point[0]) * d1p1(point[1]) * d1p0(point[2]);
          data.phi[ 3].ref_hess[0][0] = d2p1(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ 3].ref_hess[1][1] = p1(point[0]) * d2p1(point[1]) * p0(point[2]);
          data.phi[ 3].ref_hess[2][2] = p1(point[0]) * p1(point[1]) * d2p0(point[2]);
          data.phi[ 3].ref_hess[1][0] =
          data.phi[ 3].ref_hess[0][1] = d1p1(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ 3].ref_hess[2][0] =
          data.phi[ 3].ref_hess[0][2] = d1p1(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[ 3].ref_hess[2][1] =
          data.phi[ 3].ref_hess[1][2] = p1(point[0]) * d1p1(point[1]) * d1p0(point[2]);
          data.phi[ 4].ref_hess[0][0] = d2p0(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[ 4].ref_hess[1][1] = p0(point[0]) * d2p0(point[1]) * p1(point[2]);
          data.phi[ 4].ref_hess[2][2] = p0(point[0]) * p0(point[1]) * d2p1(point[2]);
          data.phi[ 4].ref_hess[1][0] =
          data.phi[ 4].ref_hess[0][1] = d1p0(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[ 4].ref_hess[2][0] =
          data.phi[ 4].ref_hess[0][2] = d1p0(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[ 4].ref_hess[2][1] =
          data.phi[ 4].ref_hess[1][2] = p0(point[0]) * d1p0(point[1]) * d1p1(point[2]);
          data.phi[ 5].ref_hess[0][0] = d2p1(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[ 5].ref_hess[1][1] = p1(point[0]) * d2p0(point[1]) * p1(point[2]);
          data.phi[ 5].ref_hess[2][2] = p1(point[0]) * p0(point[1]) * d2p1(point[2]);
          data.phi[ 5].ref_hess[1][0] =
          data.phi[ 5].ref_hess[0][1] = d1p1(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[ 5].ref_hess[2][0] =
          data.phi[ 5].ref_hess[0][2] = d1p1(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[ 5].ref_hess[2][1] =
          data.phi[ 5].ref_hess[1][2] = p1(point[0]) * d1p0(point[1]) * d1p1(point[2]);
          data.phi[ 6].ref_hess[0][0] = d2p0(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[ 6].ref_hess[1][1] = p0(point[0]) * d2p1(point[1]) * p1(point[2]);
          data.phi[ 6].ref_hess[2][2] = p0(point[0]) * p1(point[1]) * d2p1(point[2]);
          data.phi[ 6].ref_hess[1][0] =
          data.phi[ 6].ref_hess[0][1] = d1p0(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[ 6].ref_hess[2][0] =
          data.phi[ 6].ref_hess[0][2] = d1p0(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[ 6].ref_hess[2][1] =
          data.phi[ 6].ref_hess[1][2] = p0(point[0]) * d1p1(point[1]) * d1p1(point[2]);
          data.phi[ 7].ref_hess[0][0] = d2p1(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[ 7].ref_hess[1][1] = p1(point[0]) * d2p1(point[1]) * p1(point[2]);
          data.phi[ 7].ref_hess[2][2] = p1(point[0]) * p1(point[1]) * d2p1(point[2]);
          data.phi[ 7].ref_hess[1][0] =
          data.phi[ 7].ref_hess[0][1] = d1p1(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[ 7].ref_hess[2][0] =
          data.phi[ 7].ref_hess[0][2] = d1p1(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[ 7].ref_hess[2][1] =
          data.phi[ 7].ref_hess[1][2] = p1(point[0]) * d1p1(point[1]) * d1p1(point[2]);

          // edge dofs
          data.phi[ 8].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_hess[2][2] = p2(point[0]) * p0(point[1]) * d2p0(point[2]);
          data.phi[ 8].ref_hess[1][0] =
          data.phi[ 8].ref_hess[0][1] = d1p2(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_hess[2][0] =
          data.phi[ 8].ref_hess[0][2] = d1p2(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ 8].ref_hess[2][1] =
          data.phi[ 8].ref_hess[1][2] = p2(point[0]) * d1p0(point[1]) * d1p0(point[2]);
          data.phi[ 9].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_hess[2][2] = p2(point[0]) * p1(point[1]) * d2p0(point[2]);
          data.phi[ 9].ref_hess[1][0] =
          data.phi[ 9].ref_hess[0][1] = d1p2(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_hess[2][0] =
          data.phi[ 9].ref_hess[0][2] = d1p2(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[ 9].ref_hess[2][1] =
          data.phi[ 9].ref_hess[1][2] = p2(point[0]) * d1p1(point[1]) * d1p0(point[2]);
          data.phi[10].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[10].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]) * p1(point[2]);
          data.phi[10].ref_hess[2][2] = p2(point[0]) * p0(point[1]) * d2p1(point[2]);
          data.phi[10].ref_hess[1][0] =
          data.phi[10].ref_hess[0][1] = d1p2(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[10].ref_hess[2][0] =
          data.phi[10].ref_hess[0][2] = d1p2(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[10].ref_hess[2][1] =
          data.phi[10].ref_hess[1][2] = p2(point[0]) * d1p0(point[1]) * d1p1(point[2]);
          data.phi[11].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[11].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]) * p1(point[2]);
          data.phi[11].ref_hess[2][2] = p2(point[0]) * p1(point[1]) * d2p1(point[2]);
          data.phi[11].ref_hess[1][0] =
          data.phi[11].ref_hess[0][1] = d1p2(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[11].ref_hess[2][0] =
          data.phi[11].ref_hess[0][2] = d1p2(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[11].ref_hess[2][1] =
          data.phi[11].ref_hess[1][2] = p2(point[0]) * d1p1(point[1]) * d1p1(point[2]);
          data.phi[12].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[12].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]) * p0(point[2]);
          data.phi[12].ref_hess[2][2] = p0(point[0]) * p2(point[1]) * d2p0(point[2]);
          data.phi[12].ref_hess[1][0] =
          data.phi[12].ref_hess[0][1] = d1p0(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[12].ref_hess[2][0] =
          data.phi[12].ref_hess[0][2] = d1p0(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[12].ref_hess[2][1] =
          data.phi[12].ref_hess[1][2] = p0(point[0]) * d1p2(point[1]) * d1p0(point[2]);
          data.phi[13].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[13].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]) * p0(point[2]);
          data.phi[13].ref_hess[2][2] = p1(point[0]) * p2(point[1]) * d2p0(point[2]);
          data.phi[13].ref_hess[1][0] =
          data.phi[13].ref_hess[0][1] = d1p1(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[13].ref_hess[2][0] =
          data.phi[13].ref_hess[0][2] = d1p1(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[13].ref_hess[2][1] =
          data.phi[13].ref_hess[1][2] = p1(point[0]) * d1p2(point[1]) * d1p0(point[2]);
          data.phi[14].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[14].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]) * p1(point[2]);
          data.phi[14].ref_hess[2][2] = p0(point[0]) * p2(point[1]) * d2p1(point[2]);
          data.phi[14].ref_hess[1][0] =
          data.phi[14].ref_hess[0][1] = d1p0(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[14].ref_hess[2][0] =
          data.phi[14].ref_hess[0][2] = d1p0(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[14].ref_hess[2][1] =
          data.phi[14].ref_hess[1][2] = p0(point[0]) * d1p2(point[1]) * d1p1(point[2]);
          data.phi[15].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[15].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]) * p1(point[2]);
          data.phi[15].ref_hess[2][2] = p1(point[0]) * p2(point[1]) * d2p1(point[2]);
          data.phi[15].ref_hess[1][0] =
          data.phi[15].ref_hess[0][1] = d1p1(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[15].ref_hess[2][0] =
          data.phi[15].ref_hess[0][2] = d1p1(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[15].ref_hess[2][1] =
          data.phi[15].ref_hess[1][2] = p1(point[0]) * d1p2(point[1]) * d1p1(point[2]);
          data.phi[16].ref_hess[0][0] = d2p0(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[16].ref_hess[1][1] = p0(point[0]) * d2p0(point[1]) * p2(point[2]);
          data.phi[16].ref_hess[2][2] = p0(point[0]) * p0(point[1]) * d2p2(point[2]);
          data.phi[16].ref_hess[1][0] =
          data.phi[16].ref_hess[0][1] = d1p0(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[16].ref_hess[2][0] =
          data.phi[16].ref_hess[0][2] = d1p0(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[16].ref_hess[2][1] =
          data.phi[16].ref_hess[1][2] = p0(point[0]) * d1p0(point[1]) * d1p2(point[2]);
          data.phi[17].ref_hess[0][0] = d2p1(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[17].ref_hess[1][1] = p1(point[0]) * d2p0(point[1]) * p2(point[2]);
          data.phi[17].ref_hess[2][2] = p1(point[0]) * p0(point[1]) * d2p2(point[2]);
          data.phi[17].ref_hess[1][0] =
          data.phi[17].ref_hess[0][1] = d1p1(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[17].ref_hess[2][0] =
          data.phi[17].ref_hess[0][2] = d1p1(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[17].ref_hess[2][1] =
          data.phi[17].ref_hess[1][2] = p1(point[0]) * d1p0(point[1]) * d1p2(point[2]);
          data.phi[18].ref_hess[0][0] = d2p0(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[18].ref_hess[1][1] = p0(point[0]) * d2p1(point[1]) * p2(point[2]);
          data.phi[18].ref_hess[2][2] = p0(point[0]) * p1(point[1]) * d2p2(point[2]);
          data.phi[18].ref_hess[1][0] =
          data.phi[18].ref_hess[0][1] = d1p0(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[18].ref_hess[2][0] =
          data.phi[18].ref_hess[0][2] = d1p0(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[18].ref_hess[2][1] =
          data.phi[18].ref_hess[1][2] = p0(point[0]) * d1p1(point[1]) * d1p2(point[2]);
          data.phi[19].ref_hess[0][0] = d2p1(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[19].ref_hess[1][1] = p1(point[0]) * d2p1(point[1]) * p2(point[2]);
          data.phi[19].ref_hess[2][2] = p1(point[0]) * p1(point[1]) * d2p2(point[2]);
          data.phi[19].ref_hess[1][0] =
          data.phi[19].ref_hess[0][1] = d1p1(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[19].ref_hess[2][0] =
          data.phi[19].ref_hess[0][2] = d1p1(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[19].ref_hess[2][1] =
          data.phi[19].ref_hess[1][2] = p1(point[0]) * d1p1(point[1]) * d1p2(point[2]);

          // face dofs
          data.phi[20].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[20].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]) * p0(point[2]);
          data.phi[20].ref_hess[2][2] = p2(point[0]) * p2(point[1]) * d2p0(point[2]);
          data.phi[20].ref_hess[1][0] =
          data.phi[20].ref_hess[0][1] = d1p2(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[20].ref_hess[2][0] =
          data.phi[20].ref_hess[0][2] = d1p2(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[20].ref_hess[2][1] =
          data.phi[20].ref_hess[1][2] = p2(point[0]) * d1p2(point[1]) * d1p0(point[2]);
          data.phi[21].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[21].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]) * p1(point[2]);
          data.phi[21].ref_hess[2][2] = p2(point[0]) * p2(point[1]) * d2p1(point[2]);
          data.phi[21].ref_hess[1][0] =
          data.phi[21].ref_hess[0][1] = d1p2(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[21].ref_hess[2][0] =
          data.phi[21].ref_hess[0][2] = d1p2(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[21].ref_hess[2][1] =
          data.phi[21].ref_hess[1][2] = p2(point[0]) * d1p2(point[1]) * d1p1(point[2]);
          data.phi[22].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[22].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]) * p2(point[2]);
          data.phi[22].ref_hess[2][2] = p2(point[0]) * p0(point[1]) * d2p2(point[2]);
          data.phi[22].ref_hess[1][0] =
          data.phi[22].ref_hess[0][1] = d1p2(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[22].ref_hess[2][0] =
          data.phi[22].ref_hess[0][2] = d1p2(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[22].ref_hess[2][1] =
          data.phi[22].ref_hess[1][2] = p2(point[0]) * d1p0(point[1]) * d1p2(point[2]);
          data.phi[23].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[23].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]) * p2(point[2]);
          data.phi[23].ref_hess[2][2] = p2(point[0]) * p1(point[1]) * d2p2(point[2]);
          data.phi[23].ref_hess[1][0] =
          data.phi[23].ref_hess[0][1] = d1p2(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[23].ref_hess[2][0] =
          data.phi[23].ref_hess[0][2] = d1p2(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[23].ref_hess[2][1] =
          data.phi[23].ref_hess[1][2] = p2(point[0]) * d1p1(point[1]) * d1p2(point[2]);
          data.phi[24].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[24].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]) * p2(point[2]);
          data.phi[24].ref_hess[2][2] = p0(point[0]) * p2(point[1]) * d2p2(point[2]);
          data.phi[24].ref_hess[1][0] =
          data.phi[24].ref_hess[0][1] = d1p0(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[24].ref_hess[2][0] =
          data.phi[24].ref_hess[0][2] = d1p0(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[24].ref_hess[2][1] =
          data.phi[24].ref_hess[1][2] = p0(point[0]) * d1p2(point[1]) * d1p2(point[2]);
          data.phi[25].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[25].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]) * p2(point[2]);
          data.phi[25].ref_hess[2][2] = p1(point[0]) * p2(point[1]) * d2p2(point[2]);
          data.phi[25].ref_hess[1][0] =
          data.phi[25].ref_hess[0][1] = d1p1(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[25].ref_hess[2][0] =
          data.phi[25].ref_hess[0][2] = d1p1(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[25].ref_hess[2][1] =
          data.phi[25].ref_hess[1][2] = p1(point[0]) * d1p2(point[1]) * d1p2(point[2]);

          // center dof
          data.phi[26].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[26].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]) * p2(point[2]);
          data.phi[26].ref_hess[2][2] = p2(point[0]) * p2(point[1]) * d2p2(point[2]);
          data.phi[26].ref_hess[1][0] =
          data.phi[26].ref_hess[0][1] = d1p2(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[26].ref_hess[2][0] =
          data.phi[26].ref_hess[0][2] = d1p2(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[26].ref_hess[2][1] =
          data.phi[26].ref_hess[1][2] = p2(point[0]) * d1p2(point[1]) * d1p2(point[2]);
        }
      }; // class Evaluator<...,Hypercube<3>>
    } // namespace Lagrange2
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_LAGRANGE2_EVALUATOR_HPP
