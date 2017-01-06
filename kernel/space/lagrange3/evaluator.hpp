#pragma once
#ifndef KERNEL_SPACE_LAGRANGE3_EVALUATOR_HPP
#define KERNEL_SPACE_LAGRANGE3_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/space/parametric_evaluator.hpp>
#include <kernel/space/dof_mapping_common.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Lagrange3
    {
      /**
       * \brief Lagrange-3 Element Evaluator reference capabilities
       *
       * \author Peter Zajac
       */
      static constexpr SpaceTags ref_caps = SpaceTags::ref_value | SpaceTags::ref_grad | SpaceTags::ref_hess;

      static constexpr SpaceTags ref_caps_hexa = SpaceTags::ref_value | SpaceTags::ref_grad;

      /// \cond internal
      namespace Intern
      {
        // p0, p1, p2 and p3 are the 1D basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the Hypercube evaluators.

        // basis function for left vertex
        template<typename T_>
        inline T_ p0(T_ x)
        {
          return T_(0.0625) * (-T_(1) + x*( T_(1) + T_(9)*x*(T_(1) - x)));
        }

        // basis function for right vertex
        template<typename T_>
        inline T_ p1(T_ x)
        {
          return T_(0.0625) * (-T_(1) + x*(-T_(1) + T_(9)*x*(T_(1) + x)));
        }

        // basis function for first edge point
        template<typename T_>
        inline T_ p2(T_ x)
        {
          return T_(0.5625) * (T_(1) + x*(-T_(3) + x*(-T_(1) + T_(3)*x)));
        }

        // basis function for second edge point
        template<typename T_>
        inline T_ p3(T_ x)
        {
          return T_(0.5625) * (T_(1) + x*( T_(3) + x*(-T_(1) - T_(3)*x)));
        }

        // first order derivatives

        template<typename T_>
        inline T_ d1p0(T_ x)
        {
          return T_(0.0625) * ( T_(1) + x*(T_(18) - T_(27)*x));
        }

        template<typename T_>
        inline T_ d1p1(T_ x)
        {
          return T_(0.0625) * (-T_(1) + x*(T_(18) + T_(27)*x));
        }

        template<typename T_>
        inline T_ d1p2(T_ x)
        {
          return T_(0.0625) * (-T_(27) + x*(-T_(18) + T_(81)*x));
        }

        template<typename T_>
        inline T_ d1p3(T_ x)
        {
          return T_(0.0625) * ( T_(27) + x*(-T_(18) - T_(81)*x));
        }

        // second order derivatives

        template<typename T_>
        inline T_ d2p0(T_ x)
        {
          return T_(1.125) * (T_(1) - T_(3)*x);
        }

        template<typename T_>
        inline T_ d2p1(T_ x)
        {
          return T_(1.125) * (T_(1) + T_(3)*x);
        }

        template<typename T_>
        inline T_ d2p2(T_ x)
        {
          return T_(1.125) * (-T_(1) + T_(9)*x);
        }

        template<typename T_>
        inline T_ d2p3(T_ x)
        {
          return T_(1.125) * (-T_(1) - T_(9)*x);
        }
      } // namespace Intern
      /// \endcond

      /**
       * \brief Lagrange-3 Element Evaluator class template declaration.
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
       * \brief Lagrange-3 Element evaluator implementation for 1D Hypercube shape
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
        void eval_ref_values(EvalData_& data, const DomainPointType& point) const
        {
          data.phi[0].ref_value = Intern::p0(point[0]);
          data.phi[1].ref_value = Intern::p1(point[0]);
          data.phi[2].ref_value = Intern::p2(point[0]);
          data.phi[3].ref_value = Intern::p3(point[0]);
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
        void eval_ref_gradients(EvalData_& data, const DomainPointType& point) const
        {
          data.phi[0].ref_grad[0] = Intern::d1p0(point[0]);
          data.phi[1].ref_grad[0] = Intern::d1p1(point[0]);
          data.phi[2].ref_grad[0] = Intern::d1p2(point[0]);
          data.phi[3].ref_grad[0] = Intern::d1p3(point[0]);
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
        void eval_ref_hessians(EvalData_& data, const DomainPointType& point) const
        {
          data.phi[0].ref_hess[0][0] = Intern::d2p0(point[0]);
          data.phi[1].ref_hess[0][0] = Intern::d2p1(point[0]);
          data.phi[2].ref_hess[0][0] = Intern::d2p2(point[0]);
          data.phi[3].ref_hess[0][0] = Intern::d2p3(point[0]);
        }
      }; // class Evaluator<...,Hypercube<1>>

      /**
       * \brief Lagrange-3 Element evaluator implementation for Quadrilateral shape
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

      protected:
        /// edge dof indices
        int ek[3][2];

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

        void NOINLINE prepare(const TrafoEvaluator_& trafo_eval)
        {
          // compute edge orientations
          Geometry::Intern::SubIndexMapping<Shape::Simplex<2>, 1, 0> sim(
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,0>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,1>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<1,0>());

          // fetch edge dof indices
          for(int i(0); i < 3; ++i)
          {
            for(int j(0); j < 2; ++j)
            {
              ek[i][j] = 3 + 2*i + int(sim.map(i, j));
            }
          }
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
        void NOINLINE eval_ref_values(EvalData_& data, const DomainPointType& point) const
        {
          // vertex dofs
          data.phi[0].ref_value = - DataType(4.5) * (point[0] + point[1] - (DataType(1) / DataType(3))) * (point[0] + point[1] - (DataType(2) / DataType(3))) * (point[0] + point[1] - DataType(1));
          data.phi[1].ref_value = DataType(4.5) * point[0] * (point[0] - (DataType(1) / DataType(3))) * (point[0] - (DataType(2) / DataType(3)));
          data.phi[2].ref_value = DataType(4.5) * point[1] * (point[1] - (DataType(1) / DataType(3))) * (point[1] - (DataType(2) / DataType(3)));
          // edge dofs
          data.phi[ek[0][0]].ref_value = DataType(13.5) * point[0] * point[1] * (point[0] - (DataType(1) / DataType(3)));
          data.phi[ek[0][1]].ref_value = DataType(13.5) * point[0] * point[1] * (point[1] - (DataType(1) / DataType(3)));
          data.phi[ek[1][0]].ref_value = - DataType(13.5) * point[1] * (point[1] - (DataType(1) / DataType(3))) * (point[0] + point[1] - DataType(1));
          data.phi[ek[1][1]].ref_value = DataType(13.5) * point[1] * (point[0] + point[1] - DataType(1)) * (point[0] + point[1] - (DataType(2) / DataType(3)));
          data.phi[ek[2][0]].ref_value = DataType(13.5) * point[0] * (point[0] + point[1] - DataType(1)) * (point[0] + point[1] - (DataType(2) / DataType(3)));
          data.phi[ek[2][1]].ref_value = - DataType(13.5) * point[0] * (point[0] - (DataType(1) / DataType(3))) * (point[0] + point[1] - DataType(1));
          // center dofs
          data.phi[9].ref_value = - DataType(27) * point[0] * point[1] * (point[0] + point[1] - DataType(1));
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
        void NOINLINE eval_ref_gradients(EvalData_& data, const DomainPointType& point) const
        {
          // vertex dofs
          data.phi[0].ref_grad[0] = - DataType(0.5) * (DataType(27) * point[0] * point[0] + (DataType(54) * point[1] - DataType(36)) * point[0] + DataType(27) * point[1] * point[1] - DataType(36) * point[1] + DataType(11));
          data.phi[0].ref_grad[1] = - DataType(0.5) * (DataType(27) * point[1] * point[1] + (DataType(54) * point[0] - DataType(36)) * point[1] + DataType(27) * point[0] * point[0] - DataType(36) * point[0] + DataType(11));
          data.phi[1].ref_grad[0] = DataType(0.5) * (DataType(27) * point[0] * point[0] - DataType(18) * point[0] + DataType(2));
          data.phi[1].ref_grad[1] = DataType(0);
          data.phi[2].ref_grad[0] = DataType(0);
          data.phi[2].ref_grad[1] = DataType(0.5) * (DataType(27) * point[1] * point[1] - DataType(18) * point[1] + DataType(2));
          // edge dofs
          data.phi[ek[0][0]].ref_grad[0] = DataType(0.5) * (DataType(9) * point[1] * (DataType(6) * point[0] - DataType(1)));
          data.phi[ek[0][0]].ref_grad[1] = DataType(0.5) * (DataType(9) * point[0] * (DataType(3) * point[0] - DataType(1)));
          data.phi[ek[0][1]].ref_grad[0] = DataType(0.5) * (DataType(9) * point[1] * (DataType(3) * point[1] - DataType(1)));
          data.phi[ek[0][1]].ref_grad[1] = DataType(0.5) * (DataType(9) * point[0] * (DataType(6) * point[1] - DataType(1)));
          data.phi[ek[1][0]].ref_grad[0] = DataType(0.5) * (- DataType(9) * point[1] * (DataType(3) * point[1] - DataType(1)));
          data.phi[ek[1][0]].ref_grad[1] = DataType(0.5) * (- DataType(81) * point[1] * point[1] - (DataType(54) * point[0] - DataType(72)) * point[1] + DataType(9) * point[0] - DataType(9));
          data.phi[ek[1][1]].ref_grad[0] = DataType(0.5) * (DataType(9) * point[1] * (DataType(6) * point[0] + DataType(6) * point[1] - DataType(5)));
          data.phi[ek[1][1]].ref_grad[1] = DataType(0.5) * (DataType(81) * point[1] * point[1] + ( DataType(108) * point[0] - DataType(90)) * point[1] + DataType(27) * point[0] * point[0] - DataType(45) * point[0] + DataType(18));
          data.phi[ek[2][0]].ref_grad[0] = DataType(0.5) * (DataType(81) * point[0] * point[0] + ( DataType(108) * point[1] - DataType(90)) * point[0] + DataType(27) * point[1] * point[1] - DataType(45) * point[1] + DataType(18));
          data.phi[ek[2][0]].ref_grad[1] = DataType(0.5) * (DataType(9) * point[0] * (DataType(6) * point[1] + DataType(6) * point[0] - DataType(5)));
          data.phi[ek[2][1]].ref_grad[0] = DataType(0.5) * (- DataType(81) * point[0] * point[0] - (DataType(54) * point[1] - DataType(72)) * point[0] + DataType(9) * point[1] - DataType(9));
          data.phi[ek[2][1]].ref_grad[1] = DataType(0.5) * (DataType(27) * (DataType(1) / DataType(3) - point[0]) * point[0]);
          // center dofs
          data.phi[9].ref_grad[0] = - DataType(27) * point[1] * (DataType(2) * point[0] + point[1] - DataType(1));
          data.phi[9].ref_grad[1] = - DataType(27) * point[0] * (DataType(2) * point[1] + point[0] - DataType(1));
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
        void NOINLINE eval_ref_hessians(EvalData_& data, const DomainPointType& point) const
        {
          // vertex dofs
          data.phi[0].ref_hess[0][0] = - DataType(27) * (point[0] + point[1]) + DataType(18); // dxx
          data.phi[0].ref_hess[1][1] = - DataType(27) * (point[1] + point[0]) + DataType(18); // dyy
          data.phi[0].ref_hess[1][0] =
          data.phi[0].ref_hess[0][1] = - DataType(27) * (point[1] + point[0]) + DataType(18); // dxy
          data.phi[1].ref_hess[0][0] = DataType(27) * point[0] - DataType(9);
          data.phi[1].ref_hess[1][1] = DataType(0);
          data.phi[1].ref_hess[1][0] =
          data.phi[1].ref_hess[0][1] = DataType(0);
          data.phi[2].ref_hess[0][0] = DataType(0);
          data.phi[2].ref_hess[1][1] = DataType(27) * point[1] - DataType(9);
          data.phi[2].ref_hess[1][0] =
          data.phi[2].ref_hess[0][1] = DataType(0);
          // edge dofs
          data.phi[ek[0][0]].ref_hess[0][0] = DataType(27) * point[1];
          data.phi[ek[0][0]].ref_hess[1][1] = DataType(0);
          data.phi[ek[0][0]].ref_hess[0][1] =
          data.phi[ek[0][0]].ref_hess[1][0] = DataType(0.5) * (DataType(54) * point[0] - DataType(9));
          data.phi[ek[0][1]].ref_hess[0][0] = DataType(0);
          data.phi[ek[0][1]].ref_hess[1][1] = DataType(27) * point[0];
          data.phi[ek[0][1]].ref_hess[0][1] =
          data.phi[ek[0][1]].ref_hess[1][0] = DataType(0.5) * (DataType(54) * point[1] - DataType(9));
          data.phi[ek[1][0]].ref_hess[0][0] = DataType(0);
          data.phi[ek[1][0]].ref_hess[1][1] = - DataType(81) * point[1] - DataType(27) * point[0] + DataType(36);
          data.phi[ek[1][0]].ref_hess[0][1] =
          data.phi[ek[1][0]].ref_hess[1][0] = - DataType(0.5) * (DataType(54) * point[1] - DataType(9));
          data.phi[ek[1][1]].ref_hess[0][0] = DataType(27) * point[1];
          data.phi[ek[1][1]].ref_hess[1][1] = DataType(81) * point[1] + DataType(54) * point[0] - DataType(45);
          data.phi[ek[1][1]].ref_hess[0][1] =
          data.phi[ek[1][1]].ref_hess[1][0] = DataType(0.5) * (DataType(108) * point[1] + DataType(54) * point[0] - DataType(45));
          data.phi[ek[2][0]].ref_hess[0][0] = DataType(81) * point[0] + DataType(54) * point[1] - DataType(45);
          data.phi[ek[2][0]].ref_hess[1][1] = DataType(27) * point[0];
          data.phi[ek[2][0]].ref_hess[0][1] =
          data.phi[ek[2][0]].ref_hess[1][0] = DataType(0.5) * (DataType(108) * point[0] + DataType(54) * point[1] - DataType(45));
          data.phi[ek[2][1]].ref_hess[0][0] = - DataType(81) * point[0] - DataType(27) * point[1] + DataType(36);
          data.phi[ek[2][1]].ref_hess[1][1] = DataType(0);
          data.phi[ek[2][1]].ref_hess[0][1] =
          data.phi[ek[2][1]].ref_hess[1][0] = - DataType(0.5) * (DataType(54) * point[0] - DataType(9));
          // center dofs
          data.phi[9].ref_hess[0][0] = - DataType(54) * point[1];
          data.phi[9].ref_hess[1][1] = - DataType(54) * point[0];
          data.phi[9].ref_hess[0][1] =
          data.phi[9].ref_hess[1][0] = DataType(27) * (DataType(1) - DataType(2)*(point[0] + point[1]));
        }
      }; // class Evaluator<...,Simplex<2>>

      /**
       * \brief Lagrange-3 Element evaluator implementation for Quadrilateral shape
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

      protected:
        /// edge dof indices
        int ek[4][2];

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
          return 16;
        }

        void NOINLINE prepare(const TrafoEvaluator_& trafo_eval)
        {
          // compute edge orientations
          Geometry::Intern::SubIndexMapping<Shape::Hypercube<2>, 1, 0> sim(
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,0>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,1>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<1,0>());

          // fetch edge dof indices
          for(int i(0); i < 4; ++i)
          {
            for(int j(0); j < 2; ++j)
            {
              ek[i][j] = 4 + 2*i + int(sim.map(i, j));
            }
          }
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
        void NOINLINE eval_ref_values(EvalData_& data, const DomainPointType& point) const
        {
          using namespace Lagrange3::Intern;

          // vertex dofs
          data.phi[0].ref_value = p0(point[0]) * p0(point[1]);
          data.phi[1].ref_value = p1(point[0]) * p0(point[1]);
          data.phi[2].ref_value = p0(point[0]) * p1(point[1]);
          data.phi[3].ref_value = p1(point[0]) * p1(point[1]);
          // edge dofs
          data.phi[ek[0][0]].ref_value = p2(point[0]) * p0(point[1]);
          data.phi[ek[0][1]].ref_value = p3(point[0]) * p0(point[1]);
          data.phi[ek[1][0]].ref_value = p2(point[0]) * p1(point[1]);
          data.phi[ek[1][1]].ref_value = p3(point[0]) * p1(point[1]);
          data.phi[ek[2][0]].ref_value = p0(point[0]) * p2(point[1]);
          data.phi[ek[2][1]].ref_value = p0(point[0]) * p3(point[1]);
          data.phi[ek[3][0]].ref_value = p1(point[0]) * p2(point[1]);
          data.phi[ek[3][1]].ref_value = p1(point[0]) * p3(point[1]);
          // center dofs
          data.phi[12].ref_value = p2(point[0]) * p2(point[1]);
          data.phi[13].ref_value = p3(point[0]) * p2(point[1]);
          data.phi[14].ref_value = p2(point[0]) * p3(point[1]);
          data.phi[15].ref_value = p3(point[0]) * p3(point[1]);
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
        void NOINLINE eval_ref_gradients(EvalData_& data, const DomainPointType& point) const
        {
          using namespace Lagrange3::Intern;

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
          data.phi[ek[0][0]].ref_grad[0] = d1p2(point[0]) * p0(point[1]);
          data.phi[ek[0][0]].ref_grad[1] = p2(point[0]) * d1p0(point[1]);
          data.phi[ek[0][1]].ref_grad[0] = d1p3(point[0]) * p0(point[1]);
          data.phi[ek[0][1]].ref_grad[1] = p3(point[0]) * d1p0(point[1]);
          data.phi[ek[1][0]].ref_grad[0] = d1p2(point[0]) * p1(point[1]);
          data.phi[ek[1][0]].ref_grad[1] = p2(point[0]) * d1p1(point[1]);
          data.phi[ek[1][1]].ref_grad[0] = d1p3(point[0]) * p1(point[1]);
          data.phi[ek[1][1]].ref_grad[1] = p3(point[0]) * d1p1(point[1]);
          data.phi[ek[2][0]].ref_grad[0] = d1p0(point[0]) * p2(point[1]);
          data.phi[ek[2][0]].ref_grad[1] = p0(point[0]) * d1p2(point[1]);
          data.phi[ek[2][1]].ref_grad[0] = d1p0(point[0]) * p3(point[1]);
          data.phi[ek[2][1]].ref_grad[1] = p0(point[0]) * d1p3(point[1]);
          data.phi[ek[3][0]].ref_grad[0] = d1p1(point[0]) * p2(point[1]);
          data.phi[ek[3][0]].ref_grad[1] = p1(point[0]) * d1p2(point[1]);
          data.phi[ek[3][1]].ref_grad[0] = d1p1(point[0]) * p3(point[1]);
          data.phi[ek[3][1]].ref_grad[1] = p1(point[0]) * d1p3(point[1]);
          // center dofs
          data.phi[12].ref_grad[0] = d1p2(point[0]) * p2(point[1]);
          data.phi[12].ref_grad[1] = p2(point[0]) * d1p2(point[1]);
          data.phi[13].ref_grad[0] = d1p3(point[0]) * p2(point[1]);
          data.phi[13].ref_grad[1] = p3(point[0]) * d1p2(point[1]);
          data.phi[14].ref_grad[0] = d1p2(point[0]) * p3(point[1]);
          data.phi[14].ref_grad[1] = p2(point[0]) * d1p3(point[1]);
          data.phi[15].ref_grad[0] = d1p3(point[0]) * p3(point[1]);
          data.phi[15].ref_grad[1] = p3(point[0]) * d1p3(point[1]);
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
        void NOINLINE eval_ref_hessians(EvalData_& data, const DomainPointType& point) const
        {
          using namespace Lagrange3::Intern;

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
          data.phi[ek[0][0]].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]);
          data.phi[ek[0][0]].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]);
          data.phi[ek[0][0]].ref_hess[0][1] =
          data.phi[ek[0][0]].ref_hess[1][0] = d1p2(point[0]) * d1p0(point[1]);
          data.phi[ek[0][1]].ref_hess[0][0] = d2p3(point[0]) * p0(point[1]);
          data.phi[ek[0][1]].ref_hess[1][1] = p3(point[0]) * d2p0(point[1]);
          data.phi[ek[0][1]].ref_hess[0][1] =
          data.phi[ek[0][1]].ref_hess[1][0] = d1p3(point[0]) * d1p0(point[1]);
          data.phi[ek[1][0]].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]);
          data.phi[ek[1][0]].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]);
          data.phi[ek[1][0]].ref_hess[0][1] =
          data.phi[ek[1][0]].ref_hess[1][0] = d1p2(point[0]) * d1p1(point[1]);
          data.phi[ek[1][1]].ref_hess[0][0] = d2p3(point[0]) * p1(point[1]);
          data.phi[ek[1][1]].ref_hess[1][1] = p3(point[0]) * d2p1(point[1]);
          data.phi[ek[1][1]].ref_hess[0][1] =
          data.phi[ek[1][1]].ref_hess[1][0] = d1p3(point[0]) * d1p1(point[1]);
          data.phi[ek[2][0]].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]);
          data.phi[ek[2][0]].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]);
          data.phi[ek[2][0]].ref_hess[0][1] =
          data.phi[ek[2][0]].ref_hess[1][0] = d1p0(point[0]) * d1p2(point[1]);
          data.phi[ek[2][1]].ref_hess[0][0] = d2p0(point[0]) * p3(point[1]);
          data.phi[ek[2][1]].ref_hess[1][1] = p0(point[0]) * d2p3(point[1]);
          data.phi[ek[2][1]].ref_hess[0][1] =
          data.phi[ek[2][1]].ref_hess[1][0] = d1p0(point[0]) * d1p3(point[1]);
          data.phi[ek[3][0]].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]);
          data.phi[ek[3][0]].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]);
          data.phi[ek[3][0]].ref_hess[0][1] =
          data.phi[ek[3][0]].ref_hess[1][0] = d1p1(point[0]) * d1p2(point[1]);
          data.phi[ek[3][1]].ref_hess[0][0] = d2p1(point[0]) * p3(point[1]);
          data.phi[ek[3][1]].ref_hess[1][1] = p1(point[0]) * d2p3(point[1]);
          data.phi[ek[3][1]].ref_hess[0][1] =
          data.phi[ek[3][1]].ref_hess[1][0] = d1p1(point[0]) * d1p3(point[1]);
          // center dofs
          data.phi[12].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]);
          data.phi[12].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]);
          data.phi[12].ref_hess[0][1] =
          data.phi[12].ref_hess[1][0] = d1p2(point[0]) * d1p2(point[1]);
          data.phi[13].ref_hess[0][0] = d2p3(point[0]) * p2(point[1]);
          data.phi[13].ref_hess[1][1] = p3(point[0]) * d2p2(point[1]);
          data.phi[13].ref_hess[0][1] =
          data.phi[13].ref_hess[1][0] = d1p3(point[0]) * d1p2(point[1]);
          data.phi[14].ref_hess[0][0] = d2p2(point[0]) * p3(point[1]);
          data.phi[14].ref_hess[1][1] = p2(point[0]) * d2p3(point[1]);
          data.phi[14].ref_hess[0][1] =
          data.phi[14].ref_hess[1][0] = d1p2(point[0]) * d1p3(point[1]);
          data.phi[15].ref_hess[0][0] = d2p3(point[0]) * p3(point[1]);
          data.phi[15].ref_hess[1][1] = p3(point[0]) * d2p3(point[1]);
          data.phi[15].ref_hess[0][1] =
          data.phi[15].ref_hess[1][0] = d1p3(point[0]) * d1p3(point[1]);
        }
      }; // class Evaluator<...,Hypercube<2>>

      /**
       * \brief Lagrange-3 Element evaluator implementation for Hexahedron shape
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
          ref_caps_hexa>
      {
      public:
        /// base-class typedef
        typedef ParametricEvaluator<Evaluator, TrafoEvaluator_, SpaceEvalTraits_, ref_caps_hexa> BaseClass;

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

      protected:
        /// edge dof indices
        int ek[12][2];
        /// quad dof indices
        int qk[6][4];

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
          return 64;
        }

        void NOINLINE prepare(const TrafoEvaluator_& trafo_eval)
        {
          // compute edge orientations
          Geometry::Intern::SubIndexMapping<Shape::Hypercube<3>, 1, 0> sim1(
            trafo_eval.get_trafo().get_mesh().template get_index_set<3,0>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<3,1>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<1,0>());

          // compute quad orientations
          Geometry::Intern::SubIndexMapping<Shape::Hypercube<3>, 2, 0> sim2(
            trafo_eval.get_trafo().get_mesh().template get_index_set<3,0>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<3,2>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,0>());

          // fetch edge dof indices (9,...,31)
          for(int i(0); i < 12; ++i)
          {
            for(int j(0); j < 2; ++j)
            {
              ek[i][j] = 8 + 2*i + int(sim1.map(i, j));
            }
          }

          // fetch quad dof indices (32,...,55)
          for(int i(0); i < 6; ++i)
          {
            for(int j(0); j < 4; ++j)
            {
              qk[i][j] = 32 + 4*i + int(sim2.map(i, j));
            }
          }
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
        void NOINLINE eval_ref_values(EvalData_& data, const DomainPointType& point) const
        {
          using namespace Lagrange3::Intern;

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
          data.phi[ek[ 0][0]].ref_value = p2(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ek[ 0][1]].ref_value = p3(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ek[ 1][0]].ref_value = p2(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ek[ 1][1]].ref_value = p3(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ek[ 2][0]].ref_value = p2(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[ek[ 2][1]].ref_value = p3(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[ek[ 3][0]].ref_value = p2(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[ek[ 3][1]].ref_value = p3(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[ek[ 4][0]].ref_value = p0(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[ek[ 4][1]].ref_value = p0(point[0]) * p3(point[1]) * p0(point[2]);
          data.phi[ek[ 5][0]].ref_value = p1(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[ek[ 5][1]].ref_value = p1(point[0]) * p3(point[1]) * p0(point[2]);
          data.phi[ek[ 6][0]].ref_value = p0(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[ek[ 6][1]].ref_value = p0(point[0]) * p3(point[1]) * p1(point[2]);
          data.phi[ek[ 7][0]].ref_value = p1(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[ek[ 7][1]].ref_value = p1(point[0]) * p3(point[1]) * p1(point[2]);
          data.phi[ek[ 8][0]].ref_value = p0(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[ek[ 8][1]].ref_value = p0(point[0]) * p0(point[1]) * p3(point[2]);
          data.phi[ek[ 9][0]].ref_value = p1(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[ek[ 9][1]].ref_value = p1(point[0]) * p0(point[1]) * p3(point[2]);
          data.phi[ek[10][0]].ref_value = p0(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[ek[10][1]].ref_value = p0(point[0]) * p1(point[1]) * p3(point[2]);
          data.phi[ek[11][0]].ref_value = p1(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[ek[11][1]].ref_value = p1(point[0]) * p1(point[1]) * p3(point[2]);
          // face dofs
          data.phi[qk[0][0]].ref_value = p2(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[qk[0][1]].ref_value = p3(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[qk[0][2]].ref_value = p2(point[0]) * p3(point[1]) * p0(point[2]);
          data.phi[qk[0][3]].ref_value = p3(point[0]) * p3(point[1]) * p0(point[2]);
          data.phi[qk[1][0]].ref_value = p2(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[qk[1][1]].ref_value = p3(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[qk[1][2]].ref_value = p2(point[0]) * p3(point[1]) * p1(point[2]);
          data.phi[qk[1][3]].ref_value = p3(point[0]) * p3(point[1]) * p1(point[2]);
          data.phi[qk[2][0]].ref_value = p2(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[qk[2][1]].ref_value = p3(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[qk[2][2]].ref_value = p2(point[0]) * p0(point[1]) * p3(point[2]);
          data.phi[qk[2][3]].ref_value = p3(point[0]) * p0(point[1]) * p3(point[2]);
          data.phi[qk[3][0]].ref_value = p2(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[qk[3][1]].ref_value = p3(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[qk[3][2]].ref_value = p2(point[0]) * p1(point[1]) * p3(point[2]);
          data.phi[qk[3][3]].ref_value = p3(point[0]) * p1(point[1]) * p3(point[2]);
          data.phi[qk[4][0]].ref_value = p0(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[qk[4][1]].ref_value = p0(point[0]) * p3(point[1]) * p2(point[2]);
          data.phi[qk[4][2]].ref_value = p0(point[0]) * p2(point[1]) * p3(point[2]);
          data.phi[qk[4][3]].ref_value = p0(point[0]) * p3(point[1]) * p3(point[2]);
          data.phi[qk[5][0]].ref_value = p1(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[qk[5][1]].ref_value = p1(point[0]) * p3(point[1]) * p2(point[2]);
          data.phi[qk[5][2]].ref_value = p1(point[0]) * p2(point[1]) * p3(point[2]);
          data.phi[qk[5][3]].ref_value = p1(point[0]) * p3(point[1]) * p3(point[2]);
          // center dofs
          data.phi[56].ref_value = p2(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[57].ref_value = p3(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[58].ref_value = p2(point[0]) * p3(point[1]) * p2(point[2]);
          data.phi[59].ref_value = p3(point[0]) * p3(point[1]) * p2(point[2]);
          data.phi[60].ref_value = p2(point[0]) * p2(point[1]) * p3(point[2]);
          data.phi[61].ref_value = p3(point[0]) * p2(point[1]) * p3(point[2]);
          data.phi[62].ref_value = p2(point[0]) * p3(point[1]) * p3(point[2]);
          data.phi[63].ref_value = p3(point[0]) * p3(point[1]) * p3(point[2]);
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
        void NOINLINE eval_ref_gradients(EvalData_& data, const DomainPointType& point) const
        {
          using namespace Lagrange3::Intern;

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
          data.phi[ek[ 0][0]].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ek[ 0][0]].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ek[ 0][0]].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ek[ 0][1]].ref_grad[0] = d1p3(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ek[ 0][1]].ref_grad[1] = p3(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ek[ 0][1]].ref_grad[2] = p3(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ek[ 1][0]].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ek[ 1][0]].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ek[ 1][0]].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[ek[ 1][1]].ref_grad[0] = d1p3(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ek[ 1][1]].ref_grad[1] = p3(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ek[ 1][1]].ref_grad[2] = p3(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[ek[ 2][0]].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[ek[ 2][0]].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[ek[ 2][0]].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[ek[ 2][1]].ref_grad[0] = d1p3(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[ek[ 2][1]].ref_grad[1] = p3(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[ek[ 2][1]].ref_grad[2] = p3(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[ek[ 3][0]].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[ek[ 3][0]].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[ek[ 3][0]].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[ek[ 3][1]].ref_grad[0] = d1p3(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[ek[ 3][1]].ref_grad[1] = p3(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[ek[ 3][1]].ref_grad[2] = p3(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[ek[ 4][0]].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[ek[ 4][0]].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[ek[ 4][0]].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[ek[ 4][1]].ref_grad[0] = d1p0(point[0]) * p3(point[1]) * p0(point[2]);
          data.phi[ek[ 4][1]].ref_grad[1] = p0(point[0]) * d1p3(point[1]) * p0(point[2]);
          data.phi[ek[ 4][1]].ref_grad[2] = p0(point[0]) * p3(point[1]) * d1p0(point[2]);
          data.phi[ek[ 5][0]].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[ek[ 5][0]].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[ek[ 5][0]].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[ek[ 5][1]].ref_grad[0] = d1p1(point[0]) * p3(point[1]) * p0(point[2]);
          data.phi[ek[ 5][1]].ref_grad[1] = p1(point[0]) * d1p3(point[1]) * p0(point[2]);
          data.phi[ek[ 5][1]].ref_grad[2] = p1(point[0]) * p3(point[1]) * d1p0(point[2]);
          data.phi[ek[ 6][0]].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[ek[ 6][0]].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[ek[ 6][0]].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[ek[ 6][1]].ref_grad[0] = d1p0(point[0]) * p3(point[1]) * p1(point[2]);
          data.phi[ek[ 6][1]].ref_grad[1] = p0(point[0]) * d1p3(point[1]) * p1(point[2]);
          data.phi[ek[ 6][1]].ref_grad[2] = p0(point[0]) * p3(point[1]) * d1p1(point[2]);
          data.phi[ek[ 7][0]].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[ek[ 7][0]].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[ek[ 7][0]].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[ek[ 7][1]].ref_grad[0] = d1p1(point[0]) * p3(point[1]) * p1(point[2]);
          data.phi[ek[ 7][1]].ref_grad[1] = p1(point[0]) * d1p3(point[1]) * p1(point[2]);
          data.phi[ek[ 7][1]].ref_grad[2] = p1(point[0]) * p3(point[1]) * d1p1(point[2]);
          data.phi[ek[ 8][0]].ref_grad[0] = d1p0(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[ek[ 8][0]].ref_grad[1] = p0(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[ek[ 8][0]].ref_grad[2] = p0(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[ek[ 8][1]].ref_grad[0] = d1p0(point[0]) * p0(point[1]) * p3(point[2]);
          data.phi[ek[ 8][1]].ref_grad[1] = p0(point[0]) * d1p0(point[1]) * p3(point[2]);
          data.phi[ek[ 8][1]].ref_grad[2] = p0(point[0]) * p0(point[1]) * d1p3(point[2]);
          data.phi[ek[ 9][0]].ref_grad[0] = d1p1(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[ek[ 9][0]].ref_grad[1] = p1(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[ek[ 9][0]].ref_grad[2] = p1(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[ek[ 9][1]].ref_grad[0] = d1p1(point[0]) * p0(point[1]) * p3(point[2]);
          data.phi[ek[ 9][1]].ref_grad[1] = p1(point[0]) * d1p0(point[1]) * p3(point[2]);
          data.phi[ek[ 9][1]].ref_grad[2] = p1(point[0]) * p0(point[1]) * d1p3(point[2]);
          data.phi[ek[10][0]].ref_grad[0] = d1p0(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[ek[10][0]].ref_grad[1] = p0(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[ek[10][0]].ref_grad[2] = p0(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[ek[10][1]].ref_grad[0] = d1p0(point[0]) * p1(point[1]) * p3(point[2]);
          data.phi[ek[10][1]].ref_grad[1] = p0(point[0]) * d1p1(point[1]) * p3(point[2]);
          data.phi[ek[10][1]].ref_grad[2] = p0(point[0]) * p1(point[1]) * d1p3(point[2]);
          data.phi[ek[11][0]].ref_grad[0] = d1p1(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[ek[11][0]].ref_grad[1] = p1(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[ek[11][0]].ref_grad[2] = p1(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[ek[11][1]].ref_grad[0] = d1p1(point[0]) * p1(point[1]) * p3(point[2]);
          data.phi[ek[11][1]].ref_grad[1] = p1(point[0]) * d1p1(point[1]) * p3(point[2]);
          data.phi[ek[11][1]].ref_grad[2] = p1(point[0]) * p1(point[1]) * d1p3(point[2]);
          // face dofs
          data.phi[qk[0][0]].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[qk[0][0]].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[qk[0][0]].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[qk[0][1]].ref_grad[0] = d1p3(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[qk[0][1]].ref_grad[1] = p3(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[qk[0][1]].ref_grad[2] = p3(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[qk[0][2]].ref_grad[0] = d1p2(point[0]) * p3(point[1]) * p0(point[2]);
          data.phi[qk[0][2]].ref_grad[1] = p2(point[0]) * d1p3(point[1]) * p0(point[2]);
          data.phi[qk[0][2]].ref_grad[2] = p2(point[0]) * p3(point[1]) * d1p0(point[2]);
          data.phi[qk[0][3]].ref_grad[0] = d1p3(point[0]) * p3(point[1]) * p0(point[2]);
          data.phi[qk[0][3]].ref_grad[1] = p3(point[0]) * d1p3(point[1]) * p0(point[2]);
          data.phi[qk[0][3]].ref_grad[2] = p3(point[0]) * p3(point[1]) * d1p0(point[2]);
          data.phi[qk[1][0]].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[qk[1][0]].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[qk[1][0]].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[qk[1][1]].ref_grad[0] = d1p3(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[qk[1][1]].ref_grad[1] = p3(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[qk[1][1]].ref_grad[2] = p3(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[qk[1][2]].ref_grad[0] = d1p2(point[0]) * p3(point[1]) * p1(point[2]);
          data.phi[qk[1][2]].ref_grad[1] = p2(point[0]) * d1p3(point[1]) * p1(point[2]);
          data.phi[qk[1][2]].ref_grad[2] = p2(point[0]) * p3(point[1]) * d1p1(point[2]);
          data.phi[qk[1][3]].ref_grad[0] = d1p3(point[0]) * p3(point[1]) * p1(point[2]);
          data.phi[qk[1][3]].ref_grad[1] = p3(point[0]) * d1p3(point[1]) * p1(point[2]);
          data.phi[qk[1][3]].ref_grad[2] = p3(point[0]) * p3(point[1]) * d1p1(point[2]);
          data.phi[qk[2][0]].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[qk[2][0]].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[qk[2][0]].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[qk[2][1]].ref_grad[0] = d1p3(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[qk[2][1]].ref_grad[1] = p3(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[qk[2][1]].ref_grad[2] = p3(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[qk[2][2]].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p3(point[2]);
          data.phi[qk[2][2]].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p3(point[2]);
          data.phi[qk[2][2]].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p3(point[2]);
          data.phi[qk[2][3]].ref_grad[0] = d1p3(point[0]) * p0(point[1]) * p3(point[2]);
          data.phi[qk[2][3]].ref_grad[1] = p3(point[0]) * d1p0(point[1]) * p3(point[2]);
          data.phi[qk[2][3]].ref_grad[2] = p3(point[0]) * p0(point[1]) * d1p3(point[2]);
          data.phi[qk[3][0]].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[qk[3][0]].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[qk[3][0]].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[qk[3][1]].ref_grad[0] = d1p3(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[qk[3][1]].ref_grad[1] = p3(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[qk[3][1]].ref_grad[2] = p3(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[qk[3][2]].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p3(point[2]);
          data.phi[qk[3][2]].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p3(point[2]);
          data.phi[qk[3][2]].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p3(point[2]);
          data.phi[qk[3][3]].ref_grad[0] = d1p3(point[0]) * p1(point[1]) * p3(point[2]);
          data.phi[qk[3][3]].ref_grad[1] = p3(point[0]) * d1p1(point[1]) * p3(point[2]);
          data.phi[qk[3][3]].ref_grad[2] = p3(point[0]) * p1(point[1]) * d1p3(point[2]);
          data.phi[qk[4][0]].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[qk[4][0]].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[qk[4][0]].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[qk[4][1]].ref_grad[0] = d1p0(point[0]) * p3(point[1]) * p2(point[2]);
          data.phi[qk[4][1]].ref_grad[1] = p0(point[0]) * d1p3(point[1]) * p2(point[2]);
          data.phi[qk[4][1]].ref_grad[2] = p0(point[0]) * p3(point[1]) * d1p2(point[2]);
          data.phi[qk[4][2]].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p3(point[2]);
          data.phi[qk[4][2]].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p3(point[2]);
          data.phi[qk[4][2]].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p3(point[2]);
          data.phi[qk[4][3]].ref_grad[0] = d1p0(point[0]) * p3(point[1]) * p3(point[2]);
          data.phi[qk[4][3]].ref_grad[1] = p0(point[0]) * d1p3(point[1]) * p3(point[2]);
          data.phi[qk[4][3]].ref_grad[2] = p0(point[0]) * p3(point[1]) * d1p3(point[2]);
          data.phi[qk[5][0]].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[qk[5][0]].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[qk[5][0]].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[qk[5][1]].ref_grad[0] = d1p1(point[0]) * p3(point[1]) * p2(point[2]);
          data.phi[qk[5][1]].ref_grad[1] = p1(point[0]) * d1p3(point[1]) * p2(point[2]);
          data.phi[qk[5][1]].ref_grad[2] = p1(point[0]) * p3(point[1]) * d1p2(point[2]);
          data.phi[qk[5][2]].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p3(point[2]);
          data.phi[qk[5][2]].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p3(point[2]);
          data.phi[qk[5][2]].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p3(point[2]);
          data.phi[qk[5][3]].ref_grad[0] = d1p1(point[0]) * p3(point[1]) * p3(point[2]);
          data.phi[qk[5][3]].ref_grad[1] = p1(point[0]) * d1p3(point[1]) * p3(point[2]);
          data.phi[qk[5][3]].ref_grad[2] = p1(point[0]) * p3(point[1]) * d1p3(point[2]);
          // center dofs
          data.phi[56].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[56].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[56].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[57].ref_grad[0] = d1p3(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[57].ref_grad[1] = p3(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[57].ref_grad[2] = p3(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[58].ref_grad[0] = d1p2(point[0]) * p3(point[1]) * p2(point[2]);
          data.phi[58].ref_grad[1] = p2(point[0]) * d1p3(point[1]) * p2(point[2]);
          data.phi[58].ref_grad[2] = p2(point[0]) * p3(point[1]) * d1p2(point[2]);
          data.phi[59].ref_grad[0] = d1p3(point[0]) * p3(point[1]) * p2(point[2]);
          data.phi[59].ref_grad[1] = p3(point[0]) * d1p3(point[1]) * p2(point[2]);
          data.phi[59].ref_grad[2] = p3(point[0]) * p3(point[1]) * d1p2(point[2]);
          data.phi[60].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p3(point[2]);
          data.phi[60].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p3(point[2]);
          data.phi[60].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p3(point[2]);
          data.phi[61].ref_grad[0] = d1p3(point[0]) * p2(point[1]) * p3(point[2]);
          data.phi[61].ref_grad[1] = p3(point[0]) * d1p2(point[1]) * p3(point[2]);
          data.phi[61].ref_grad[2] = p3(point[0]) * p2(point[1]) * d1p3(point[2]);
          data.phi[62].ref_grad[0] = d1p2(point[0]) * p3(point[1]) * p3(point[2]);
          data.phi[62].ref_grad[1] = p2(point[0]) * d1p3(point[1]) * p3(point[2]);
          data.phi[62].ref_grad[2] = p2(point[0]) * p3(point[1]) * d1p3(point[2]);
          data.phi[63].ref_grad[0] = d1p3(point[0]) * p3(point[1]) * p3(point[2]);
          data.phi[63].ref_grad[1] = p3(point[0]) * d1p3(point[1]) * p3(point[2]);
          data.phi[63].ref_grad[2] = p3(point[0]) * p3(point[1]) * d1p3(point[2]);
        }
      }; // class Evaluator<...,Hypercube<3>>
    } // namespace Lagrange3
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_LAGRANGE3_EVALUATOR_HPP
