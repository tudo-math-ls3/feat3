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

        void prepare(const TrafoEvaluator_& trafo_eval)
        {
          // compute edge orientations
          Geometry::Intern::SubIndexMapping<Shape::Hypercube<2>, 1, 0> sim(
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,0>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<2,1>()[trafo_eval.get_cell_index()],
            trafo_eval.get_trafo().get_mesh().template get_index_set<1,0>());

          // fetch edge dof indices
          ek[0][0] =  4 + int(sim.map(0, 0));
          ek[0][1] =  5 - int(sim.map(0, 0));
          ek[1][0] =  6 + int(sim.map(1, 0));
          ek[1][1] =  7 - int(sim.map(1, 0));
          ek[2][0] =  8 + int(sim.map(2, 0));
          ek[2][1] =  9 - int(sim.map(2, 0));
          ek[3][0] = 10 + int(sim.map(3, 0));
          ek[3][1] = 11 - int(sim.map(3, 0));
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
    } // namespace Lagrange3
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_LAGRANGE3_EVALUATOR_HPP
