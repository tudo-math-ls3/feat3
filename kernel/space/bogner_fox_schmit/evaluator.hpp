#pragma once
#ifndef KERNEL_SPACE_BOGNER_FOX_SCHMIT_EVALUATOR_HPP
#define KERNEL_SPACE_BOGNER_FOX_SCHMIT_EVALUATOR_HPP 1

// includes, FEAST
#include <kernel/space/parametric_evaluator.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace BognerFoxSchmit
    {
      /**
       * \brief Bogner-Fox-Schmit Element Evaluator reference capabilities
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
          can_ref_grad = 1,
          /// can compute reference hessians
          can_ref_hess = 1,
        };
      };

      namespace Intern
      {
        template<typename T_>
        inline T_ p1(T_ x)
        {
          return T_(0.5) + T_(0.25) * x * (x*x - T_(3));
        }

        template<typename T_>
        inline T_ p2(T_ x)
        {
          return T_(0.5) - T_(0.25) * x * (x*x - T_(3));
        }

        template<typename T_>
        inline T_ q1(T_ x)
        {
          return T_(0.25) * (x + T_(1)) * Math::sqr(x - T_(1));
        }

        template<typename T_>
        inline T_ q2(T_ x)
        {
          return T_(0.25) * (x - T_(1)) * Math::sqr(x + T_(1));
        }

        // first order derivatives follow
        template<typename T_>
        inline T_ dp1(T_ x)
        {
          return +T_(0.75) * (x*x - T_(1));
        }

        template<typename T_>
        inline T_ dp2(T_ x)
        {
          return -T_(0.75) * (x*x - T_(1));
        }

        template<typename T_>
        inline T_ dq1(T_ x)
        {
          return T_(0.25) * (T_(3)*x + T_(1)) * (x - T_(1));
        }

        template<typename T_>
        inline T_ dq2(T_ x)
        {
          return T_(0.25) * (T_(3)*x - T_(1)) * (x + T_(1));
        }

        // second order derivatives follow
        template<typename T_>
        inline T_ ddp1(T_ x)
        {
          return +T_(1.5) * x;
        }

        template<typename T_>
        inline T_ ddp2(T_ x)
        {
          return -T_(1.5) * x;
        }

        template<typename T_>
        inline T_ ddq1(T_ x)
        {
          return T_(0.5) * (T_(3)*x - T_(1));
        }

        template<typename T_>
        inline T_ ddq2(T_ x)
        {
          return T_(0.5) * (T_(3)*x + T_(1));
        }
      } // namespace Intern

      /**
       * \brief Bogner-Fox-Schmit Element Evaluator class template declaration.
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
       * \brief Bogner-Fox-Schmit Element Evaluator implementation for Hypercube<1> shape
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
          ReferenceCapabilities>
      {
      public:
        /// base-class typedef
        typedef ParametricEvaluator<Evaluator, TrafoEvaluator_, SpaceEvalTraits_, ReferenceCapabilities> BaseClass;

        /// space type
        typedef Space_ SpaceType;

        /// space evaluation traits
        typedef SpaceEvalTraits_ SpaceEvalTraits;

        /// trafo evaluator type
        typedef TrafoEvaluator_ TrafoEvaluator;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

      protected:
        /// transform coefficient
        DataType _coeff;

        struct CoeffTrafoConfig :
          public Trafo::ConfigBase
        {
          enum
          {
            need_jac_det = 1
          };
        };

        typedef typename TrafoEvaluator::template ConfigTraits<CoeffTrafoConfig>::EvalDataType CoeffEvalData;

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
        Index get_num_local_dofs() const
        {
          return 4;
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& trafo_eval)
        {
          // domain point
          DomainPointType dom_point(DataType(0));

          // evaluate trafo in interval midpoint
          CoeffEvalData coeff_data;
          trafo_eval(coeff_data, dom_point);

          // store coefficient
          _coeff = coeff_data.jac_det;
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
        template<typename EvalData_>
        void eval_ref_values(
          EvalData_& data,
          const DomainPointType& point) const
        {
          // P1(x) = 1/2 + 1/4 * x * (x^2 - 3)
          data.phi[0].ref_value = Intern::p1(DataType(point[0]));
            //DataType(0.5) + DataType(0.25) * point[0] * (Math::sqr(point[0]) - DataType(3));
          // Q1(x) = c/4 * (x + 1) * (x - 1)^2
          data.phi[1].ref_value = _coeff * Intern::q1(DataType(point[0]));
            //DataType(0.25) * (point[0] + DataType(1)) * Math::sqr(point[0] - DataType(1));
          // P2(x) = 1/2 - 1/4 * x * (x^2 - 3)
          data.phi[2].ref_value = Intern::p2(DataType(point[0]));
            //DataType(0.5) - DataType(0.25) * point[0] * (Math::sqr(point[0]) - DataType(3));
          // Q2(x) = c/4 * (x - 1) * (x + 1)^2
          data.phi[3].ref_value = _coeff * Intern::q2(DataType(point[0]));
            //DataType(0.25) * (point[0] - DataType(1)) * Math::sqr(point[0] + DataType(1));
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
          // dx P1(x) =  +3/4 * (x^2 - 1)
          data.phi[0].ref_grad[0] = Intern::dp1(DataType(point[0]));
            //+DataType(0.75) * (Math::sqr(point[0]) - DataType(1));
          // dx Q1(x) = c/4 * (3*x + 1) * (x - 1)
          data.phi[1].ref_grad[0] = _coeff * Intern::dq1(DataType(point[0]));
            //DataType(0.25) * (DataType(3)*point[0] + DataType(1)) * (point[0] - DataType(1));
          // dx P2(x) = -3/4 * (x^2 - 1)
          data.phi[2].ref_grad[0] = Intern::dp2(DataType(point[0]));
            //-DataType(0.75) * (Math::sqr(point[0]) - DataType(1));
          // dx Q2(x) = c/4 * (3*x - 1) * (x + 1)
          data.phi[3].ref_grad[0] = _coeff * Intern::dq2(DataType(point[0]));
            //DataType(0.25) * (DataType(3)*point[0] - DataType(1)) * (point[0] + DataType(1));
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
          data.phi[0].ref_hess[0][0] = Intern::ddp1(DataType(point[0]));
          data.phi[1].ref_hess[0][0] = _coeff * Intern::ddq1(DataType(point[0]));
          data.phi[2].ref_hess[0][0] = Intern::ddp2(DataType(point[0]));
          data.phi[3].ref_hess[0][0] = _coeff * Intern::ddq2(DataType(point[0]));
        }
      }; // class Evaluator<...,Hypercube<1>>
    } // namespace BognerFoxSchmit
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_BOGNER_FOX_SCHMIT_EVALUATOR_HPP
