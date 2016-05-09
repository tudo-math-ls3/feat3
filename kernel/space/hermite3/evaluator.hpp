#pragma once
#ifndef KERNEL_SPACE_HERMITE3_EVALUATOR_HPP
#define KERNEL_SPACE_HERMITE3_EVALUATOR_HPP 1

// includes, FEAST
#include <kernel/space/parametric_evaluator.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace Hermite3
    {
      /**
       * \brief Hermite-3 Element Evaluator reference capabilities
       *
       * \author Peter Zajac
       */
      static constexpr SpaceTags ref_caps = SpaceTags::ref_value | SpaceTags::ref_grad | SpaceTags::ref_hess;

      /// \cond internal
      namespace Intern
      {
        // p1, p2, q1 and q2 are the 1D basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the evaluators.

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
      /// \endcond

      /**
       * \brief Hermite-3 Element Evaluator class template declaration.
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
       * \brief Hermite-3 Element Evaluator implementation for Hypercube<1> shape
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

        /// trafo evaluator type
        typedef TrafoEvaluator_ TrafoEvaluator;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

      protected:
        /// first-order derivative transform coefficient
        DataType _coeff;

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
          typedef typename TrafoEvaluator::template ConfigTraits<TrafoTags::jac_mat>::EvalDataType CoeffEvalData;
          CoeffEvalData coeff_data;
          trafo_eval(coeff_data, dom_point);

          // store coefficient
          _coeff = coeff_data.jac_mat(0,0);
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
          // Q1(x) = c/4 * (x + 1) * (x - 1)^2
          data.phi[1].ref_value = _coeff * Intern::q1(DataType(point[0]));
          // P2(x) = 1/2 - 1/4 * x * (x^2 - 3)
          data.phi[2].ref_value = Intern::p2(DataType(point[0]));
          // Q2(x) = c/4 * (x - 1) * (x + 1)^2
          data.phi[3].ref_value = _coeff * Intern::q2(DataType(point[0]));
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
          // dx Q1(x) = c/4 * (3*x + 1) * (x - 1)
          data.phi[1].ref_grad[0] = _coeff * Intern::dq1(DataType(point[0]));
          // dx P2(x) = -3/4 * (x^2 - 1)
          data.phi[2].ref_grad[0] = Intern::dp2(DataType(point[0]));
          // dx Q2(x) = c/4 * (3*x - 1) * (x + 1)
          data.phi[3].ref_grad[0] = _coeff * Intern::dq2(DataType(point[0]));
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
          // dxx P1(x) = +3/2 * x
          data.phi[0].ref_hess[0][0] = Intern::ddp1(DataType(point[0]));
          // dxx Q1(x) = 1/2 * (3*x - 1)
          data.phi[1].ref_hess[0][0] = _coeff * Intern::ddq1(DataType(point[0]));
          // dxx P2(x) = -3/2 * x
          data.phi[2].ref_hess[0][0] = Intern::ddp2(DataType(point[0]));
          // dxx Q2(x) = 1/2 * (3*x + 1)
          data.phi[3].ref_hess[0][0] = _coeff * Intern::ddq2(DataType(point[0]));
        }
      }; // class Evaluator<...,Hypercube<1>>

      /**
       * \brief Hermite-3 Element Evaluator implementation for Hypercube<2> shape
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

        /// trafo evaluator type
        typedef TrafoEvaluator_ TrafoEvaluator;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

      protected:
        /// first-order derivative transform coefficients
        Tiny::Matrix<DataType, 2, 2> _coeff_fod[4];

        /**
         * \brief Transforms the X-derivative basis function
         *
         * \param[in] i
         * The index of the vertex of the basis function.
         *
         * \param[in] rx, ry
         * The values of the X- and Y-derivative basis functions on the reference element.
         *
         * \returns
         * The value of the X-derivative basis function on the real element.
         */
        DataType _trans_dx(Index i, DataType rx, DataType ry) const
        {
          return _coeff_fod[i](0,0) * rx + _coeff_fod[i](0,1) * ry;
        }

        /**
         * \brief Transforms the Y-derivative basis function
         *
         * \param[in] i
         * The index of the vertex of the basis function.
         *
         * \param[in] rx, ry
         * The values of the X- and Y-derivative basis functions on the reference element.
         *
         * \returns
         * The value of the Y-derivative basis function on the real element.
         */
        DataType _trans_dy(Index i, DataType rx, DataType ry) const
        {
          return _coeff_fod[i](1,0) * rx + _coeff_fod[i](1,1) * ry;
        }

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
          typedef typename TrafoEvaluator::template ConfigTraits<TrafoTags::jac_mat>::EvalDataType CoeffEvalData;
          CoeffEvalData coeff_data;

          // loop over all four vertices of the quad
          for(int i(0); i < 4; ++i)
          {
            // compute reference vertex coordinates
            dom_point[0] = DataType(-1 + ((i << 1) & 2));
            dom_point[1] = DataType(-1 + ((i     ) & 2));

            // evaluate trafo data
            trafo_eval(coeff_data, dom_point);

            // store jacobian matrix
            _coeff_fod[i] = coeff_data.jac_mat;
          }
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
          const DataType x(point[0]), y(point[1]);
          DataType d0, d1;

          // Vertex 1: (-1,-1)
          data.phi[ 0].ref_value = Intern::p1(x) * Intern::p1(y);
          d0 = Intern::q1(x) * Intern::p1(y);
          d1 = Intern::p1(x) * Intern::q1(y);
          data.phi[ 1].ref_value = _trans_dx(0, d0, d1);
          data.phi[ 2].ref_value = _trans_dy(0, d0, d1);

          // Vertex 2: (+1,-1)
          data.phi[ 3].ref_value = Intern::p2(x) * Intern::p1(y);
          d0 = Intern::q2(x) * Intern::p1(y);
          d1 = Intern::p2(x) * Intern::q1(y);
          data.phi[ 4].ref_value = _trans_dx(1, d0, d1);
          data.phi[ 5].ref_value = _trans_dy(1, d0, d1);

          // Vertex 3: (-1,+1)
          data.phi[ 6].ref_value = Intern::p1(x) * Intern::p2(y);
          d0 = Intern::q1(x) * Intern::p2(y);
          d1 = Intern::p1(x) * Intern::q2(y);
          data.phi[ 7].ref_value = _trans_dx(2, d0, d1);
          data.phi[ 8].ref_value = _trans_dy(2, d0, d1);

          // Vertex 4: (+1,+1)
          data.phi[ 9].ref_value = Intern::p2(x) * Intern::p2(y);
          d0 = Intern::q2(x) * Intern::p2(y);
          d1 = Intern::p2(x) * Intern::q2(y);
          data.phi[10].ref_value = _trans_dx(3, d0, d1);
          data.phi[11].ref_value = _trans_dy(3, d0, d1);

          // inner basis functions
          data.phi[12].ref_value = Intern::q1(x) * Intern::q1(y);
          data.phi[13].ref_value = Intern::q2(x) * Intern::q1(y);
          data.phi[14].ref_value = Intern::q1(x) * Intern::q2(y);
          data.phi[15].ref_value = Intern::q2(x) * Intern::q2(y);
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
          const DataType x(point[0]), y(point[1]);
          DataType d0x, d0y, d1x, d1y;

          // Vertex 1: (-1,-1)
          data.phi[ 0].ref_grad[0] = Intern::dp1(x) * Intern::p1(y);
          data.phi[ 0].ref_grad[1] = Intern::p1(x) * Intern::dp1(y);
          d0x = Intern::dq1(x) * Intern::p1(y);
          d0y = Intern::q1(x) * Intern::dp1(y);
          d1x = Intern::dp1(x) * Intern::q1(y);
          d1y = Intern::p1(x) * Intern::dq1(y);
          data.phi[ 1].ref_grad[0] = _trans_dx(0, d0x, d1x);
          data.phi[ 1].ref_grad[1] = _trans_dx(0, d0y, d1y);
          data.phi[ 2].ref_grad[0] = _trans_dy(0, d0x, d1x);
          data.phi[ 2].ref_grad[1] = _trans_dy(0, d0y, d1y);

          // Vertex 2: (+1,-1)
          data.phi[ 3].ref_grad[0] = Intern::dp2(x) * Intern::p1(y);
          data.phi[ 3].ref_grad[1] = Intern::p2(x) * Intern::dp1(y);
          d0x = Intern::dq2(x) * Intern::p1(y);
          d0y = Intern::q2(x) * Intern::dp1(y);
          d1x = Intern::dp2(x) * Intern::q1(y);
          d1y = Intern::p2(x) * Intern::dq1(y);
          data.phi[ 4].ref_grad[0] = _trans_dx(1, d0x, d1x);
          data.phi[ 4].ref_grad[1] = _trans_dx(1, d0y, d1y);
          data.phi[ 5].ref_grad[0] = _trans_dy(1, d0x, d1x);
          data.phi[ 5].ref_grad[1] = _trans_dy(1, d0y, d1y);

          // Vertex 3: (-1,+1)
          data.phi[ 6].ref_grad[0] = Intern::dp1(x) * Intern::p2(y);
          data.phi[ 6].ref_grad[1] = Intern::p1(x) * Intern::dp2(y);
          d0x = Intern::dq1(x) * Intern::p2(y);
          d0y = Intern::q1(x) * Intern::dp2(y);
          d1x = Intern::dp1(x) * Intern::q2(y);
          d1y = Intern::p1(x) * Intern::dq2(y);
          data.phi[ 7].ref_grad[0] = _trans_dx(2, d0x, d1x);
          data.phi[ 7].ref_grad[1] = _trans_dx(2, d0y, d1y);
          data.phi[ 8].ref_grad[0] = _trans_dy(2, d0x, d1x);
          data.phi[ 8].ref_grad[1] = _trans_dy(2, d0y, d1y);

          // Vertex 4: (+1,+1)
          data.phi[ 9].ref_grad[0] = Intern::dp2(x) * Intern::p2(y);
          data.phi[ 9].ref_grad[1] = Intern::p2(x) * Intern::dp2(y);
          d0x = Intern::dq2(x) * Intern::p2(y);
          d0y = Intern::q2(x) * Intern::dp2(y);
          d1x = Intern::dp2(x) * Intern::q2(y);
          d1y = Intern::p2(x) * Intern::dq2(y);
          data.phi[10].ref_grad[0] = _trans_dx(3, d0x, d1x);
          data.phi[10].ref_grad[1] = _trans_dx(3, d0y, d1y);
          data.phi[11].ref_grad[0] = _trans_dy(3, d0x, d1x);
          data.phi[11].ref_grad[1] = _trans_dy(3, d0y, d1y);

          // inner basis functions
          data.phi[12].ref_grad[0] = Intern::dq1(x) * Intern::q1(y);
          data.phi[12].ref_grad[1] = Intern::q1(x) * Intern::dq1(y);
          data.phi[13].ref_grad[0] = Intern::dq2(x) * Intern::q1(y);
          data.phi[13].ref_grad[1] = Intern::q2(x) * Intern::dq1(y);
          data.phi[14].ref_grad[0] = Intern::dq1(x) * Intern::q2(y);
          data.phi[14].ref_grad[1] = Intern::q1(x) * Intern::dq2(y);
          data.phi[15].ref_grad[0] = Intern::dq2(x) * Intern::q2(y);
          data.phi[15].ref_grad[1] = Intern::q2(x) * Intern::dq2(y);
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
          const DataType x(point[0]), y(point[1]);
          DataType d0xx, d0xy, d0yy, d1xx, d1xy, d1yy;

          // Vertex 1: (-1,-1)
          data.phi[ 0].ref_hess[0][0] = Intern::ddp1(x) * Intern::p1(y);
          data.phi[ 0].ref_hess[0][1] =
          data.phi[ 0].ref_hess[1][0] = Intern::dp1(x) * Intern::dp1(y);
          data.phi[ 0].ref_hess[1][1] = Intern::p1(x) * Intern::ddp1(y);
          d0xx = Intern::ddq1(x) * Intern::p1(y);
          d0xy = Intern::dq1(x) * Intern::dp1(y);
          d0yy = Intern::q1(x) * Intern::ddp1(y);
          d1xx = Intern::ddp1(x) * Intern::q1(y);
          d1xy = Intern::dp1(x) * Intern::dq1(y);
          d1yy = Intern::p1(x) * Intern::ddq1(y);
          data.phi[ 1].ref_hess[0][0] = _trans_dx(0, d0xx, d1xx);
          data.phi[ 1].ref_hess[0][1] =
          data.phi[ 1].ref_hess[1][0] = _trans_dx(0, d0xy, d1xy);
          data.phi[ 1].ref_hess[1][1] = _trans_dx(0, d0yy, d1yy);
          data.phi[ 2].ref_hess[0][0] = _trans_dy(0, d0xx, d1xx);
          data.phi[ 2].ref_hess[0][1] =
          data.phi[ 2].ref_hess[1][0] = _trans_dy(0, d0xy, d1xy);
          data.phi[ 2].ref_hess[1][1] = _trans_dy(0, d0yy, d1yy);

          // Vertex 2: (+1,-1)
          data.phi[ 3].ref_hess[0][0] = Intern::ddp2(x) * Intern::p1(y);
          data.phi[ 3].ref_hess[0][1] =
          data.phi[ 3].ref_hess[1][0] = Intern::dp2(x) * Intern::dp1(y);
          data.phi[ 3].ref_hess[1][1] = Intern::p2(x) * Intern::ddp1(y);
          d0xx = Intern::ddq2(x) * Intern::p1(y);
          d0xy = Intern::dq2(x) * Intern::dp1(y);
          d0yy = Intern::q2(x) * Intern::ddp1(y);
          d1xx = Intern::ddp2(x) * Intern::q1(y);
          d1xy = Intern::dp2(x) * Intern::dq1(y);
          d1yy = Intern::p2(x) * Intern::ddq1(y);
          data.phi[ 4].ref_hess[0][0] = _trans_dx(1, d0xx, d1xx);
          data.phi[ 4].ref_hess[0][1] =
          data.phi[ 4].ref_hess[1][0] = _trans_dx(1, d0xy, d1xy);
          data.phi[ 4].ref_hess[1][1] = _trans_dx(1, d0yy, d1yy);
          data.phi[ 5].ref_hess[0][0] = _trans_dy(1, d0xx, d1xx);
          data.phi[ 5].ref_hess[0][1] =
          data.phi[ 5].ref_hess[1][0] = _trans_dy(1, d0xy, d1xy);
          data.phi[ 5].ref_hess[1][1] = _trans_dy(1, d0yy, d1yy);

          // Vertex 3: (-1,+1)
          data.phi[ 6].ref_hess[0][0] = Intern::ddp1(x) * Intern::p2(y);
          data.phi[ 6].ref_hess[0][1] =
          data.phi[ 6].ref_hess[1][0] = Intern::dp1(x) * Intern::dp2(y);
          data.phi[ 6].ref_hess[1][1] = Intern::p1(x) * Intern::ddp2(y);
          d0xx = Intern::ddq1(x) * Intern::p2(y);
          d0xy = Intern::dq1(x) * Intern::dp2(y);
          d0yy = Intern::q1(x) * Intern::ddp2(y);
          d1xx = Intern::ddp1(x) * Intern::q2(y);
          d1xy = Intern::dp1(x) * Intern::dq2(y);
          d1yy = Intern::p1(x) * Intern::ddq2(y);
          data.phi[ 7].ref_hess[0][0] = _trans_dx(2, d0xx, d1xx);
          data.phi[ 7].ref_hess[0][1] =
          data.phi[ 7].ref_hess[1][0] = _trans_dx(2, d0xy, d1xy);
          data.phi[ 7].ref_hess[1][1] = _trans_dx(2, d0yy, d1yy);
          data.phi[ 8].ref_hess[0][0] = _trans_dy(2, d0xx, d1xx);
          data.phi[ 8].ref_hess[0][1] =
          data.phi[ 8].ref_hess[1][0] = _trans_dy(2, d0xy, d1xy);
          data.phi[ 8].ref_hess[1][1] = _trans_dy(2, d0yy, d1yy);

          // Vertex 4: (+1,+1)
          data.phi[ 9].ref_hess[0][0] = Intern::ddp2(x) * Intern::p2(y);
          data.phi[ 9].ref_hess[0][1] =
          data.phi[ 9].ref_hess[1][0] = Intern::dp2(x) * Intern::dp2(y);
          data.phi[ 9].ref_hess[1][1] = Intern::p2(x) * Intern::ddp2(y);
          d0xx = Intern::ddq2(x) * Intern::p2(y);
          d0xy = Intern::dq2(x) * Intern::dp2(y);
          d0yy = Intern::q2(x) * Intern::ddp2(y);
          d1xx = Intern::ddp2(x) * Intern::q2(y);
          d1xy = Intern::dp2(x) * Intern::dq2(y);
          d1yy = Intern::p2(x) * Intern::ddq2(y);
          data.phi[10].ref_hess[0][0] = _trans_dx(3, d0xx, d1xx);
          data.phi[10].ref_hess[0][1] =
          data.phi[10].ref_hess[1][0] = _trans_dx(3, d0xy, d1xy);
          data.phi[10].ref_hess[1][1] = _trans_dx(3, d0yy, d1yy);
          data.phi[11].ref_hess[0][0] = _trans_dy(3, d0xx, d1xx);
          data.phi[11].ref_hess[0][1] =
          data.phi[11].ref_hess[1][0] = _trans_dy(3, d0xy, d1xy);
          data.phi[11].ref_hess[1][1] = _trans_dy(3, d0yy, d1yy);

          // inner basis functions
          data.phi[12].ref_hess[0][0] = Intern::ddq1(x) * Intern::q1(y);
          data.phi[12].ref_hess[0][1] =
          data.phi[12].ref_hess[1][0] = Intern::dq1(x) * Intern::dq1(y);
          data.phi[12].ref_hess[1][1] = Intern::q1(x) * Intern::ddq1(y);
          data.phi[13].ref_hess[0][0] = Intern::ddq2(x) * Intern::q1(y);
          data.phi[13].ref_hess[0][1] =
          data.phi[13].ref_hess[1][0] = Intern::dq2(x) * Intern::dq1(y);
          data.phi[13].ref_hess[1][1] = Intern::q2(x) * Intern::ddq1(y);
          data.phi[14].ref_hess[0][0] = Intern::ddq1(x) * Intern::q2(y);
          data.phi[14].ref_hess[0][1] =
          data.phi[14].ref_hess[1][0] = Intern::dq1(x) * Intern::dq2(y);
          data.phi[14].ref_hess[1][1] = Intern::q1(x) * Intern::ddq2(y);
          data.phi[15].ref_hess[0][0] = Intern::ddq2(x) * Intern::q2(y);
          data.phi[15].ref_hess[0][1] =
          data.phi[15].ref_hess[1][0] = Intern::dq2(x) * Intern::dq2(y);
          data.phi[15].ref_hess[1][1] = Intern::q2(x) * Intern::ddq2(y);
        }
      }; // class Evaluator<...,Hypercube<2>>

      /**
       * \brief Hermite-3 Element Evaluator implementation for Simplex<2> shape
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

        /// trafo evaluator type
        typedef TrafoEvaluator_ TrafoEvaluator;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

      protected:
        /// first-order derivative transform coefficients
        Tiny::Matrix<DataType, 2, 2> _coeff_fod;

        /**
         * \brief Transforms the X-derivative basis function
         *
         * \param[in] rx, ry
         * The values of the X- and Y-derivative basis functions on the reference element.
         *
         * \returns
         * The value of the X-derivative basis function on the real element.
         */
        DataType _trans_dx(DataType rx, DataType ry) const
        {
          return _coeff_fod(0,0) * rx + _coeff_fod(0,1) * ry;
        }

        /**
         * \brief Transforms the Y-derivative basis function
         *
         * \param[in] rx, ry
         * The values of the X- and Y-derivative basis functions on the reference element.
         *
         * \returns
         * The value of the Y-derivative basis function on the real element.
         */
        DataType _trans_dy(DataType rx, DataType ry) const
        {
          return _coeff_fod(1,0) * rx + _coeff_fod(1,1) * ry;
        }

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
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& trafo_eval)
        {
          // domain point: barycentre
          DomainPointType dom_point;
          dom_point[0] = dom_point[1] = DataType(1) / DataType(3);

          // evaluate trafo in interval midpoint
          typedef typename TrafoEvaluator::template ConfigTraits<TrafoTags::jac_mat>::EvalDataType CoeffEvalData;
          CoeffEvalData coeff_data;

          // evaluate trafo data
          trafo_eval(coeff_data, dom_point);

          // store jacobian matrix
          _coeff_fod = coeff_data.jac_mat;
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
          const DataType x(point[0]), y(point[1]);
          DataType d0, d1;

          static const DataType T1(DataType(1));
          static const DataType T2(DataType(2));
          static const DataType T3(DataType(3));
          static const DataType T7(DataType(7));
          static const DataType T13(DataType(13));
          static const DataType T27(DataType(27));

          // Vertex 1: (0,0)
          data.phi[ 0].ref_value = T13*x*y*(x + y - T1) + y*y*(T2*y - T3) + x*x*(T2*x - T3) + T1;
          d0 = x*(T1 + x*(x - T2) + y*(T3*(x - T1) + T2*y));
          d1 = y*(T1 + y*(y - T2) + x*(T3*(y - T1) + T2*x));
          data.phi[ 1].ref_value = _trans_dx(d0, d1);
          data.phi[ 2].ref_value = _trans_dy(d0, d1);

          // Vertex 2: (1,0)
          data.phi[ 3].ref_value = x*(x*(T3 - T2*x) + T7*y*(x + y - T1));
          d0 = x*(x*(x - T1) - T2*y*(x + y - T1));
          d1 = x*y*(T2*x + y - T1);
          data.phi[ 4].ref_value = _trans_dx(d0, d1);
          data.phi[ 5].ref_value = _trans_dy(d0, d1);

          // Vertex 3: (0,1)
          data.phi[ 6].ref_value = y*(y*(T3 - T2*y) + T7*x*(y + x - T1));
          d0 = y*x*(T2*y + x - T1);
          d1 = y*(y*(y - T1) - T2*x*(y + x - T1));
          data.phi[ 7].ref_value = _trans_dx(d0, d1);
          data.phi[ 8].ref_value = _trans_dy(d0, d1);

          // barycentre
          data.phi[ 9].ref_value = T27*x*y*(T1 - x - y);
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
          const DataType x(point[0]), y(point[1]);
          DataType d0x, d0y, d1x, d1y;

          static const DataType T1(DataType(1));
          static const DataType T2(DataType(2));
          static const DataType T3(DataType(3));
          static const DataType T4(DataType(4));
          static const DataType T6(DataType(6));
          static const DataType T7(DataType(7));
          static const DataType T13(DataType(13));
          static const DataType T27(DataType(27));

          // Vertex 0: (0, 0)
          data.phi[ 0].ref_grad[0] = T6*x*(x - T1) + T13*y*(T2*x + y - T1);
          data.phi[ 0].ref_grad[1] = T6*y*(y - T1) + T13*x*(T2*y + x - T1);
          d0x = x*(T3*x - T4) + y*(T6*x + T2*y - T3) + T1;
          d0y = x*(T4*y + T3*(x - T1));
          d1x = y*(T4*x + T3*(y - T1));
          d1y = y*(T3*y - T4) + x*(T6*y + T2*x - T3) + T1;
          data.phi[ 1].ref_grad[0] = _trans_dx(d0x, d1x);
          data.phi[ 1].ref_grad[1] = _trans_dx(d0y, d1y);
          data.phi[ 2].ref_grad[0] = _trans_dy(d0x, d1x);
          data.phi[ 2].ref_grad[1] = _trans_dy(d0y, d1y);

          // Vertex 1: (1, 0)
          data.phi[ 3].ref_grad[0] = T6*x*(T1 - x) + T7*y*(T2*x + y - T1);
          data.phi[ 3].ref_grad[1] = T7*x*(T2*y + x - T1);
          d0x = x*(T3*x - T2) + T2*y*(T1 - y - T2*x);
          d0y = T2*x*(T1 - x - T2*y);
          d1x = y*(T4*x + y - T1);
          d1y = x*(T2*(y + x) - T1);
          data.phi[ 4].ref_grad[0] = _trans_dx(d0x, d1x);
          data.phi[ 4].ref_grad[1] = _trans_dx(d0y, d1y);
          data.phi[ 5].ref_grad[0] = _trans_dy(d0x, d1x);
          data.phi[ 5].ref_grad[1] = _trans_dy(d0y, d1y);

          // Vertex 2: (0, 1)
          data.phi[ 6].ref_grad[0] = T7*y*(T2*x + y - T1);
          data.phi[ 6].ref_grad[1] = T6*y*(T1 - y) + T7*x*(T2*y + x - T1);
          d0x = y*(T2*(x + y) - T1);
          d0y = x*(T4*y + x - T1);
          d1x = T2*y*(T1 - y - T2*x);
          d1y = y*(T3*y - T2) + T2*x*(T1 - x - T2*y);
          data.phi[ 7].ref_grad[0] = _trans_dx(d0x, d1x);
          data.phi[ 7].ref_grad[1] = _trans_dx(d0y, d1y);
          data.phi[ 8].ref_grad[0] = _trans_dy(d0x, d1x);
          data.phi[ 8].ref_grad[1] = _trans_dy(d0y, d1y);

          // barycentre
          data.phi[ 9].ref_grad[0] = T27*y*(T1 - T2*x - y);
          data.phi[ 9].ref_grad[1] = T27*x*(T1 - T2*y - x);
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
          const DataType x(point[0]), y(point[1]);
          DataType d0xx, d0xy, d0yy, d1xx, d1xy, d1yy;

          static const DataType T1(DataType(1));
          static const DataType T2(DataType(2));
          static const DataType T3(DataType(3));
          static const DataType T4(DataType(4));
          static const DataType T6(DataType(6));
          static const DataType T7(DataType(7));
          static const DataType T13(DataType(13));
          static const DataType T14(DataType(14));
          static const DataType T26(DataType(26));
          static const DataType T27(DataType(27));
          static const DataType T54(DataType(54));

          // Vertex 1: (0,0)
          data.phi[ 0].ref_hess[0][0] = T26*y + T6*(T2*x - T1);
          data.phi[ 0].ref_hess[0][1] =
          data.phi[ 0].ref_hess[1][0] = T13*(T2*(x + y) - T1);
          data.phi[ 0].ref_hess[1][1] = T26*x + T6*(T2*y - T1);
          d0xx = T6*(x + y) - T4;
          d0xy = T4*y + T3*(T2*x - T1);
          d0yy = T4*x;
          d1xx = T4*y;
          d1xy = T4*x + T3*(T2*y - T1);
          d1yy = T6*(x + y) - T4;
          data.phi[ 1].ref_hess[0][0] = _trans_dx(d0xx, d1xx);
          data.phi[ 1].ref_hess[0][1] =
          data.phi[ 1].ref_hess[1][0] = _trans_dx(d0xy, d1xy);
          data.phi[ 1].ref_hess[1][1] = _trans_dx(d0yy, d1yy);
          data.phi[ 2].ref_hess[0][0] = _trans_dy(d0xx, d1xx);
          data.phi[ 2].ref_hess[0][1] =
          data.phi[ 2].ref_hess[1][0] = _trans_dy(d0xy, d1xy);
          data.phi[ 2].ref_hess[1][1] = _trans_dy(d0yy, d1yy);

          // Vertex 2: (1,0)
          data.phi[ 3].ref_hess[0][0] = T14*y + T6*(T1 - T2*x);
          data.phi[ 3].ref_hess[0][1] =
          data.phi[ 3].ref_hess[1][0] = T7*(T2*(x + y) - T1);
          data.phi[ 3].ref_hess[1][1] = T14*x;
          d0xx = T6*x - T2*(T2*y + T1);
          d0xy = T2 - T4*(x + y);
          d0yy = -T4*x;
          d1xx = T4*y;
          d1xy = T2*(T2*x + y) - T1;
          d1yy = T2*x;
          data.phi[ 4].ref_hess[0][0] = _trans_dx(d0xx, d1xx);
          data.phi[ 4].ref_hess[0][1] =
          data.phi[ 4].ref_hess[1][0] = _trans_dx(d0xy, d1xy);
          data.phi[ 4].ref_hess[1][1] = _trans_dx(d0yy, d1yy);
          data.phi[ 5].ref_hess[0][0] = _trans_dy(d0xx, d1xx);
          data.phi[ 5].ref_hess[0][1] =
          data.phi[ 5].ref_hess[1][0] = _trans_dy(d0xy, d1xy);
          data.phi[ 5].ref_hess[1][1] = _trans_dy(d0yy, d1yy);

          // Vertex 3: (0,1)
          data.phi[ 6].ref_hess[0][0] = T14*y;
          data.phi[ 6].ref_hess[0][1] =
          data.phi[ 6].ref_hess[1][0] = T7*(T2*(x + y) - T1);
          data.phi[ 6].ref_hess[1][1] = T14*x + T6*(T1 - T2*y);
          d0xx = T2*y;
          d0xy =  T2*(T2*y + x) - T1;
          d0yy = T4*x;
          d1xx = -T4*y;
          d1xy = T2 - T4*(y + x);
          d1yy = T6*y - T2*(T2*x + T1);
          data.phi[ 7].ref_hess[0][0] = _trans_dx(d0xx, d1xx);
          data.phi[ 7].ref_hess[0][1] =
          data.phi[ 7].ref_hess[1][0] = _trans_dx(d0xy, d1xy);
          data.phi[ 7].ref_hess[1][1] = _trans_dx(d0yy, d1yy);
          data.phi[ 8].ref_hess[0][0] = _trans_dy(d0xx, d1xx);
          data.phi[ 8].ref_hess[0][1] =
          data.phi[ 8].ref_hess[1][0] = _trans_dy(d0xy, d1xy);
          data.phi[ 8].ref_hess[1][1] = _trans_dy(d0yy, d1yy);

          // barycentre
          data.phi[ 9].ref_hess[0][0] = -T54*y;
          data.phi[ 9].ref_hess[0][1] =
          data.phi[ 9].ref_hess[1][0] = T27*(T1 - T2*(x + y));
          data.phi[ 9].ref_hess[1][1] = -T54*x;
        }
      }; // class Evaluator<...,Simplex<2>>
    } // namespace Hermite3
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_HERMITE3_EVALUATOR_HPP
